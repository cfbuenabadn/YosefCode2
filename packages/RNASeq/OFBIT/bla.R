rm(list=ls())
setwd("/data/yosef/users/allonwag/YosefCode/packages/RNASeq/summary/RCODE")
library(sva)
library(optparse)
library(preprocessCore)
library(gdata)

source("defaultCommandLineArguments.R")
option_list <- list(
  make_option("--collect", default=default_cmd_line_args[["--collect"]], type="character",
              help="Directory containing your RNA Seq results from the preproc pipeline."),
  make_option("--config", default=default_cmd_line_args[["--config"]], type="character",
              help="Config file for your project (.xls or .xlsx)."),
  make_option("--qcfields", default=default_cmd_line_args[["--qcfields"]], type="character"),
  make_option("--genefields", default=default_cmd_line_args[["--genefields"]] , type="character",
              help=""),
  make_option("--out", default=default_cmd_line_args[["--out"]], type="character",
              help=""),
  make_option("--lib", default=default_cmd_line_args[["--lib"]], type="character",
              help=""),
  make_option("--sigfile", default=default_cmd_line_args[["--sigfile"]], type="character",
              help=""),
  make_option("--housekeeping", default=default_cmd_line_args[["--housekeeping"]], type="character",
              help=""),
  make_option("--combat", action="store_true", default=default_cmd_line_args[["--combat"]],
              help="This will run the ComBat package for batch correction on your data."),
  make_option("--multiple_collect", type="character", default=default_cmd_line_args[["--multiple_collect"]],
              help="If you need to load multiple collect directories and config files, please supply a text file listing them here.."),
  make_option("--save_intermediate_files", action="store_true", default=default_cmd_line_args[["--save_intermediate_files"]],
              help="If the flag is set, more intermediate files are saved during the operation, which allows tracing and debugging of the script")
  
)

## ----- Parse Arguments -----
opt <- parse_args(OptionParser(option_list=option_list))
BC = opt$combat
collect_dir = opt$collect
config_file = opt$config 
qc_fields_file = opt$qcfields
gene_fields_file = opt$genefields
out_dir = opt$out
lib_dir = opt$lib
sig_file = opt$sigfile
# Point to all sig files
housekeeping_list = opt$housekeeping
save_intermediate_files = opt$save_intermediate_files


## ---- Source Summary Library
source(paste0(lib_dir,"/util.R"))
#source(paste0(lib_dir,"/loadRSEM.R"))
source(paste0(lib_dir,"/GeneFilter.R"))
source(paste0(lib_dir,"/SampleFilter.R"))
source(paste0(lib_dir,"/TechFilter.R"))
source(paste0(lib_dir,"/QuantileNormalization.R"))
source(paste0(lib_dir,"/TechCorrect.R"))


## ----- Load Expression Set -----
#collect_dir = "/data/yosef/TFH/processed/collect"
collect_dir = "/data/yosef/BRAIN/processed_June2015_b/collect"
load(file=file.path(collect_dir, "collectedRNASeqStudy.RData"))
#collect_dir = "/data/yosef/BRAIN/processed_July2015/collect"
#load(file=file.path(collect_dir, "collectedRNASeqStudy_withRSEM.RData"))
eSet = collectedRNASeqStudy$cuff_eSet

out_dir = file.path(collect_dir, "pipe_olfa3_HK1b")
housekeeping_list = "/home/eecs/allonwag/data/BRAIN/sources/metadata/olfactory_housekeeping_lists/OE_HKlist1b.txt"
## ----- Produce Output Directory
if (file.exists(out_dir)){
} else {
  dir.create(out_dir)
}


# olfactory experiments. Russell:
# Allon, since C1_run_id is one good way to track the experiments, to avoid any confusion, I've listed the C1 runs that should be run through the pipeline together for each experiment.  Alternatively the Expt_condition can be used, and I've listed those below the relevant C1_run_id.
# 
# Experiment 3 (WT vs. Sox cKO, UI and regeneration); so far across OEL01 through OEL12:  
#   C1_run_id: Y01, Y02, Y03, Y04, Y05, Y05a, Y05b, Y06, Y07, Y07a, Y07b, Y08, S01, S02, S03, S04, S04a, S04b, S05, S05a, S05b, S06.  
# Expt_condition:  K5ERRY_* and K5ERSox2cKO_*
#   
#   Experiment 4 (WT vs. p63 cKO, UI, time-points post KO); so far across OEL01 to OEL12:
#   C1_run_id: Y01, Y04, P01, P02, P03, P03a, P03b, P04, P05.
# Expt_condition: K5ERRY_UI_* and K5ERp63cKO_*


EXPERIMENT = "olfa_experiment3"
if(EXPERIMENT == "cortical") {
  
} else if(EXPERIMENT == "olfa_experiment3" || EXPERIMENT == "olfa_experiment4") {
  #out_dir = file.path(collect_dir, EXPERIMENT)  
  if (!file.exists(out_dir)){
    dir.create(out_dir)
  }
  
  if(EXPERIMENT == "olfa_experiment3") {
    desiredSamples = startsWith(phenoData(eSet)$MD_expt_condition, "K5ERRY_", ignore.case=T) | startsWith(phenoData(eSet)$MD_expt_condition, "K5ERSox2cKO_", ignore.case=T)
  } else if(EXPERIMENT == "olfa_experiment4") {
    desiredSamples = startsWith(phenoData(eSet)$MD_expt_condition, "K5ERRY_UI_", ignore.case=T) | startsWith(phenoData(eSet)$MD_expt_condition, "K5ERp63cKO_", ignore.case=T)
  }else {
    stop("code should never be reached")
  }
  
  eSet = eSet[, desiredSamples]  
} else {
  stop("unrecognized experiment")
}


## ----- Pre-Filtering of Failed Samples -----
# Remove all samples failing the preprocessing step, based on config file. 
# Note: These should have all NA/0 TPM from collect OR all NaNs as their quality metrics

is.failed = #grepl("Failure",phenoData(eSet)$Preproc_Code) |
  apply(is.na(exprs(eSet)) | (exprs(eSet) == 0), 2, all) |
  apply(is.na(pData(protocolData(eSet))), 1, all)
print(paste(sum(is.failed),"samples failed pipeline:"),quote = F)
write.table(matrix(colnames(eSet)[is.failed],ncol = 1),file = paste0(out_dir,"/failed_preproc_list.txt"),row.names=F, col.names=F,quote = F)

prefilt.eSet = eSet[,!is.failed]
print(sprintf("Removed %d failed samples.", sum(is.failed)))

## ----- Pre-Filtering of Transcripts: Coding + Detected ----
# Select only type 1 transcripts (Coding)
type1.eSet = prefilt.eSet[featureData(eSet)$Transcript_Type == "protein_coding",]
# Select only detected transcripts
is.expressed.sc = rowMeans(exprs(type1.eSet)) > 0

sc.eSet = type1.eSet[which(is.expressed.sc),]
print("Removed non-Type1 transcripts and undetected transcripts.")

## ----- Checkpoint: No NA shall pass! -----
stopifnot(!any(is.na(exprs(sc.eSet))))

## ----- Preliminary Gene Filtering
init.gf.vec = GeneFilter(sc.eSet,
                         count.cutoff = 10,
                         prop.failed = .80,
                         verbose = T,
                         plot.dir = paste0(out_dir,"/genefilter"),
                         plot.prefix = "pre_cell_filtering")

## ----- Create Sample Filtering Directories
if (!file.exists(paste0(out_dir,"/samplefilter/"))){dir.create(paste0(out_dir,"/samplefilter/"))}
if (!file.exists(paste0(out_dir,"/samplefilter_tech/"))){dir.create(paste0(out_dir,"/samplefilter_tech/"))}

## ----- Hard-Cutoff (Cell) Sample Filtering + Correlation with Technical Features
# The Hard-Cutoff module is for diagnostic purposes only - output is not used.
hard.sf.sc.eSet = SampleFilter(eSet = sc.eSet,gene.filter.vec = init.gf.vec, housekeeping_list = housekeeping_list,mixture = T,verbose = T, MIN_RALIGN = 5, Z_CUTOFF=NULL, plot.dir = paste0(out_dir,"/samplefilter/hard_cutoff") )
hard.gf.vec = GeneFilter(hard.sf.sc.eSet,
                         count.cutoff = 10,
                         prop.failed = .80,
                         verbose = T,
                         plot.dir = paste0(out_dir,"/genefilter"),
                         plot.prefix = "post_hard_cell_filtering")
hard.tf.sc.eSet = TechFilter(hard.sf.sc.eSet,gf.vec = hard.gf.vec,Z_CUTOFF = NULL,PROP_CUTOFF = .9,plot.dir = paste0(out_dir,"/samplefilter_tech/hard_cutoff"),plot.prefix = "")


TURN_ON_ADAPTIVE_FILTERING = F #I found that the adaptive filtering is often too strict
if(TURN_ON_ADAPTIVE_FILTERING) {
  ## ----- Adaptive Sample and Gene Filtering
  sf.sc.eSet = SampleFilter(eSet = sc.eSet,gene.filter.vec = init.gf.vec, housekeeping_list = housekeeping_list,mixture = T,verbose = T,MIN_RALIGN = 5, plot.dir = paste0(out_dir,"/samplefilter/adaptive") )
  # Update gf.vec
  gf.vec = GeneFilter(sf.sc.eSet,
                      count.cutoff = 10,
                      prop.failed = .80,
                      verbose = T,
                      plot.dir = paste0(out_dir,"/genefilter"),
                      plot.prefix = "post_adapt_cell_filtering_1")
  tf.sc.eSet = TechFilter(sf.sc.eSet,gf.vec = gf.vec,PROP_CUTOFF = .9,plot.dir = paste0(out_dir,"/samplefilter_tech/adaptive"),plot.prefix = "")
  
}else {
  #skip adaptive filtering and just use the hard filtering
  sf.sc.eSet = hard.sf.sc.eSet
  gf.vec = hard.gf.vec
  tf.sc.eSet = hard.tf.sc.eSet
}


# Update gf.vec
gf.vec = GeneFilter(sf.sc.eSet,
                    count.cutoff = 10,
                    prop.failed = .80,
                    verbose = T,
                    plot.dir = paste0(out_dir,"/genefilter"),"post_adapt_cell_filtering_2")

if(save_intermediate_files) 
{
  write.table(exprs(tf.sc.eSet), file=paste0(out_dir,"/exprsAfterTechFilter.txt"), sep = "\t", col.names = NA, quote=F)
  save(sf.sc.eSet, gf.vec, file=file.path(out_dir, "imageAfterFiltering.RData"))
}

##----- Normalization
if (!file.exists(paste0(out_dir,"/normalization/"))){dir.create(paste0(out_dir,"/normalization/"))}

nrm.sc.matrix = QuantileNormalization(exprs(tf.sc.eSet), gf.vec = gf.vec, plot.dir = paste0(out_dir,"/normalization/qn_sc"))
print("Working with quantile normalization!")

#nrm.sc.matrix = QuartileNormalization(exprs(tf.sc.eSet), gf.vec = gf.vec, plot.dir = paste0(out_dir,"/normalization/qrn_sc"))
#print("Working with upper quartile normalization!")

nrm.sc.eSet = tf.sc.eSet
exprs(nrm.sc.eSet) = nrm.sc.matrix
gf.vec = gf.vec & (apply(exprs(nrm.sc.eSet),1,sd) > 10^(-6))

if(save_intermediate_files) 
{
  write.table(exprs(nrm.sc.eSet), file=paste0(out_dir,"/exprsAfterScaling.txt"), sep = "\t", col.names = NA, quote=F)
}

## Now do RUV
source(paste0(lib_dir,"/runRUVg.R"))
RUV_K = 1
hk_genes =  as.character(as.matrix(read.table(housekeeping_list)))
# RUVg.nrm.sc.matrix = runRUVg(nrm.sc.matrix, rownames(nrm.sc.matrix) %in% hk_genes, K = RUV_K);
# 
# if(save_intermediate_files) 
# {
#   write.table(RUVg.nrm.sc.matrix, file=paste0(out_dir,"/exprsAfterRUVg.txt"), sep = "\t", col.names = NA, quote=F)
# }
# 
# RUVg.nrm.sc.matrix.only.filtered.genes = runRUVg(nrm.sc.matrix[gf.vec,], rownames(nrm.sc.matrix[gf.vec,]) %in% hk_genes, K = RUV_K);
# 
# if(save_intermediate_files) 
# {
#   write.table(RUVg.nrm.sc.matrix.only.filtered.genes, file=paste0(out_dir,"/exprsAfterRUVg_onlyFilteredGenes.txt"), sep = "\t", col.names = NA, quote=F)
# }
# 
RUVg.nrm.sc.matrix.only.hk.filtered.genes = runRUVg(nrm.sc.matrix, gf.vec & (rownames(nrm.sc.matrix) %in% hk_genes), K = RUV_K);

if(save_intermediate_files) 
{
  write.table(RUVg.nrm.sc.matrix.only.hk.filtered.genes, file=paste0(out_dir,"/exprsAfterRUVg_onlyHKFilteredGenes.txt"), sep = "\t", col.names = NA, quote=F)
  save(RUVg.nrm.sc.matrix.only.hk.filtered.genes, gf.vec, file=file.path(out_dir, "imageAfterNormalization.RData"))
  
}

# # as expected, RUVg computes the residuals for each gene separately (after deducing the W from the HK genes)
# #so as long as the HK genes have only filtered genes it doesn't matter if the extrapolation is done on all the genes or only on the filtered subset
# max(abs(RUVg.nrm.sc.matrix.only.hk.filtered.genes["Smap1",] - RUVg.nrm.sc.matrix.only.filtered.genes["Smap1",]))

# #as.character(unlist((read.table(housekeeping_list)))) 
# 
# ##----- Put COMBAT here
# ##----- Technical Adjustment
# 
# tc.sc.matrix = TechCorrect(nrm.sc.eSet,ignore.zeroes = F,QTHRESH = .01, gf.vec = gf.vec,PROP_CUTOFF = .9, plot.dir = paste0(out_dir,"/normalization/tech"))
# tc.sc.eSet = sf.sc.eSet
# exprs(tc.sc.eSet) = tc.sc.matrix
# if(save_intermediate_files) 
# {
#   write.table(tc.sc.matrix, file=paste0(out_dir,"/exprsAfterTechCorrect.txt"), sep = "\t", col.names = NA, quote=F)
# }
# 
# ##----- Projections: Weighted PCA
# 
# # Weights
# fnr_weights = FNRw(tc.sc.eSet,tc.sc.eSet,gf.vec = gf.vec, FN_thresh = 0,housekeeping_list = housekeeping_list)
# 
# ###### Combat Batch Correction #######
# # Maybe wrap up the signature analysis into a module, and then you can do before and after.
# # This is still in development, please check to make sure the parameters make sense for your study.
# if(BC){
#   
#   dfPheno <- pData(tc.sc.eSet)
#   #vCompare <- c("Condition_Code","Time_Code")
#   
#   exprs(tc.sc.eSet) <- log10(exprs(tc.sc.eSet)+1)
#   combat.eSet = tc.sc.eSet[which(gf.vec),]
#   modcombat = model.matrix(~ as.factor(Condition_Code), data=dfPheno)
#   combat_edata = ComBat(dat=exprs(combat.eSet), batch=pData(combat.eSet)$Batch_Code, mod=modcombat, par.prior=T, prior.plots=T)
#   
#   dim(combat.eSet)
#   dim(combat_edata)
#   
#   vAdjGenes <- rownames(combat_edata)
#   combat.eSet <- tc.sc.eSet
#   
#   if (all(rownames(exprs(combat.eSet)[vAdjGenes,])==rownames(combat_edata[vAdjGenes,]))){
#     exprs(combat.eSet)[vAdjGenes,] <- combat_edata[vAdjGenes,]
#   }
#   tc.sc.eSet <- combat.eSet
# }
# 
# # wPCA on filtered genes (to save time)
# source(paste0(lib_dir,"/wPCA.R"))
# x = log10(exprs(tc.sc.eSet)[gf.vec,]+1)
# w = fnr_weights[gf.vec,]
# w[w != 1] = 1 - w[w != 1]
# wpc = wPCA(x,w)
# 
# conditions = phenoData(tc.sc.eSet)$Condition_Code
# batch = phenoData(tc.sc.eSet)$Batch_Code
# 
# pdf(paste0(out_dir,"/condvswpc.pdf"))
# condition_colors = rainbow(length(levels(conditions)))
# plot(wpc$x[,1],wpc$x[,2],col = condition_colors[conditions], xlab = "wPC1", ylab = "wPC2" )
# legend(x = "topleft",legend = levels(conditions),pch = 1,col = condition_colors )
# dev.off()
# 
# pdf(paste0(out_dir,"/batchvswpc.pdf"))
# batch_colors = rainbow(length(levels(batch)))
# plot(wpc$x[,1],wpc$x[,2],col = batch_colors[batch], xlab = "wPC1", ylab = "wPC2" )
# legend(x = "topleft",legend = levels(batch),pch = 1,col = batch_colors )
# dev.off()
# 
# ### STOP HERE FOR NEXT WEEK
# 
# ## ----- Signature Analysis
# 
# source(paste0(lib_dir,"/CalcSig.R"))
# 
# x = log10(exprs(tc.sc.eSet)[gf.vec,]+1)
# w = fnr_weights[gf.vec,]
# s = CalcSig(x,w,table = sig_file,scale = T)
# Q = na.omit(t(processQf(pData(protocolData(tc.sc.eSet)),rownames(protocolData(tc.sc.eSet)))))
# 
# qs = rbind(Q,s)
# qs = qs[apply(qs,1,sd) > 0,]
# 
# MAX_EXP_PCS = 5
# cors = cor(t(qs),wpc$x,method = 'spearman')
# Fish = (1/2)*log((1+cors)/(1-cors))
# z = sqrt((dim(x)[2]-3)/1.06)*Fish
# p = 2*pnorm(-abs(z))
# is.top = 1:dim(p)[2] %in% 1:MAX_EXP_PCS
# p = p[,is.top]
# q = matrix(p.adjust(unlist(p),method = 'fdr'),nrow = dim(p)[1])
# is.sig.assoc.pc = colSums(q < .01) > 0
# is.sig.assoc.sig = rowSums(q < .01) > 0
# 
# library(gplots)
# colnames(cors) = paste("wPC",1:dim(cors)[2])
# pdf(paste0(out_dir,"/sigsvwpc.pdf"),width = 12)
# svpc = heatmap.2(cors[is.sig.assoc.sig,][,is.top][,is.sig.assoc.pc],margins = c(5,40),key.title = "",key.xlab  = "Spearman Correlation",density.info = 'none',trace = 'none',col = colorRampPalette(c("purple","black","yellow")),cexRow = .5,cexCol = .5)
# dev.off()
# 
# sig_table = cbind(rownames(cors),q)[order(apply(q,1,min)),]
# write.table(sig_table,quote = F,row.names = F, col.names = F,paste0(out_dir,"/sigsvwpc.txt"))
# 
# ## Arrow Plots
# NUM_LEVELS = 2
# 
# pdf(paste0(out_dir,"/arrowsigsvwpc.pdf"))
# plot(wpc$x[,1]/max(abs((wpc$x[,1]))),wpc$x[,2]/max(abs((wpc$x[,2]))), xlim = c(-1,1),ylim = c(-1,1),pch = 1, col = condition_colors[conditions], xlab = "wPC1", ylab = "wPC2")
# legend(x = "topleft",legend = levels(conditions),pch = 1,col = condition_colors )
# 
# # Sigs Associated with PC1
# for (i in 1:NUM_LEVELS){
#   sname = rownames(cors)[order(q[,1])][i]
#   lm1 = lm(qs[sname,] ~ wpc$x[,1] + wpc$x[,2] )
#   vec = c(lm1$coefficients[2],lm1$coefficients[3])
#   nvec = vec/sqrt(sum(vec^2))
#   arrows(0,0,summary(lm1)$r.squared*nvec[1],summary(lm1)$r.squared*nvec[2], lwd = 2, lty = i, col = 'black')
# }
# 
# # Sigs Associated with PC2
# for (i in 1:NUM_LEVELS){
#   sname = rownames(cors)[order(q[,2])][i]
#   lm1 = lm(qs[sname,] ~ wpc$x[,1] + wpc$x[,2] )
#   vec = c(lm1$coefficients[2],lm1$coefficients[3])
#   nvec = vec/sqrt(sum(vec^2))
#   arrows(0,0,summary(lm1)$r.squared*nvec[1],summary(lm1)$r.squared*nvec[2], lwd = 2, lty = i, col = 'purple')
# }
# 
# all.snames = c(rownames(cors)[order(q[,1])][1:NUM_LEVELS],rownames(cors)[order(q[,2])][1:NUM_LEVELS])
# all.cols = c(rep("black",NUM_LEVELS),rep("purple",NUM_LEVELS))
# all.types = rep(c(1:NUM_LEVELS),2)
# legend(x = "bottomleft",legend = all.snames, lty = all.types, col = all.cols ,cex = .5)
# dev.off()
# 
# ## Yosef Plots
# source(paste0(lib_dir,"/YosefPlot.R"))
# x = exprs(tc.sc.eSet)
# YosefPlot(x,gf.vec = gf.vec,as.character(unlist(read.table(housekeeping_list))),draw.lines = T,out.dir = paste0(out_dir,"/yplots"),plot.name = "hk_genes.pdf")
# 
# for (sname in c(rownames(cors)[q[,1] == min(q[,1])],rownames(cors)[q[,1] == sort(q[,1])[2]],rownames(cors)[q[,2] == min(q[,2])],rownames(cors)[q[,2] == sort(q[,2])[2]])){
#   if (!sname %in% rownames(Q)){
#     s1 = GetSigSet(table = sig_file,name = sname)
#     YosefPlot(x,gf.vec = gf.vec, s1[,1],draw.lines = T,out.dir = paste0(out_dir,"/yplots"),plot.name = paste0(sname,".pdf"))
#   }
# }
# 
