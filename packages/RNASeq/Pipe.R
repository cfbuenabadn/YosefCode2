rm(list=ls())
library(sva)
library(optparse)
library(CCA)
library(ppls)
option_list <- list(
  make_option("--collect", default="/data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates1_d7_Arm/rsem", type="character",
              help="Directory containing your RNA Seq results from the preproc pipeline."),
  make_option("--config", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/TestConfig-All.xls", type="character",
              help="Config file for your project."),
  make_option("--qcfields", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/qc_fields.txt", type="character"),
  make_option("--genefields", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/gene_fields.txt", type="character",
              help=""),
  make_option("--out", default="/data/yosef/CD8_effector_diff/out/SingleCell-RNA-Seq/04-14-2015_NoBC_40bins_MedianVersionOfTechCorrect1_SLOWLY2", type="character",
              help=""),
  make_option("--lib", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline", type="character",
              help=""),
  make_option("--sigfile", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/immsig-ProperCase.txt", type="character",
              help=""),
  make_option("--housekeeping", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/house_keeping_mouse_TitleCase.txt", type="character",
              help=""),
  make_option("--combat", action="store_true", default=FALSE,
              help="This will run the ComBat package for batch correction on your data."),
  make_option("--multiple_collect", type="character", default="/data/yosef/CD8_effector_diff/src/SummaryPipeline/CollectFolders.txt",
              help="If you need to load multiple collect directories and config files, please supply a text file listing them here..")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

BC = opt$combat
collect_dir = opt$collect
config_file = opt$config 
qc_fields_file = opt$qcfields
gene_fields_file = opt$genefields
out_dir = opt$out
lib_dir = opt$lib
sig_file = opt$sigfile

fileSample<-file(paste0(out_dir,"/SampleLog.tab"))
write(("Sample Information"), fileSample)


# Point to all sig files
housekeeping_list = opt$housekeeping

## ----- Produce Output Directory
if (file.exists(out_dir)){
  stop("Output directory already exists! rm -r and try again...")
} else {
  dir.create(out_dir)
}

source(paste0(lib_dir,"/loadRSEM.R"))
## ----- LoadRSEMStudy - function to load multiple collect directories
if(opt$multiple_collect != ""){
  dfCollect <- read.table(opt$multiple_collect,header=T,sep="\t")
}
apply(dfCollect,1,function(x) print(x[2]))
li_eSet <-apply(dfCollect,1,function(x) loadRSEM(collect_dir = x[2],config_file =x[1],qc_fields_file = qc_fields_file,gene_fields_file = gene_fields_file))

eSet = li_eSet[[1]]
if (length(li_eSet) > 1){
  for (i in 2:length(li_eSet)){
    eSet <- Biobase::combine(eSet,li_eSet[[i]])
  }
}

summary(eSet)


## ----- Load Data and Pre-Filtering of Failed Samples
source(paste0(lib_dir,"/loadRSEM.R"))
eSet = loadRSEM(collect_dir = collect_dir,config_file = config_file,qc_fields_file = qc_fields_file,gene_fields_file = gene_fields_file)

# Remove failed samples: is.na(TPM)
is.failed = grepl("Failure",phenoData(eSet)$Preproc_Code)
print(paste(sum(is.failed),"samples failed pipeline:"),quote = F)
write.table(matrix(colnames(eSet)[is.failed],ncol = 1),row.names=F, col.names=F,quote = F)
print("Removed failed samples.")
prefilt.eSet = eSet[,!is.failed]

## ----- Filter Reads
# Select only type 1 transcripts (Coding)
type1.eSet = prefilt.eSet[featureData(eSet)$Transcript_Type == 1,]

# Select only detected transcripts
is.expressed.sc = rowMeans(exprs(type1.eSet)) > 0
sc.eSet = type1.eSet[which(is.expressed.sc),]

# No nan shall pass
stopifnot(!any(is.na(exprs(sc.eSet))))

source(paste0(lib_dir,"/GeneFilter.R"))
gf.vec = GeneFilter(sc.eSet,
                    count.cutoff = 10,
                    prop.failed = .85,
                    verbose = T,
                    plot.dir = paste0(out_dir,"/genefilter"))

## ----- Sample Filtering
source(paste0(lib_dir,"/SampleFilter.R"))
sf.sc.eSet = SampleFilter(eSet = sc.eSet,gene.filter.vec = gf.vec, housekeeping_list = housekeeping_list,mixture = T,verbose = T,plot.dir = paste0(out_dir,"/samplefilter_sc") )
gf.vec = gf.vec & (apply(exprs(sf.sc.eSet),1,sd) > 0)


##----- Normalization
source(paste0(lib_dir,"/QuantileNormalization.R"))
qn.sc.matrix = QuantileNormalization(exprs(sf.sc.eSet), gf.vec = gf.vec, plot.dir = paste0(out_dir,"/qn_sc"))
qn.sc.eSet = sf.sc.eSet
exprs(qn.sc.eSet) = qn.sc.matrix
gf.vec = gf.vec & (apply(exprs(qn.sc.eSet),1,sd) > 10^(-10))

##----- Technical Adjustment

source(paste0(lib_dir,"/TechCorrect.R"))
tc.sc.matrix = TechCorrect(qn.sc.eSet,ignore.zeroes = F, gf.vec = gf.vec, maxnumbins = 10,Z_CUTOFF = .5,PROP_CUTOFF = .9,plot.dir = paste0(out_dir,"/tech_sc_10"))
tc.sc.eSet = sf.sc.eSet
exprs(tc.sc.eSet) = tc.sc.matrix


##----- Projections: Weighted PCA

# Weights
fnr_weights = FNRw(tc.sc.eSet,tc.sc.eSet,gf.vec = gf.vec, FN_thresh = 0,housekeeping_list = housekeeping_list)

# wPCA on filtered genes (to save time)
source(paste0(lib_dir,"/wPCA.R"))
x = log10(exprs(tc.sc.eSet)[gf.vec,]+1)
w = fnr_weights[gf.vec,]
w[w != 1] = 1 - w[w != 1]
wpc = wPCA(x,w)

conditions = phenoData(tc.sc.eSet)$Condition_Code
batch = phenoData(tc.sc.eSet)$Batch_Code

pdf(paste0(out_dir,"/condvswpc.pdf"))
condition_colors = rainbow(length(levels(conditions)))
plot(wpc$x[,1],wpc$x[,2],col = condition_colors[conditions], xlab = "wPC1", ylab = "wPC2" )
legend(x = "topleft",legend = levels(conditions),pch = 1,col = condition_colors )
dev.off()

pdf(paste0(out_dir,"/batchvswpc.pdf"))
batch_colors = rainbow(length(levels(batch)))
plot(wpc$x[,1],wpc$x[,2],col = batch_colors[batch], xlab = "wPC1", ylab = "wPC2" )
legend(x = "topleft",legend = levels(batch),pch = 1,col = batch_colors )
dev.off()

## ----- Signature Analysis

source(paste0(lib_dir,"/CalcSig.R"))

x = log10(exprs(tc.sc.eSet)[gf.vec,]+1)
w = fnr_weights[gf.vec,]
s = CalcSig(x,w,table = sig_file,scale = T)
Q = na.omit(t(processQf(pData(protocolData(tc.sc.eSet)),rownames(protocolData(tc.sc.eSet)))))

qs = rbind(Q,s)
qs = qs[apply(qs,1,sd) > 0,]

cors = cor(t(qs),wpc$x,method = 'spearman')
Fish = (1/2)*log((1+cors)/(1-cors))
z = sqrt((dim(x)[2]-3)/1.06)*Fish
p = 2*pnorm(-abs(z))
is.top = wpc$sdev^2 >= sum(wpc$sdev^2)/(dim(x)[2]-1)
p = p[,is.top]
q = matrix(p.adjust(unlist(p),method = 'fdr'),nrow = dim(p)[1])
is.sig.assoc.pc = colSums(q < .01) > 0
is.sig.assoc.sig = rowSums(q < .01) > 0

library(gplots)
colnames(cors) = paste("wPC",1:dim(cors)[2])
pdf(paste0(out_dir,"/sigsvwpc.pdf"),width = 12)
svpc = heatmap.2(cors[is.sig.assoc.sig,][,is.top][,is.sig.assoc.pc],margins = c(5,40),key.title = "",key.xlab  = "Spearman Correlation",density.info = 'none',trace = 'none',col = colorRampPalette(c("purple","black","yellow")),cexRow = .5,cexCol = .5)
dev.off()

sig_table = cbind(rownames(cors),q)[order(apply(q,1,min)),]
write.table(sig_table,quote = F,row.names = F, col.names = F,paste0(out_dir,"/sigsvwpc.txt"))

## Arrow Plots
NUM_LEVELS = 2

pdf(paste0(out_dir,"/arrowsigsvwpc.pdf"))
plot(wpc$x[,1]/max(abs((wpc$x[,1]))),wpc$x[,2]/max(abs((wpc$x[,2]))), xlim = c(-1,1),ylim = c(-1,1),pch = 1, col = condition_colors[conditions], xlab = "wPC1", ylab = "wPC2")
legend(x = "topleft",legend = levels(conditions),pch = 1,col = condition_colors )

# Sigs Associated with PC1
for (i in 1:NUM_LEVELS){
  sname = rownames(cors)[order(q[,1])][i]
  lm1 = lm(qs[sname,] ~ wpc$x[,1] + wpc$x[,2] )
  vec = c(lm1$coefficients[2],lm1$coefficients[3])
  nvec = vec/sqrt(sum(vec^2))
  arrows(0,0,summary(lm1)$r.squared*nvec[1],summary(lm1)$r.squared*nvec[2], lwd = 2, lty = i, col = 'black')
}

# Sigs Associated with PC2
for (i in 1:NUM_LEVELS){
  sname = rownames(cors)[order(q[,2])][i]
  lm1 = lm(qs[sname,] ~ wpc$x[,1] + wpc$x[,2] )
  vec = c(lm1$coefficients[2],lm1$coefficients[3])
  nvec = vec/sqrt(sum(vec^2))
  arrows(0,0,summary(lm1)$r.squared*nvec[1],summary(lm1)$r.squared*nvec[2], lwd = 2, lty = i, col = 'purple')
}

all.snames = c(rownames(cors)[order(q[,1])][1:NUM_LEVELS],rownames(cors)[order(q[,2])][1:NUM_LEVELS])
all.cols = c(rep("black",NUM_LEVELS),rep("purple",NUM_LEVELS))
all.types = rep(c(1:NUM_LEVELS),2)
legend(x = "bottomleft",legend = all.snames, lty = all.types, col = all.cols ,cex = .5)
dev.off()

## Yosef Plots
source(paste0(lib_dir,"/YosefPlot.R"))
x = exprs(tc.sc.eSet)
YosefPlot(x,gf.vec = gf.vec,as.character(unlist(read.table(housekeeping_list))),draw.lines = T,out.dir = paste0(out_dir,"/yplots"),plot.name = "hk_genes.pdf")

for (sname in c(rownames(cors)[q[,1] == min(q[,1])],rownames(cors)[q[,1] == sort(q[,1])[2]],rownames(cors)[q[,2] == min(q[,2])],rownames(cors)[q[,2] == sort(q[,2])[2]])){
  s1 = GetSigSet(table = sig_file,name = sname)
  YosefPlot(x,gf.vec = gf.vec, s1[,1],draw.lines = T,out.dir = paste0(out_dir,"/yplots"),plot.name = paste0(sname,".pdf"))
}


