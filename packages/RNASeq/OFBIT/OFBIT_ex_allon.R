rm(list=ls())
setwd("~/data/users/allonwag/YosefCode/packages/RNASeq/OFBIT")
source("OFBIT.R")

housekeeping_list = "/data/yosef/CD8_effector_diff/src/SummaryPipeline/house_keeping_mouse_TitleCase.txt"

EXPERIMENT = "olfa_experiment3"
if(EXPERIMENT == "cortical") {

collect_dir="/data/yosef/BRAIN/processed_July2015/collect"
out_dir = file.path(collect_dir, "pipe")
load(file=paste0(out_dir,"/imageAfterFiltering.RData"))

# cellType = as.character(phenoData(collectedRNASeqStudy$cuff_eSet)$MD_cell_type)
# cellType[cellType == "SCNN1A_(LAYER_IV)"] = "SCNN1A+_(LAYER_IV)"
# cellType[cellType == "RBP4_(LAYER_V)"] = "RBP4+_(LAYER_V)"
# 
# allCellList =  rownames(phenoData(collectedRNASeqStudy$cuff_eSet))
# filteredCellList = colnames(sf.sc.eSet)
# include = match(filteredCellList, allCellList)

# cellType = cellType[include]

cellType = as.character(phenoData(sf.sc.eSet)$MD_cell_type)
cellType[cellType == "SCNN1A_(LAYER_IV)"] = "SCNN1A+_(LAYER_IV)"
cellType[cellType == "RBP4_(LAYER_V)"] = "RBP4+_(LAYER_V)"
biobatch = cellType

} else if (EXPERIMENT == "olfa_experiment3") {
  collect_dir="/data/yosef/BRAIN/processed_June2015_b/collect/olfa_experiment3"
  load(file=paste0(collect_dir,"/imageAfterFiltering.RData"))
  
  experimental_condition = as.character(phenoData(sf.sc.eSet)$MD_expt_condition)
  biobatch = experimental_condition
  
} else if (EXPERIMENT == "olfa_experiment4") {
  collect_dir="/data/yosef/BRAIN/processed_June2015_b/collect/olfa_experiment4"
  load(file=paste0(collect_dir,"/imageAfterFiltering.RData"))
  
  experimental_condition = as.character(phenoData(sf.sc.eSet)$MD_expt_condition)
  biobatch = experimental_condition
} else {
  stop("unrecognized experiment")
}

e = exprs(sf.sc.eSet)
q = pData(protocolData(sf.sc.eSet))

c1_run = as.character(phenoData(sf.sc.eSet)$MD_c1_run_id)
techbatch = c1_run
hk_names = as.character(as.matrix(read.table(housekeeping_list)))

out.dir = "~/archive/users/allonwag/OFBIT"
out.file = file.path(out.dir, "normreport2.txt")
error.file = file.path(out.dir, "normerrors.txt")
log.file = file.path(out.dir, "normlog.txt")
plot.dir = file.path(out.dir, "normplots")
if (!file.exists(out.dir)){
  dir.create(out.dir)
}
if (!file.exists(plot.dir)){
  dir.create(plot.dir)
}


# important note: OFBIT always gets linear scales (i.e., not log scale) data even for TPMs!
# give it the full list of genes - OFBIT has a strong and weak filtering

#I removed combat from the norm.methods because the technical covariates so much correlates with the biological ones there's no solution...
if (file.exists(out.file)) file.remove(out.file)
if (file.exists(error.file)) file.remove(error.file)
if (file.exists(log.file)) file.remove(log.file)
OFBIT(e = e,type = "TPM",q = q,techbatch = techbatch,biobatch = biobatch ,hk_genes =  rownames(e) %in% hk_names ,out.file=out.file, plot.dir=plot.dir,     
      combat.methods="no",
      error.file=error.file,
      log.file=log.file)


#filtering.methods=c("weak","strong"),
#binned.norm.methods=list(), 
#regression.norm.methods=list(),