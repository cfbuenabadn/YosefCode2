rm(list=ls())
setwd("~/data/users/allonwag/YosefCode/packages/RNASeq/OFBIT")
source("OFBIT.R")

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
c1_run = as.character(phenoData(sf.sc.eSet)$MD_c1_run_id)


e = exprs(sf.sc.eSet)
q = pData(protocolData(sf.sc.eSet))
biobatch = cellType
techbatch = c1_run
hk_names = as.character(as.matrix(read.table(housekeeping_list)))

out.dir = "~/archive/users/allonwag/OFBIT"
out.file= file.path(out.dir, "normreport2.txt")
plot.dir = file.path(out.dir, "normplots")
if (!file.exists(out.dir)){
  dir.create(out.dir)
}
if (!file.exists(plot.dir)){
  dir.create(plot.dir)
}




OFBIT(e = e,type = "TPM",q = q,techbatch = techbatch,biobatch = biobatch ,hk_genes = hk_names %in% rownames(e),out.file=out.file, plot.dir=plot.dir,
      binned.norm.methods=c("LocScale","GlobScale","AdjEig"),
      combat.methods="no") #I removed combat from the norm.methods because the technical covariates so much correlates with the biological ones there's no solution...


