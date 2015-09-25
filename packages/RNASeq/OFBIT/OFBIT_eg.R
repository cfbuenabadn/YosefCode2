source("~/YosefCode/packages/RNASeq/OFBIT/OFBIT.R")
e = exprs(sf.sc.eSet)
q = pData(protocolData(sf.sc.eSet))
biobatch = phenoData(sf.sc.eSet)$Condition_Code
techbatch = NULL
hk_names = as.character(as.matrix(read.table(housekeeping_list)))

OFBIT(e = e,type = "TPM",q = q,techbatch = techbatch,biobatch = biobatch,estimate_fnr = T,tf.vec = gf.vec,
      error.file = "/data/yosef/users/mbcole/Fluidigm/error.txt",
      hk_genes =  rownames(e) %in% hk_names,
      out.file = "/data/yosef/users/mbcole/Fluidigm/normreport22.txt",
      plot.dir = "/data/yosef/users/mbcole/Fluidigm/normplots22"    )
