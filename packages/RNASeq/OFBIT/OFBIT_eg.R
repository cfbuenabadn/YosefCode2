source("~/OFBIT/OFBIT.R")
e = exprs(sf.sc.eSet)
q = pData(protocolData(sf.sc.eSet))
biobatch = phenoData(sf.sc.eSet)$Condition_Code
techbatch = NULL
hk_names = as.character(as.matrix(read.table(housekeeping_list)))

OFBIT(e = e,type = "TPM",q = q,techbatch = techbatch,biobatch = biobatch ,hk_genes = hk_names %in% rownames(e),out.file = "/data/yosef/users/mbcole/Fluidigm/normreport2.txt",plot.dir = "/data/yosef/users/mbcole/Fluidigm/normplots2")
