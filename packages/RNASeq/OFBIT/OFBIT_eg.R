source("~/YosefCode/packages/RNASeq/OFBIT/OFBIT.R")
load("~/yosef/users/mbcole/WTP63/WTp63_counts_QCscores_someMetaD_GeneLists.rda")
e = WTp63_THcounts
q = t(qualityScores)

biobatch = expt_condition
techbatch = NULL
hk_names = as.character(as.matrix(HKlistD1))

gf.vec = rep(T,dim(e)[1])
OFBIT(e = e,type = "TPM",q = q,techbatch = techbatch,biobatch = biobatch,estimate_fnr = T,tf.vec = gf.vec,
      error.file = "/data/yosef/users/mbcole/WTP63/error_fnr.txt",
      hk_genes =  rownames(e) %in% hk_names,
      out.file = "/data/yosef/users/mbcole/WTP63/normreport_fnr.txt",
      plot.dir = "/data/yosef/users/mbcole/WTP63/normplots_fnr"    )
