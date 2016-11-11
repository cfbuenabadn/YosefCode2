#for debug
rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/archive/users/allonwag/temp/big_pipe_out/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_file.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=F, LOAD_CUFF=T, LOAD_KALLISTO=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))



rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
#collect_dir="/data/yosef/BRAIN/processed_Sep2015/collect"
collect_dir="/data/yosef2/BRAIN/processed_olfactory_Jun2016/collect_20160818/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                              config_file=file.path(collect_dir, "config_olfactory.xlsx"),
                             qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                             gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                             LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef/BRAIN/processed_July2015/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_cortical.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=T, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))

rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef2/BRAIN/processed_cortical_Oct2016/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_cortical.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef/BRAIN/processed_Bateup_Aug2015/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_bateup.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=T, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef/BRAIN/processed_Zebrafish_Oct2015/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_samisrael.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=T, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data2/Published_Data/TH17/processed_aw20160712/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_th17.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))



rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data2/Published_Data/Shalek_DC/co/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path("~/data2/Published_Data/Shalek_DC/config_shalek2014_aw.txt"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=FALSE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data2/TFH/processed_20160720/collect/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_tfh.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))




rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data2/TFH/processed_20160720/collectWithFateMapping/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_tfh.csv"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))






rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data/TFH/processed2/collect/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_FC_01930.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=FALSE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))


rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data/TFH/processed_20161012/FC_01930/collect/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_FC_01930.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=FALSE, LOAD_CUFF=FALSE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))



rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="~/data2/BASF/Nutraceuticals/processed_RNASeq_20160826/collect/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_basf.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=TRUE, LOAD_CUFF=TRUE, LOAD_KALLISTO=TRUE)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))
