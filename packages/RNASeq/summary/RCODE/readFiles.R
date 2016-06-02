rm(list=ls())
setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef/BRAIN/processed_Sep2015/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                              config_file=file.path(collect_dir, "config_olfactory.xlsx"),
                             qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                             gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                             LOAD_RSEM=T, LOAD_CUFF=T)
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


setwd("/data/yosef/users/allonwag//YosefCode//packages//RNASeq//summary//RCODE")
source("loadProcessedRNASeq_NG.R")
collect_dir="/data/yosef2/TFH/processed/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_madeup.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=T, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))