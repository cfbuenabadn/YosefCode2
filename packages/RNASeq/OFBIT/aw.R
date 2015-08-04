collect_dir = "/data/yosef/TFH/processed/collect/"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_tfh.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=F, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path("/data/yosef/TFH/processed/collect", "collectedRNASeqStudy.RData"))

collect_dir="/data/yosef/BRAIN/processed_June2015_b/collect2"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_brain.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=F, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy_olfactory.RData"))

collect_dir="/data/yosef/BRAIN/processed_July2015/collect"
collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
                                              config_file=file.path(collect_dir, "config_cortical.xlsx"),
                                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
                                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
                                              LOAD_RSEM=T, LOAD_CUFF=T)
save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy_withRSEM.RData"))
