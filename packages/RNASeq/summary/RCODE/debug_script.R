rm(list=ls())
system("Rscript Pipe.R -a -b -c")


# option_list <- list(
#   make_option("--collect", default="~/data/BRAIN/processed_15_05/150515_HS3A/rsem", type="character",
#               help="Directory containing your RNA Seq results from the preproc pipeline."),
#   make_option("--config", default="~/data/BRAIN/processed_15_05/150515_HS3A/configOlfactory.xls", type="character",
#               help="Config file for your project (.xls or .xlsx)."),
#   make_option("--qcfields", default="~/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt", type="character"),
#   make_option("--genefields", default="~/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt", type="character",
#               help=""),
#   make_option("--out", default="~/data/BRAIN/processed_15_05/150515_HS3A/rsem/summary", type="character",
#               help=""),
#   make_option("--lib", default="~/YosefCode/packages/RNASeq/summary/RCODE", type="character",
#               help=""),
#   make_option("--sigfile", default="~/YosefCode/packages/RNASeq/summary/EXAMPLE/reference_files/immsig.txt", type="character",
#               help=""),
#   make_option("--housekeeping", default="/home/eecs/allonwag/archive/users/allonwag/temp/house_keeping.txt", type="character",
#               help=""),
#   make_option("--combat", action="store_true", default=FALSE,
#               help="This will run the ComBat package for batch correction on your data."),
#   make_option("--multiple_collect", type="character", default=FALSE,
#               help="If you need to load multiple collect directories and config files, please supply a text file listing them here..")
# )
# 
# ## ----- Parse Arguments -----
