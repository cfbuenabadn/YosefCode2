# Data Cleaning Module 0: Load TPM Data and QC Summary from rsem Preproc Pipeline
# Michael Cole, March 2015
# -------------------------

library(Biobase)
library(gdata)
library(xlsx)

# Function: loadRSEM
# Usage: Load tpm data from collect script output <- RNASeq preprocessing pipeline w/rsem
# Params:
# -------------------------
# collect_dir = path to collect output directory
#   Must contain 3 files:
#     qc_table.txt, gene_list.txt, rsem_tpmTable.txt
# config_file = path to config file containing info on all samples
# qc_fields_file = path to names of qc fields (one per line - no tabs!)
# gene_fields_file = path to names of gene info fields (one per line - no tabs!)

loadRSEM = function(collect_dir,config_file,qc_fields_file,gene_fields_file){

  ##----- Loading Sample Quality Data
  # Note - No QC fields available!
  qc_table = as.data.frame(t(as.matrix(read.table(paste0(collect_dir,"/qc_table.txt"), sep = "\t",header = F))))
  colnames(qc_table) = as.character(unlist(read.table(qc_fields_file, header = F,sep = "\t")))
  rownames(qc_table) = as.character(unlist(read.table(paste0(collect_dir,"/cell_list.txt"), sep = "\t",header = F)))

  ##----- Loading Gene Info
  # Note - No Gene Fields available!
  gene_info = read.table(paste0(collect_dir,"/gene_list.txt"), header = F, sep = "\t")
  colnames(gene_info) = as.character(unlist(read.table(gene_fields_file, header = F, sep = "\t")))
  rownames(gene_info) = gene_info$Gene_Symbol

  ##----- Loading TPM Data
  tpm_table = read.table(paste0(collect_dir,"/rsem_tpmTable.txt"), header = F, sep = "\t")
  rownames(tpm_table) = rownames(gene_info)
  colnames(tpm_table) = rownames(qc_table)

  ##----- Phenotype Data
  suff = gsub(".*\\.","",config_file)
  if (suff == "xls"){
    config_table = read.xls(config_file,sheet = 1)
    config_key = read.xls(config_file,sheet = 2)
  }else if (suff == "xlsx"){
    config_table = read.xlsx(config_file,sheetIndex = 1)
    config_key = read.xlsx(config_file,sheetIndex = 2)
  } else{
    stop("Configuration file must be .xls or .xlsx.")
  }

  # Decode Configuration
  for (level in levels(config_key$Code_Class)){
    codes = config_table[,colnames(config_table) == level]
    map = config_key[config_key$Code_Class == level,]
    rownames(map) = map$Code
    desc = map[as.character(codes),]$Description
    config_table[,colnames(config_table) == level] = as.factor(as.character(desc))
  }

  # Names assigned by collect are based on preprocessing directory names
  rownames(config_table) = gsub("(.*/)([^/]*/[^/]*$)","\\2",config_table$Preproc_Dir)
  sample.info = config_table[colnames(tpm_table),]
  
  ##----- Generate eSet
  protoDat = new("AnnotatedDataFrame", data = qc_table)
  featureDat = new("AnnotatedDataFrame", data = gene_info)
  phenoData = new("AnnotatedDataFrame", sample.info)
  eMatrix = data.matrix(tpm_table)
  eSet = new("ExpressionSet", exprs = eMatrix, featureData = featureDat, protocolData = protoDat, phenoData = phenoData)
  return(eSet)
}

# Function: loadRSEMStudy
# Usage: Load tpm data from multiple collect script outputs
# Params:
# -------------------------
# muliple_collect = path to tab-delim. table containing paths to config files in first column and collect directories in second column
# qc_fields_file = path to names of qc fields (one per line - no tabs!)
# gene_fields_file = path to names of gene info fields (one per line - no tabs!)

loadRSEMStudy = function(muliple_collect,qc_fields_file,gene_fields_file){
  dfCollect <- read.table(multiple_collect,header=T,sep="\t")
  
  print("Combining multiple collects:")
  apply(dfCollect,1,function(x) print(x[2]))
  
  # Call loadRSEM on each collect
  li_eSet <-apply(dfCollect,1,function(x) loadRSEM(collect_dir = x[2],config_file =x[1],qc_fields_file = qc_fields_file,gene_fields_file = gene_fields_file))
  
  eSet = li_eSet[[1]]
  if (length(li_eSet) > 1){
    for (i in 2:length(li_eSet)){
      eSet <- Biobase::combine(eSet,li_eSet[[i]])
    }
  }
  return(eSet)
}