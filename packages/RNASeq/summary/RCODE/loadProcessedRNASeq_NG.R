# revision history
# Allon Wagner, July 2015
# Michael Cole, March 2015
# -------------------------

library(Biobase)
library(gdata)
library(xlsx)

loadProcessedRNASeq_NG = function(collect_dir, config_file=file.path(collect_dir, "config_brain.xlsx"), 
                                  qc_fields_file, gene_fields_file, 
                                  LOAD_RSEM=T, LOAD_CUFF=T, LOAD_KALLISTO=T)
{
  if(LOAD_RSEM)
  {
    print("Loading RSEM results...")
    rsem_eSet = loadRSEM(file.path(collect_dir, "rsem"), config_file, qc_fields_file, gene_fields_file);
    print("RSEM results loaded successfully")
  }
  else
  {
    rsem_eSet = NULL
  }
  
  if(LOAD_CUFF)
  {
    print("Loading cufflinks results...")
    cuff_eSet = loadCuff(file.path(collect_dir, "cuff"), config_file, qc_fields_file, gene_fields_file);
    print("Cufflinks results loaded successfully")
  }
  else
  {
    cuff_eSet = NULL
  }
    
  #Kallisto pipeline is not implemented yet
  if(LOAD_KALLISTO)
  {
    print("Loading Kallisto results...")
    kallisto_eSet = loadKallisto(file.path(collect_dir, "kallisto"), config_file, qc_fields_file, gene_fields_file);
    print("Kallisto results loaded successfully")
  }
  else
  {
    kallisto_eSet = NULL
  }
  
  print("All data loaded successfully")
  return(list("rsem_eSet"=rsem_eSet, "cuff_eSet"=cuff_eSet, "kallisto_eSet"=kallisto_eSet))
}

loadKallisto = function(collect_dir ,config_file, qc_fields_file, gene_fields_file)
{
  commonOutput = loadCommonPreprocOutput(collect_dir ,config_file, qc_fields_file, gene_fields_file)
  
  ##----- Loading TPM Data
  tpm_table = loadExpressionMatrix(file.path(collect_dir, "kallisto_tpmTable.txt"), commonOutput)
  print("Kallisto TPM table loaded successfully")
  
  ##----- Loading estimated counts Data
  estimatedCounts_table = loadExpressionMatrix(file.path(collect_dir, "kallisto_readCountsTable.txt"), commonOutput)
  print("Kallisto expected counts table loaded successfully")
  
  
  ##----- Generate ExpressionSet
  protoDat = new("AnnotatedDataFrame", data = commonOutput$qc_table)
  featureDat = new("AnnotatedDataFrame", data = commonOutput$gene_info)
  phenoData = new("AnnotatedDataFrame", commonOutput$config_table)
  assayData <- new.env(parent = emptyenv())
  assayData$tpm_table = data.matrix(tpm_table)
  assayData$estimatedCounts_table = data.matrix(estimatedCounts_table)
  #by default in Kallisto: use TPMs for the expression
  assayData$exprs = assayData$tpm_table
  eSet = ExpressionSet(assayData, featureData = featureDat, protocolData = protoDat, phenoData = phenoData)
  
  return(eSet) 
}

loadCuff = function(collect_dir ,config_file, qc_fields_file, gene_fields_file)
{
  commonOutput = loadCommonPreprocOutput(collect_dir ,config_file, qc_fields_file, gene_fields_file)
  
  ##----- Loading FPKM Data
  fpkm_table = loadExpressionMatrix(file.path(collect_dir, "cuff_fpkmTable.txt"), commonOutput)
  print("Cufflinks FPKM table loaded successfully")
  
  ##compute TPM data for cufflinks based on the FPKM data
  ##see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
  ##for the formula: TPM_i = (FPKM_i / (sum_j{FPKM_j})) * 10^6
  tpm_table = sweep(fpkm_table, 2, colSums(fpkm_table, na.rm=T), '/') * (10^6)
  
  #load counts (real counts, produced by featureCounts, in contrast to cufflink's expected counts)
  counts_table = loadExpressionMatrix(file.path(collect_dir, "tophat2_featureCountsTable.txt"), commonOutput)
  print("FeatureCounts table loaded successfully")
  
  ##----- Generate ExpressionSet
  protoDat = new("AnnotatedDataFrame", data = commonOutput$qc_table)
  featureDat = new("AnnotatedDataFrame", data = commonOutput$gene_info)
  phenoData = new("AnnotatedDataFrame", commonOutput$config_table)
  assayData <- new.env(parent = emptyenv())
  assayData$fpkm_table = data.matrix(fpkm_table)
  assayData$tpm_table = data.matrix(tpm_table)
  assayData$counts_table = data.matrix(counts_table)
  #by default: use FPKMs for the expression
  #assayData$exprs = assayData$fpkm_table
  #decided to change the default to TPM!
  assayData$exprs = assayData$tpm_table
  eSet = ExpressionSet(assayData, featureData = featureDat, protocolData = protoDat, phenoData = phenoData)
  
}

loadRSEM = function(collect_dir ,config_file, qc_fields_file, gene_fields_file)
{
  commonOutput = loadCommonPreprocOutput(collect_dir ,config_file, qc_fields_file, gene_fields_file)
    
  ##----- Loading TPM Data
  tpm_table = loadExpressionMatrix(file.path(collect_dir, "rsem_tpmTable.txt"), commonOutput)
  print("RSEM TPM table loaded successfully")
  
  ##----- Loading FPKM Data
  fpkm_table = loadExpressionMatrix(file.path(collect_dir, "rsem_fpkmTable.txt"), commonOutput)
  print("RSEM FPKM table loaded successfully")
  
  ##----- Loading expected counts Data
  expectedCounts_table = loadExpressionMatrix(file.path(collect_dir, "rsem_readCountsTable.txt"), commonOutput)
  print("RSEM expected counts table loaded successfully")
  
  
  ##----- Generate ExpressionSet
  protoDat = new("AnnotatedDataFrame", data = commonOutput$qc_table)
  featureDat = new("AnnotatedDataFrame", data = commonOutput$gene_info)
  phenoData = new("AnnotatedDataFrame", commonOutput$config_table)
  assayData <- new.env(parent = emptyenv())
  assayData$tpm_table = data.matrix(tpm_table)
  assayData$fpkm_table = data.matrix(fpkm_table)
  assayData$expectedCounts_table = data.matrix(expectedCounts_table)
  #by default in RSEM: use TPMs for the expression
  assayData$exprs = assayData$tpm_table
  eSet = ExpressionSet(assayData, featureData = featureDat, protocolData = protoDat, phenoData = phenoData)
  
  return(eSet)
}

loadExpressionMatrix = function(fileName, commonOutput)
{
  expressionTable = read.table(fileName, header = F, sep = "\t")
  
  #exclude samples that did not appear in the config file:
  expressionTable = expressionTable[,!commonOutput$samplesInCollectedDataButNotInConfig]
  
  rownames(expressionTable) = rownames(commonOutput$gene_info)
  colnames(expressionTable) = rownames(commonOutput$config_table)#commonOutput$cell_list
  
  return(expressionTable)
}

#Loads QC and gene info that have the same format in all branches of the pipeline
loadCommonPreprocOutput = function(collect_dir ,config_file, qc_fields_file, gene_fields_file)
{
  
  #Load list of cells in the collected data
  cell_list = as.character(unlist(read.table(file.path(collect_dir, "cell_list.txt"), header = F)))
  
  
  ##----- Loading Sample Quality Data
  # Note - No QC fields available!
  qc_table = as.data.frame(t(as.matrix(read.table(paste0(collect_dir,"/qc_table.txt"), sep = "\t",header = F))))
  colnames(qc_table) = as.character(unlist(read.table(qc_fields_file, header = F,sep = "\t")))
  #observe that the rownames of qc_table are set below, after we read the config_table and sort it to be in the same order as the collected data
  
  
  ##----- Loading Gene Info
  # Observe that there are 3 comment line before the beginning of the actual table
  gene_info = read.table(paste0(collect_dir,"/gene_list.txt"), header = F, sep = "\t",skip=3)
  colnames(gene_info) = as.character(unlist(read.table(gene_fields_file, header = F, sep = "\t")))
  
  #at first I thought that the rownames must be the Gene_ID and not the gene symbol because the gene symbol is not unique... (there are gene symbols that are associated with more than one gene ID)
  #ALAS, it turns out that (at least in mm10's cufflinks dictionary) gene IDs are not unique as well
  #so I resorted to using gene
  rownames(gene_info) = make.unique(as.character(gene_info$Gene_Symbol))
  
  
  
  ##----- Phenotype Data
  suff = gsub(".*\\.","",config_file)
  if (suff == "xls"){
    config_table = read.xls(config_file, sheet=1)
  }else if (suff == "xlsx"){
    #Allon: I use read.xlsx2 instead of read.xlsx because read.xlsx may be  more general but it's soooooo slow
    config_table = read.xlsx2(config_file, sheetIndex=1)
  } else{
    stop("Configuration file must be .xls or .xlsx.")
  }
  
  #name the rows after the sequencing ID of each sample
  rownames(config_table) = config_table$sample_sequencing_id;
  
  #there are rows in the collected data that do not appear in the config --> discard them
  samplesInCollectedDataButNotInConfig = !(cell_list %in% config_table$output_name);
  if (any(samplesInCollectedDataButNotInConfig))
  {
    print(sprintf("%d / %d samples appear in the collected QC table but not in the config file - excluding them...", sum(samplesInCollectedDataButNotInConfig), length(samplesInCollectedDataButNotInConfig)))
    cell_list = cell_list[!samplesInCollectedDataButNotInConfig]
    qc_table = qc_table[!samplesInCollectedDataButNotInConfig, ]
  }

  #there are rows in the config that do not appear in the collected data --> a show-stopper!
  #this indicates an error that the user would probably like to fix before proceeding
  samplesInConfigButNotInCollectedData = !(config_table$output_name %in% cell_list);
  if(any(samplesInConfigButNotInCollectedData))
  {
    stop(sprintf("There are %d samples that appear in the config file but do not appear either in the collected data!", sum(samplesInConfigButNotInCollectedData)))
  }
  

  # order the config file in the same order as the collected data - we've already verified above they have the exact same list of cells
  config_table = config_table[match(cell_list, config_table$output_name),]
  
  #once the config_table is in the same order of the collected data (which is also the order of the QC table)
  #we can make the uniqueIDs also the sample names of the QC table
  rownames(qc_table) = config_table$sample_sequencing_id;
  
  
  KEEP_ONLY_MD_COLS = F #at first I kept only metadata cols, but Russell asked to keep all cols
  if(KEEP_ONLY_MD_COLS) {
  #keep only metadata (MD) fields from the config file.
  #You don't need the unique IDs --> these are the row names of the table
  #you don't need the output file name --> we've already ordered the config table to be at the same order as the collected data (i.e., the order of cell_list)
  is_metadata_col = startsWith(colnames(config_table), "MD_", trim=FALSE, ignore.case=FALSE)
  config_table = config_table[,is_metadata_col]
  }
  
  #the samplesInCollectedDataButNotInConfig logical vector is also returned so that the calling function will be able to exclude cells as well from the rsem-specific, cuff-specific etc. matrices that it reads 
  return(list("config_table"=config_table, "qc_table"=qc_table, "gene_info"=gene_info, "samplesInCollectedDataButNotInConfig"=samplesInCollectedDataButNotInConfig))
}

# example for debug:
# collect_dir="/data/yosef/BRAIN/processed_June2015_b/collect"
# collectedRNASeqStudy = loadProcessedRNASeq_NG(collect_dir=collect_dir,
#                               config_file=file.path(collect_dir, "config_brain.xlsx"),
#                              qc_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt",
#                              gene_fields_file="/data/yosef/CD8_effector_diff/src/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt",
#                              LOAD_RSEM=F, LOAD_CUFF=T)
# save(collectedRNASeqStudy, file=file.path(collect_dir, "collectedRNASeqStudy.RData"))
