# Removes all data in R environment
rm(list=ls())

## ----- Arguments -----
# Consider CRAN library "optparse" to define arguments in command line:
# https://cran.r-project.org/web/packages/optparse/index.html

# Directory containing your RNA Seq results from the preproc pipeline, collected.
# Example: /data/yosef/users/mbcole/Fluidigm/collect2/
collect_dir = "/data/yosef/users/mbcole/Fluidigm/collect2/"

# Config file for your project (.xls or .xlsx).
# This file contains all of the phenotype / experimental details for the study at hand.
# Example file: /data/yosef/users/mbcole/Fluidigm/Fluidigm_config_NG.xlsx
# Required fields:
#         sample_id = id unique to the sample
#         sample_sequencing_id = id unique to the library ( a sample may be sequenced multiple times )
#           These are the names of the columns in the resulting expression matrix.
#         output_name = unique library names the come out of the collect script (in cell_list.txt)
#           This field allows the load method (see below) to map all additional fields in the config file to the expression values for that library.
# Note that these are the first three fields in the example excel document, but order does not matter.
config_file = "/data/yosef/users/mbcole/Fluidigm/Fluidigm_config_NG.xlsx"

# Output directory
# Example: /data/yosef/users/mbcole/Fluidigm/out_test
out_dir = "/data/yosef/users/mbcole/Fluidigm/out_test"

# Constant: Names of QC fields evaluated in preproc pipeline.
qc_fields_file = "/home/eecs/mbcole/YosefCode/packages/RNASeq/summary/TEXT/qc_fields.txt"

# Constant: Names of gene feature fields evaluated in preproc pipeline.
gene_fields_file = "/home/eecs/mbcole/YosefCode/packages/RNASeq/summary/TEXT/gene_fields.txt"

# Constant: Path to R Script Files
lib_dir = "/home/eecs/mbcole/YosefCode/packages/RNASeq/summary/RCODE"

# Other: Path to housekeeping gene list: List of gene symbols corresponding to transcripts that show evidence of
# "robust" expression in population of interest: constitutive expression / minimal expression heterogeneity
# Human example: /home/eecs/mbcole/YosefCode/packages/RNASeq/summary/EXAMPLE//reference_files//house_keeping_human_names.txt
# Mouse example: /data/yosef/CD8_effector_diff/src/SummaryPipeline/house_keeping_mouse_TitleCase.txt
housekeeping_list = "/home/eecs/mbcole/YosefCode/packages/RNASeq/summary/EXAMPLE//reference_files/house_keeping_human_names.txt"

## ----- Produce Output Directory -----
if (file.exists(out_dir)){
  # Do Nothing
} else {
  dir.create(out_dir)
}

## ----- Load All Expression Data (RSEM, Cufflinks, Kallisto?) -----
# Loads all of the expression results into R
source(paste0(lib_dir,"/loadProcessedRNASeq_NG.R"))
dataSet = loadProcessedRNASeq_NG(collect_dir, config_file, qc_fields_file, gene_fields_file, LOAD_RSEM=T, LOAD_CUFF=T, LOAD_KALLISTO=F)

## ----- Select RSEM TPM results - in an eSet object -----
eSet = dataSet$rsem_eSet
# Rsem results contain expectedCounts_table, exprs, fpkm_table, tpm_table 
# Sets the default expression matrix to the tpm table
exprs(eSet) <- assayData(eSet)$tpm_table

## ----- Pre-Filtering of Failed Trancripts -----
# Remove all trancripts that fail in all samples (NA in all samples)
# These can occur when the reference identifies a transcript by name, but does not provide the necessary intervals for alignment.
is.bad_transcript = apply(is.na(exprs(eSet)), 1, all)

print(paste(sum(is.bad_transcript ),"transcripts failed pre-processing:"),quote = F)
print(matrix(rownames(eSet)[is.bad_transcript],ncol = 1),quote = F)

prefilt.eSet = eSet[!is.bad_transcript,]
rm(eSet)

## ----- Pre-Filtering of Failed Samples -----
# Remove all samples failing the preprocessing step:
#   1) Sample has all NA (or 0) values - all transcripts are quantified as NA for a failed run
#   2) All quality info is missing for the sample - all measures are NA for that sample
is.failed = apply(is.na(exprs(prefilt.eSet)) | (exprs(prefilt.eSet) == 0), 2, all) |
  apply(is.na(pData(protocolData(prefilt.eSet))), 1, all)

print(paste(sum(is.failed),"samples failed pre-processing:"),quote = F)
print(matrix(colnames(prefilt.eSet)[is.failed],ncol = 1),quote = F)

prefilt.eSet = prefilt.eSet[,!is.failed]

print("=========================",quote = F)
print("Removed failed samples.",quote = F)
print(paste(sum(!is.failed),"samples retained."),quote = F)

## ----- Study-Specific Sample Pre-Filtering -----
# Here you can use fields supplied in the config file to do additional filtering necessary for your study.
# Example: is.study = phenoData(prefilt.eSet)$Coverage_Type == "High"
is.study = T
prefilt.eSet = prefilt.eSet[,is.study]
print(paste(dim(prefilt.eSet)[2],"samples retained for study."),quote = F)

## ----- Study-Specific Transcript Pre-Filtering
# Select only type 1 transcripts (Coding), that are detected in at least one sample
is.study.transcript = (featureData(prefilt.eSet)$Transcript_Type == "protein_coding") & (apply((exprs(prefilt.eSet)) > 0, 1, any))
prefilt.eSet = prefilt.eSet[is.study.transcript,]
print(paste(sum(is.study.transcript ),"transcripts retained for study."),quote = F)

## ----- Checkpoint: No NA shall pass! -----
stopifnot(!any(is.na(exprs(prefilt.eSet))))

sc.eSet = prefilt.eSet
rm(prefilt.eSet)

#################################
## ----- Filtering Stage ----- ##
#################################

source(paste0(lib_dir,"/GeneFilter.R")) # Gene Filtering
source(paste0(lib_dir,"/SampleFilter.R")) # Basic Sample Filtering
source(paste0(lib_dir,"/TechFilter.R")) # Data-Driven Sample Filtering

count.cutoff = 10 # Genes falling BELOW this level...
prop.failed = .85 # In MORE than this fraction of cells...
                  # Will be removed.

## ----- Preliminary Gene Filtering
# Initialize the gene filter vector, add plots to subdirectory
init.gf.vec = GeneFilter(sc.eSet,
                         count.cutoff = count.cutoff,
                         prop.failed = prop.failed,
                         verbose = T,
                         plot.dir = paste0(out_dir,"/genefilter"),
                         plot.prefix = "pre_cell_filtering")

## ----- Create Sample Filtering Directories
if (!file.exists(paste0(out_dir,"/samplefilter/"))){dir.create(paste0(out_dir,"/samplefilter/"))}
if (!file.exists(paste0(out_dir,"/samplefilter_tech/"))){dir.create(paste0(out_dir,"/samplefilter_tech/"))}

## ----- Hard-Cutoff (Cell) Sample Filtering + Correlation with Technical Features
# The Hard-Cutoff module is for diagnostic purposes only - output is not used. Take a look at the graphical output.
# Basic Sample Filtering: Hard mode means Z_CUTOFF=NULL -> Hard cutoffs will be used.
hard.sf.sc.eSet = SampleFilter(eSet = sc.eSet,gene.filter.vec = init.gf.vec, housekeeping_list = housekeeping_list,mixture = T,verbose = T, MIN_RALIGN = 5, Z_CUTOFF=NULL, plot.dir = paste0(out_dir,"/samplefilter/hard_cutoff") )
# Update gene vector
hard.gf.vec = GeneFilter(hard.sf.sc.eSet,
                         count.cutoff = count.cutoff,
                         prop.failed = prop.failed,
                         verbose = T,
                         plot.dir = paste0(out_dir,"/genefilter"),
                         plot.prefix = "post_hard_cell_filtering")
# Data-Driven Sample Filtering: Hard mode means Z_CUTOFF=NULL -> No filtering is performed.
hard.tf.sc.eSet = TechFilter(hard.sf.sc.eSet,gf.vec = hard.gf.vec,Z_CUTOFF = NULL,PROP_CUTOFF = .9,plot.dir = paste0(out_dir,"/samplefilter_tech/hard_cutoff"),plot.prefix = "")

## ----- Adaptive Sample and Gene Filtering
# The Adaptive module is the true run of all filtering schemes.
sf.sc.eSet = SampleFilter(eSet = sc.eSet,gene.filter.vec = init.gf.vec, housekeeping_list = housekeeping_list,mixture = T,verbose = T,MIN_RALIGN = 5, plot.dir = paste0(out_dir,"/samplefilter/adaptive") )
# Update gene vector
gf.vec = GeneFilter(sf.sc.eSet,
                    count.cutoff = count.cutoff,
                    prop.failed = prop.failed,
                    verbose = T,
                    plot.dir = paste0(out_dir,"/genefilter"),
                    plot.prefix = "post_adapt_cell_filtering_1")
tf.sc.eSet = TechFilter(sf.sc.eSet,gf.vec = gf.vec,PROP_CUTOFF = .9,plot.dir = paste0(out_dir,"/samplefilter_tech/adaptive"),plot.prefix = "")

# Final update of gene vector
gf.vec = GeneFilter(tf.sc.eSet,
                    count.cutoff = count.cutoff,
                    prop.failed = prop.failed,
                    verbose = T,
                    plot.dir = paste0(out_dir,"/genefilter"),"post_adapt_cell_filtering_2")

#####################################
## ----- Normalization Stage ----- ##
#####################################

# ...

