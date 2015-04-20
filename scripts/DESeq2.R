###############################################################
# Jim Kaminski
# 2/25/2015
rm(list=ls())
.libPaths('/home/eecs/jimkaminski/R_libs')
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("pasilla")

library(optparse)
library(gridExtra)
library(stringr)
library(ggplot2)
library(reshape)
library(DESeq2)

#library("pasilla")
library("Biobase")


option_list <- list(
  make_option("--data", default="/data/yosef/CD8_effector_diff/out/ATAC-Seq/Peaks/Chronic_vs_Acute/FinalMerged.tab", type="character",
              help="Enter the count data. Program assumes first three columns are chr,start,end, and
              remaining columns are samples, with cells containing the overlap."),
  make_option("--conditions", default="/data/yosef/CD8_effector_diff/data/PeakComparisons/ChronicVsAcuteConditions.tab", type="character",
              help="Tab file of conditions for comparison."),
  make_option("--out", default="~/DESeq2Results.tab", type="character",
              help="Output of DESeq2 Results."),
  make_option("--out_bed_baseline", default="~/BedBaseline.bed", type="character",
              help="Bed file of peaks -upregulated- in condition 1."),
  make_option("--out_bed_alt", default="~/BedAlt.bed", type="character",
                          help="Bed file of peaks -upregulated in condition2."),
  make_option("--baseline", default="Acute", type="character",
              help="The condition you wish to use the as denominator in the logfold comparisons.")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

###############################################################################
# Load data, make into DESeq Matrix
print("Loading Count Data into R...")
dfCountData <- read.table(opt$data,header=T)
dfConditions <- read.table(opt$conditions,header=T)

print(dim(dfConditions))
print(dfConditions$Sample )
colnames(dfCountData) <- gsub("X","_",colnames(dfCountData))

# Limit list of samples in dfConditions to those that are in dfCountData
dfConditions <- dfConditions[dfConditions$Sample %in% colnames(dfCountData),]

rownames(dfCountData) <- paste(dfCountData$Chr,dfCountData$Start,dfCountData$End,sep="-")
dfCountData <- dfCountData[,c(4:dim(dfCountData)[2])]


print("Making DESeq Dataset...")
dedsSEQ <- DESeqDataSetFromMatrix(countData = dfCountData,
                              colData = dfConditions,
                              design = ~ Condition)

###############################################################################

# Set the baseline conditions
dedsSEQ$Condition <- relevel(dedsSEQ$Condition, opt$baseline)

# Run the analysis
dedsSEQResults <- DESeq(dedsSEQ)
deRes <- results(dedsSEQResults)
deRes <- deRes[order(deRes$padj),]

head(deRes)
table(complete.cases(deRes))
liSplit <- str_split(rownames(deRes), pattern="-", n = Inf)
dfRes <- as.data.frame(matrix(unlist(liSplit),nrow=dim(deRes)[1],byrow=T))
colnames(dfRes) <-c("chr","start","end")
dim(deRes)


dfResults <- cbind(dfRes,deRes)
dfInvestigate <- dfResults[!complete.cases(dfResults),]
dfResults <- dfResults[complete.cases(dfResults),]

vPvals <- sort(unique(dfResults$padj))
vCount <- lapply(vPvals, function(x) sum(dfResults$padj < x))
vCount <- unlist(vCount)
dfPValPlot <- cbind(vPvals,vCount)
plot(vPvals,vCount)

dfResults <- dfResults[dfResults$padj<0.05,]
write.table(dfResults,opt$out,sep="\t",row.names=F,quote=FALSE)

dfResults <-dfResults[with(dfResults, order(chr, start,end)), ]
write.table(dfResults[dfResults$log2FoldChange>0,c(1,2,3)],opt$out_bed_alt,sep="\t",row.names=F,quote=FALSE,col.names=F)
write.table(dfResults[dfResults$log2FoldChange<0,c(1,2,3)],opt$out_bed_baseline,sep="\t",row.names=F,quote=FALSE,col.names=F)
print("Wrote tables.")


