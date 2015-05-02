##############################################################################
# Jim Kaminski
# Andre King assisted with the config_file option
# (based in part on some code originally written by Faraz Tavakoli)
# 2/27/2015
# Yosef Lab
##############################################################################
"""
This program is designed to take a file of peaks from MACS2
and several "bed" files of annotations on the genome. It identifies overlap
between the peaks and annotations, and conducts a binomial test on each
annotation to determine if it is "enriched" in the peaks.

Input:
* strFileConfig - path to config file. This is a two column tab delimited file.
The first column lists paths to data files (peaks or annotations), the second
column can contain "peaks","0","1". 
    peaks - MACS2 output
    0 - Annotation file, do not include
    1 - Annotation file, include

Output
* strOverlapResults - This is a tab delimited matrix. Each row is a peak, each
column is an annotation. "1" indicates overlap.

* strTestResults -  This is a tab delimited matrix. Each row is a peak, each
column is an annotation. "1" indicates overlap.

"""

"""
Remaining Tasks:
* Finish setting it up for gff files.
* Write a function to convert Nir's db's to bed files.
* Print out pvalues from binomial test.

Note: Some of our gff files are missing a newline, and have "chr" capitalized.
 You can fix them by using sed:
 sed -i 's/.Chr/\nchr/g' /data/yosef/index_files/mm9/other_annot/regulatory_features.gff
 sed -i 's/ESNT_/ES\nNT_/g' /data/yosef/index_files/mm9/other_annot/regulatory_features.gff
"""

import numpy as np
import subprocess as sp
import scipy.stats as stats
import re
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os
import argparse
import src.Overlaps_Functions as odb

parser = argparse.ArgumentParser(description='')
grpParam = parser.add_argument_group('')
grpParam.add_argument('--config', type=str, dest='strConfig',default="", 
                      help='Please enter the full path to the config file. This is a two column tab delimited file. \
                      The first column lists paths to data files (peaks or annotations), the second \
                      column can contain "peaks","0","1". \
                          peaks - MACS2 output \
                          0 - Annotation file, do not include \
                          1 - Annotation file, include \
                          NirDB - Nir\'s matlab style database.')
grpParam.add_argument('--out', type=str, dest='dirOut', help='Output folder',default=os.getcwd()+os.sep+"out")
grpParam.add_argument('--bed', type=str, dest='dirBed', help='Folder for bed files',default= "")
grpParam.add_argument('--peaks', type=str, dest='astrPeaks', nargs='*',help='Enter one or more peak files, saved as bed files.',default= "")
grpParam.add_argument('--bedtools', type=str, dest='cmdBedtools', help='Path to bedtools program',default= "/home/eecs/jimkaminski/tools/bedtools2/bin/bedtools")
args = parser.parse_args()

################################################
# Create Output Folders

dirOut = args.dirOut
dirTmp = args.dirOut + os.sep + "tmp"
odb.check_create_dir( dirOut )
odb.check_create_dir( dirTmp )
if args.dirBed == "":
    dirBed = dirOut + os.sep + "bed"
else:
    dirBed = args.dirBed
odb.check_create_dir( dirBed )

############################################################################
# Read in the list of peak, bed, and gff files from the config file.


astrPeaks = args.astrPeaks
adfFullResults = []
astrBeds = []
astrGFFs = []



astrNirDBs = []
astrUnknown = []

with open(args.strConfig, 'rb') as fileConfig:
    for astrLine in fileConfig:
        
        if astrLine.strip() != "" and astrLine[0]!="#":
            print astrLine
            strPath = astrLine.split('\t')[0].strip()
            strType = astrLine.split('\t')[1].strip()            
        
            mtchBed = re.search(r'.bed',strPath)
            mtchGff = re.search(r'(.gff)|(.gtf)',strPath)
        
            
            if strType=="peaks":
                astrPeaks.append(strPath)
            elif strType=="1" and mtchBed:
                astrBeds.append(strPath)
            elif strType=="1" and mtchGff:
                astrGFFs.append(strPath)
            elif strType=="1" and not mtchBed:
                astrUnknown.append(strPath)
            elif strType=="NirDB":
                astrNirDBs.append(strPath)

print astrBeds

##############################################################################
# Make bed files from files in "Nir DB" format.

for strDir in astrNirDBs:
    for strBed in odb.ReadMatlabDB(strDir,dirBed):
        astrBeds.append(strBed)
    


##############################################################################
# Run each annotation file against peak files.

# At the end, each dfFullTable is a pandas dataframe with each row as a peak,
# and chr,start, and end as the first three columns, and each bed/gff with 
# overlap as additional columns. Each cell contains the bp of overlap with the
# peak.

dictPeakdfs = {}
dictCutPeaks = {}

for strPeaks in astrPeaks:
    strCutPeaks = dirTmp + os.sep + os.path.basename(strPeaks).replace(".bed","cut.bed")
    strStem = strPeaks
    
    odb.cut_file(strPeaks,strCutPeaks)
    dictCutPeaks[strPeaks] = strCutPeaks
    dfPeaks = pd.read_table(strCutPeaks,sep="\t",names=["chr","start","end"])
    dfPeaks.columns = ["chr","start","end"]

    #We want to reflect this line of unix code:
    # cut -f1-3 /data/yosef/index_files/mm9/other_annot/Gata/CD8.bed | 
    # bedtools intersect -filenames -a reduced.bed -b stdin | sed -e 's/$/\tCD8/'
    

    strAllHits = dirTmp + os.sep + "AllHits.bed"    
    with open(strAllHits,"w") as fileHits:
        for strBed in astrBeds:
            odb.IntersectBED( strBed, strCutPeaks,fileHits,args.cmdBedtools,strBed,"bed against peaks" )
        for strGFF in astrGFFs:
            odb.IntersectGFF( strGFF, strCutPeaks,fileHits,args.cmdBedtools,strGFF,"gff against peaks" )
                
    dictPeakdfs[strPeaks] = odb.ConvertOverlapHits(strAllHits,dirTmp,dfPeaks,os.path.basename(strPeaks))
    os.remove(strAllHits)
    
##############################################################################
# Find overlap of annotation with peaks.

dictAnnotOL = {}


print astrBeds
for strAnnot in astrBeds:
    strAllAnnotHits = dirTmp + os.sep + os.path.basename(strAnnot) +"_AllHits.bed"   
    dfBed = pd.read_table(strAnnot,sep="\t",usecols=[0,1,2])
    dfBed.columns = ["chr","start","end"]
    with open(strAllAnnotHits,"w") as fileHits:
        for strPeaks in astrPeaks:
            strCutPeaks = dirTmp + os.sep + os.path.basename(strPeaks).replace(".bed","cut.bed")
            odb.IntersectBED( strAnnot,strCutPeaks,fileHits,args.cmdBedtools,strPeaks,"peaks against bed")
    strStem = os.path.basename(strAnnot)
    dfAnnotOL = odb.ConvertOverlapHits(strAllAnnotHits,dirTmp,dfBed,strStem)
    os.remove(strAllAnnotHits)
    dictAnnotOL[strAnnot] = dfAnnotOL

for strAnnot in astrGFFs:
    strAllAnnotHits = dirTmp + os.sep + os.path.basename(strAnnot) +"_AllHits.bed"   
    dfGFF = pd.read_table(strAnnot,sep="\t",usecols=[0,3,4])
    dfGFF.columns = ["chr","start","end"]
    with open(strAllAnnotHits,"w") as fileHits:
        for strPeaks in astrPeaks:
            strCutPeaks = dirTmp + os.sep + os.path.basename(strPeaks).replace(".bed","cut.bed")
            odb.IntersectGFF( strAnnot,strCutPeaks,fileHits,args.cmdBedtools,strPeaks,"peaks against gff")
    strStem = os.path.basename(strAnnot)
    dfAnnotOL = odb.ConvertOverlapHits(strAllAnnotHits,dirTmp,dfGFF,strStem)
    os.remove(strAllAnnotHits)
    dictAnnotOL[strAnnot] = dfAnnotOL





# At the end, each dfAnnotOverlap is a pandas dataframe with each row as an interval,
# and chr,start, and end as the first three columns, and each peak with 
# overlap as additional columns. Each cell contains the bp of overlap with the
# bed/gff interval.


###################################################################
# Get effective genome size for the binomial tests.
# This section merges all GFF and BED files used, and sums their
# overlaps to get the effective genome size.


astrBedsForN = []
for strBed in astrBeds:
    if "repeats.bed" not in strBed:
        astrBedsForN.append(strBed)


strUnionBed = args.dirOut + os.sep + "UnsortedEffectiveGenome.bed"
strBedEffGenome = args.dirOut + os.sep + "EffectiveGenome.bed"


# Cat all bed files, append to strBedEffGenome
with open(strUnionBed,'w') as fileUnionBed: 

    pCat = sp.Popen(["cat"] +astrBedsForN, stdout=sp.PIPE)
    stdoutCat = pCat.communicate()[0]    
    pCut = sp.Popen(["cut","-f1-3"],stdout=fileUnionBed,stdin=sp.PIPE)
    pCut.communicate(stdoutCat)
    
        
# Run merge on gff files to get them into bed format, append to strUnionBed
with open(strUnionBed,'a') as fileUnionBed:
    for strGFF in astrGFFs:
        pBedtools = sp.Popen([args.cmdBedtools,"merge", "-i",strGFF],
                            stdout=fileUnionBed)
        pBedtools.communicate()
    
# Sort the big combined bed file
strUnionSortedBed = args.dirOut + os.sep + "SortedEffectiveGenome.bed"
with open(strUnionSortedBed,'w') as fileUnionSortedBed: 
    sp.check_call([args.cmdBedtools,"sort","-i",strUnionBed],
                    stdout=fileUnionSortedBed)

# Merge the intervals in the sorted union bed file
odb.MergeIntervals( strUnionSortedBed,strBedEffGenome,args.cmdBedtools)


dfEffectiveGenome = pd.read_table(strBedEffGenome,sep="\t")
dfEffectiveGenome.columns = ["chr","start","end"]
iSizeEffGenome = (dfEffectiveGenome['end'] - dfEffectiveGenome['start']).sum()

print "Effective Genome Size:",iSizeEffGenome

dictIntervalsSize = {}
dictAnnotSize ={}
for strAnnot in astrBeds:
    dictAnnotSize[strAnnot] = odb.GetSizeBed(strAnnot)
    dictIntervalsSize[strAnnot] = odb.GetLinesFile(strAnnot)
for strAnnot in astrGFFs:
    dictAnnotSize[strAnnot] = odb.GetSizeGFF(strAnnot)
    dictIntervalsSize[strAnnot] = odb.GetLinesFile(strAnnot)
        



###################################################################
# Step 5: Run the binomial tests.
###################################################################

astrAnnot = astrBeds + astrGFFs


for strPeaks in dictPeakdfs.keys():
    print "Processing", strCutPeaks
    dfPV = pd.DataFrame(0, index=[], columns=['Annotation','Fold Enrichment','p-value',
    'Peaks with Overlap','Overlap(bp)','Total Peaks','Pct Peaks Overlapping','AnnotationPeaks','Annotation Peaks Overlapping','AnnotationSize(bp)',
    'AnnotationSize / EffectiveGenome','Effective Genome Size'])
    dfResults =dictPeakdfs[strPeaks]
    i=0
    for strAnnot in astrAnnot:
        print "Calculating results for", strAnnot
        
        if "repeats.bed" in strAnnot or "multiz30way_score_over0-70.bed" in strAnnot :
            iGSize = 2716965481 
        else: 
            iGSize = iSizeEffGenome

        if strAnnot not in (dfResults.columns.values):
            dfResults[strAnnot] = 0
            
            iOverlappingPeaks = (dfResults[strAnnot] >= 1).sum()
            iOverlapBP = dfResults[strAnnot].sum()
            iTotalPeaks = len(dfResults.index)
            dPctPeaksOverlap = np.nan
            dAnnotToGenomeRatio = dictAnnotSize[strAnnot]/float(iGSize)
            dEnrichment = np.nan
            
            dfOLData = dictAnnotOL[strAnnot]
            iAnnotPeaksWithOL = 0
            dPVal = np.nan
            
        else:
        
            iOverlappingPeaks = (dfResults[strAnnot] >= 1).sum()
            iOverlapBP = dfResults[strAnnot].sum()
            iTotalPeaks = len(dfResults.index)
            dPctPeaksOverlap = iOverlappingPeaks/float(iTotalPeaks)
            dAnnotToGenomeRatio = dictAnnotSize[strAnnot]/float(iGSize)
            dEnrichment = dPctPeaksOverlap/dAnnotToGenomeRatio            
            
            dfOLData = dictAnnotOL[strAnnot]
            iAnnotPeaksWithOL = (dfOLData[strPeaks] >= 1).sum() 
        
            dPVal = stats.binom_test(n = iTotalPeaks,x=iOverlappingPeaks, p=dAnnotToGenomeRatio)
        # Append rows to dataframe here.
        dfPV.loc[i] = [strAnnot,dEnrichment,dPVal,iOverlappingPeaks,iOverlapBP,iTotalPeaks,dPctPeaksOverlap,dictIntervalsSize[strAnnot],iAnnotPeaksWithOL,dictAnnotSize[strAnnot],dAnnotToGenomeRatio,iGSize]
        i+=1    
                        
            
            #print iAnnotOverlap,iTotalPeakLen,dAnnotToGenomeRatio     
            #print dictAnnotSize[strAnnot],iSizeEffGenome
            #print stats.binom_test(n = iTotalPeakLen,x=iAnnotOverlap, p=dAnnotToGenomeRatio)
    iTotalCoverage = (dfResults["end"]-dfResults["start"]).sum()
    dfPV['Annotation'] =  dfPV['Annotation'].map(lambda x: os.path.basename(x))
    dfPV.to_csv(dirOut+os.sep+os.path.basename(strPeaks+"_Statistics.tab"),sep="\t",index=False)
    with open(dirOut+os.sep+os.path.basename(strPeaks+"_Statistics.tab"),"a") as fileStats:
        fileStats.write("\n\n#Total bp covered by ATAC-Seq:\t" + str(iTotalCoverage)+ "\n")
    
    for strCol in list(dfResults.columns.values)[3:]:
        dfResults.loc[dfResults[strCol] > 0, strCol] = 1
    dfResults.columns = [os.path.basename(x) for x in list(dfResults.columns.values)]
    print "Finished calculating statistics. Will now merge on information on nearby genomic features."
    
    """
    for i in range(len(list(dfResults.columns.values)[3:])):
        strAnnot = list(dfResults.columns.values)[3:][i]
        print "Calculating results for", strAnnot
        iTotalCoverage = (dfResults["end"]-dfResults["start"]).sum()
        
        if strAnnot not in ["chr","start","end"] and strAnnot != ".":
            if "repeats.bed" in strAnnot or "multiz30way_score_over0-70.bed" in strAnnot :
                iGSize = 2716965481 
            else: 
                iGSize = iSizeEffGenome
            iOverlappingPeaks = (dfResults[strAnnot] >= 1).sum()
            iOverlapBP = dfResults[strAnnot].sum()
            iTotalPeaks = len(dfResults.index)
            dPctPeaksOverlap = iOverlappingPeaks/float(iTotalPeaks)
            dAnnotToGenomeRatio = dictAnnotSize[strAnnot]/float(iGSize)
            dEnrichment = dPctPeaksOverlap/dAnnotToGenomeRatio            
            
            dfOLData = dictAnnotOL[strAnnot]
            iAnnotPeaksWithOL = (dfOLData[strPeaks] >= 1).sum() 
            
            dPVal = stats.binom_test(n = iTotalPeaks,x=iOverlappingPeaks, p=dAnnotToGenomeRatio)
            # Append rows to dataframe here.
            dfPV.loc[i] = [strAnnot,dEnrichment,dPVal,iOverlappingPeaks,iOverlapBP,iTotalPeaks,dPctPeaksOverlap,dictIntervalsSize[strAnnot],iAnnotPeaksWithOL,dictAnnotSize[strAnnot],dAnnotToGenomeRatio,iGSize]
            
                        
            
            #print iAnnotOverlap,iTotalPeakLen,dAnnotToGenomeRatio     
            #print dictAnnotSize[strAnnot],iSizeEffGenome
            #print stats.binom_test(n = iTotalPeakLen,x=iAnnotOverlap, p=dAnnotToGenomeRatio)
    
    dfPV['Annotation'] =  dfPV['Annotation'].map(lambda x: os.path.basename(x))
    dfPV.to_csv(dirOut+os.sep+os.path.basename(strPeaks+"_Statistics.tab"),sep="\t",index=False)
    with open(dirOut+os.sep+os.path.basename(strPeaks+"_Statistics.tab"),"a") as fileStats:
        fileStats.write("\n\n#Total bp covered by ATAC-Seq:\t" + str(iTotalCoverage)+ "\n")
    
    for strCol in list(dfResults.columns.values)[3:]:
        dfResults.loc[dfResults[strCol] > 0, strCol] = 1
    dfResults.columns = [os.path.basename(x) for x in list(dfResults.columns.values)]
    """
###############################################################################
#  Step 6: Add genes

    
    dfGeneResults = odb.GetGeneTable(dictCutPeaks[strPeaks],"/data/yosef/index_files/mm9/genes/tight_genes.gtf",args.cmdBedtools,dirTmp,"Gene")    
    dfGenePromoter = odb.GetGeneTable(dictCutPeaks[strPeaks],"/data/yosef/index_files/mm9/genes/tight_genes.3p.gtf",args.cmdBedtools,dirTmp,"3prime of Gene")
    dfGenesThreePrime = odb.GetGeneTable(dictCutPeaks[strPeaks],"/data/yosef/index_files/mm9/genes/tight_genes.prom.gtf",args.cmdBedtools,dirTmp,"Promoter")
    dfClosestGeneResults = odb.GetClosestGTFFeature(dictCutPeaks[strPeaks],"/data/yosef/index_files/mm9/genes/tight_genes.gtf",args.cmdBedtools,dirTmp,"ClosestGene")
    
    dfResults = pd.merge(dfResults,dfGeneResults,how="left",on=["chr","start","end"])    
    dfResults = pd.merge(dfResults,dfGenePromoter,how="left",on=["chr","start","end"]) 
    dfResults = pd.merge(dfResults,dfGenesThreePrime,how="left",on=["chr","start","end"]) 
    dfResults = pd.merge(dfResults,dfClosestGeneResults,how="left",on=["chr","start","end"])     
    
    dfResults.to_csv(dirOut+os.sep+os.path.basename(strPeaks+"_Overlaps.tab"),sep="\t",index=False)
    #dfGeneResults.to_csv(dirOut+os.sep+os.path.basename(astrPeaks[iDF]+"_Genes.tab"),sep="\t",index=False)


# Clean up at end of run
os.remove(strUnionBed)
os.remove(strBedEffGenome)
os.remove(strUnionSortedBed) 
for strCutPeaks in dictCutPeaks.values():
    os.remove(strCutPeaks)
    
"""
if ".bed" in strAnnot:
    print odb.GetSizeBed(strAnnot)
if re.search(r'(.gff)|(.gtf)',strAnnot):
    print odb.GetSizeGFF(strAnnot)
#strSorted = dirOut+os.sep+os.path.basename(strBed)+"_sorted.bed"
#strMerged = dirOut+os.sep+os.path.basename(strBed)+"_merged.bed"
#odb.SortBed(strBed,strSorted,args.cmdBedtools)            
#odb.MergeIntervals(strSorted,strMerged,args.cmdBedtools)
"""
"""
# Think about if this should be based on length of peaks
iAnnotCov = odb.GetSizeBed(strMerged)
p = float(iAnnotCov)/float(iSizeEffGenome)
N = (dfPeaks["end"] - dfPeaks["start"]).sum()
x = dfPeaks[strBed].sum()

print strBed
print N,x,p
print stats.binom_test(x,N,p)
"""








"""
# Cut all bed files to first 3 columns.
astrBeds = odb.cut_files(astrBeds,dirTmp)
astrPeaks = odb.cut_files(astrPeaks,dirTmp)



adfPeaks = []
aiTotalPeakLength = []
for strPeaks in astrPeaks:
    
    strResults = args.dirOut + os.sep + os.path.basename(strPeaks)+".tab"
    strCutResults = args.dirOut + os.sep + os.path.basename(strPeaks)+"_cut.tab"
    fileOut = open(strResults,"w")        
    sp.check_call([args.cmdBedtools,"intersect", "-filenames","-wao","-a",strPeaks,
                   "-b"] + astrBeds,
                    stdout=fileOut)
    fileOut.close()
    
    fileCut = open(strCutResults,"w")
    sp.check_call(["cut", "-f1-4,8" ,strResults],stdout = fileCut)
    fileCut.close()
    
    dfMerged = pd.read_table(strCutResults,sep="\t",
                             names=['chr', 'start','end','bed_file','bp_overlap'],
                             engine="python")
    #dfMerged['Indicator'] = np.array([1]*len(dfMerged.index))
    #dfMerged.columns = ['chr', 'start','end','bed_file','bp_overlap']
    #dfMerged = dfMerged[dfMerged.bed_file != "."]
    dfMerged.to_csv("Data.tab")    
    dfMerged = pd.DataFrame(pd.pivot_table(dfMerged,index=['chr','start','end'],columns=['bed_file'],values='bp_overlap'))
    dfMerged.to_csv(dirOut + os.sep+ "FullResults.tab")
    dfMerged = pd.read_csv(dirOut + os.sep+ "FullResults.tab",sep=",",header=0)
    
    
    adfPeaks.append(dfMerged)
"""
    
    
