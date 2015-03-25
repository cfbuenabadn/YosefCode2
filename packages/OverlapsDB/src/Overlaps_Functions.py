##############################################################################
# Jim Kaminski
# Andre King assisted with the config_file option
# (based on code originally written by Faraz Tavakoli)
# 2/27/2015
# Yosef Lab
##############################################################################
"""
This program is designed to take a file of peaks from MACS2
and several bed files of annotations on the genome. It identifies overlap
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
import os
import subprocess as sp
import pandas as pd
import re
import numpy as np

def cut_file(strIn,strOut):
    # Gets first 3 columns of tab-delimited file. This is used to cut out
    # down bed files to (chr,start,end) and avoid extra material.
    
    with open(strOut,'w') as fileCutBed:
        pCut = sp.Popen(["cut","-f1-3",strIn],stdout=sp.PIPE)
        stdoutCut = pCut.communicate()[0]
        pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,
                            stdout=fileCutBed)
        pSort.communicate(stdoutCut)[0]
        
    return strOut


def check_create_dir( strDir ):
    if not os.path.exists( strDir ):
        os.makedirs( strDir )
    return strDir
    
def IntersectBED( strBed, strCutPeaks,fileHits,cmdBedtools,strField,strOperation):
    # This function gets the intersection of your peaks file (strCutPeaks)
    # and your bed file (strBed), and writes the rows of overlap to fileHits.
    
    strBedName = os.path.basename(strBed).strip()
    strSedFN = r"s/$/\tNAME/".replace("NAME",re.escape(strField))
    
    print "Processing",strBedName,"..."
    print "Cutting",strBedName,"..."
    pCut = sp.Popen(["cut","-f1-3",strBed],stdout=sp.PIPE)
    stdoutCut = pCut.communicate()[0]


    
    print "Sorting", strBedName,"..."
    pSort = sp.Popen(["sort","-k","1,1","-k2,2n"],stdin=sp.PIPE,
                            stdout=sp.PIPE)
    stdOutSort = pSort.communicate(stdoutCut)[0]    
    
    
    print "Running bedtools intersect on",strBedName,"..."
    
    if  strOperation == "bed against peaks":    
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-wo","-sorted","-a",strCutPeaks,
                               "-b","stdin"],stdin=sp.PIPE,
                                stdout=sp.PIPE)
    elif  strOperation == "peaks against bed":    
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-wo","-sorted","-a","stdin",
                               "-b",strCutPeaks],stdin=sp.PIPE,
                                stdout=sp.PIPE)
    stdoutBedtools = pBedtools.communicate(stdOutSort)[0]
    
    print "Appending to hits file...",strBed,"..."
    pSed = sp.Popen(["sed","-e",strSedFN],stdin=sp.PIPE,
                    stdout=fileHits)
    pSed.communicate(stdoutBedtools)
    
    return 
    

    
def IntersectGFF( strGFF, strCutPeaks,fileHits,cmdBedtools,strField,strOperation):
    # This function gets the intersection of your peaks file (strCutPeaks)
    # and your GFF file (strGFF), and writes the rows of overlap to fileHits.
    
    strGFFName = os.path.basename(strGFF).strip()
    strSedFN = r"s/$/\tNAME/".replace("NAME",re.escape(strField))
    print "Processing",strGFFName,"..."

    
    print "Running bedtools intersect on",strGFFName,"(GFF file) ..."
    if  strOperation == "gff against peaks":
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-wo","-bed","-a",strCutPeaks,
                           "-b",strGFF],stdout=sp.PIPE)
        stdoutBedtools = pBedtools.communicate()[0]
        pCut = sp.Popen(["cut","-f1,2,3,4,7,8,13"],stdout=sp.PIPE,stdin=sp.PIPE)
        outCut = pCut.communicate(stdoutBedtools)[0]
        print "Appending to hits file...",strGFFName,"..."
        pSed = sp.Popen(["sed","-e",strSedFN],stdin=sp.PIPE,
                        stdout=fileHits)
        pSed.communicate(outCut)
    elif  strOperation == "peaks against gff":
        pCut = sp.Popen(["cut","-f","1,4,5",strGFF],stdout=sp.PIPE)
        print "Started cut..."
        outCut = pCut.communicate()[0]        
        print "Finished cut..."
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-wo","-bed","-a","stdin",
                           "-b",strCutPeaks],stdin=sp.PIPE,stdout=sp.PIPE)
        stdoutBedtools = pBedtools.communicate(outCut)[0]
        
        print "Appending to hits file...",strGFFName,"..."
        pSed = sp.Popen(["sed","-e",strSedFN],stdin=sp.PIPE,
                        stdout=fileHits)
        pSed.communicate(stdoutBedtools)
    
    return 

def MergeIntervals( strIn,strOut,cmdBedtools):    
    fileOut = open(strOut,'w')
    sp.check_call([cmdBedtools,"merge", "-i",strIn],
                        stdout=fileOut)
    fileOut.close()

def GetLinesFile(strFile):
    with open(strFile) as fileData:
        iTotal =sum(1 for strLine in fileData)
    return iTotal

def GetSizeBed(strBed):
    
    dfBed = pd.read_table(strBed,sep="\t",usecols=[0,1,2])
    dfBed.columns = ["chr","start","end"]
    iSize = (dfBed['end'] - dfBed['start']).sum()
    return iSize
def GetSizeGFF(strGFF):
    
    dfGFF = pd.read_table(strGFF,sep="\t",usecols=[3,4])
    dfGFF.columns = ["start","end"]
    iSize = (dfGFF['end'] - dfGFF['start']).sum()
    return iSize
   
def SortBed(strIn,strOut,cmdBedtools):
    with open(strOut,'w') as fileSortedBed: 
        sp.check_call([cmdBedtools,"sort","-i",strIn],
                        stdout=fileSortedBed)
                        
def GetGeneTable(strPeaks,strGTF,cmdBedtools,dirTmp,strTitle):
    strTmpResults = dirTmp + os.sep + "GeneTable.bed"
        
    
    with open(strTmpResults,'w') as fileTmpResults:
        pSed = sp.Popen(["sed",'/chrNT/d',strGTF],stdout=sp.PIPE)
        outSed = pSed.communicate()[0]
        pBedtools = sp.Popen([cmdBedtools,"intersect", "-wb","-bed","-a",strPeaks,
                          "-b","stdin"],stdout=fileTmpResults,stdin=sp.PIPE)
        pBedtools.communicate(outSed)
        
    dfGenes = pd.read_table(strTmpResults,sep="\t",usecols=[0,1,2,8],names=["chr","start","end",strTitle])
    dfGenes = dfGenes.drop_duplicates()
    dfGenes[strTitle] = dfGenes.groupby(['chr','start','end'])[strTitle].transform(lambda x: ','.join(x))
    dfGenes = dfGenes.drop_duplicates()
    return dfGenes
    
    
def ReadMatlabDB(strDir,dirTmp):
    """ Based on Faraz's Code"""
    # mlab_col_bin.txt - Dict of Gene Names - 
    # mlab_dat_bin.txt - Gene Columns, each column has an ID, then indicates 0,1,2 for region
    # Chr,start,end - mlab_gname.txt

    # Get Gene column names.    
    astrBeds = []
    dictGeneNames = {}
    with open(strDir+os.sep+"mlab_col_bin.txt",'r') as fileColNames:
        for strLine in fileColNames:
            astrLine = strLine.split("\t")
            dictGeneNames[astrLine[1].strip()] = astrLine[0].strip()
            
    
            
    # Load in data
    dfIntervals = pd.read_csv(os.path.join(strDir, "mlab_gname.txt"), delim_whitespace = True,header=None)
    dfIntervals.columns = ["chr","start","end","crm","gene"]
    dfGeneColumns = pd.read_csv(os.path.join(strDir, "mlab_dat_bin.txt"),delim_whitespace = True,header=0)

    # Rename the columns
    dfCombined = dfIntervals.join(dfGeneColumns)
    dfCombined.columns.values[5:] = [dictGeneNames[x] for x in dfCombined.columns.values[5:]]
    dfCombined.to_csv(dirTmp+os.sep+"Test.tab")
    
    for strAnnot in dfCombined.columns.values[5:]:
        strBedName = dirTmp+os.sep+strAnnot+".bed"
        dfBed = dfCombined[dfCombined[strAnnot] != 0]
        dfBed = dfBed[["chr","start","end"]]
        dfBed.to_csv(strBedName,sep="\t",header=False,index=False)
        astrBeds.append(strBedName)
    
    return astrBeds
                                
    # Concatenate dfIntervals and dfGeneColumns.
    # Then loop through, filter, and export to bed.
def ConvertOverlapHits(strAllHits,dirTmp,dfPeaks,strStem):
    dfAllHits = pd.read_table(strAllHits,sep="\t",usecols=[0,1,2,6,7],names=["chr","start","end","bp_overlap","bed"],
                              dtype={'chr': object, 'start': np.int64,'end':np.int64,'bp_overlap':np.int64,'bed':object} )
    iRows =  dfAllHits.shape[0]
    print "Total rows in dataframe:",iRows
    dfAllHits.to_csv(dirTmp+os.sep+strStem+"_HitTable.tab",sep="\t")
    
    strReshaped = dirTmp+os.sep+strStem+"_ReshapedHits.tab"
    strFullTable = dirTmp+os.sep+strStem+"_FullTable.tab"
    if iRows>0:
    # Reshape, and write to dirTmp+os.sep+"RehshapedHits.csv"
        pd.pivot_table(dfAllHits, index=['chr','start','end'],columns=['bed'],values='bp_overlap').to_csv(strReshaped,sep="\t")
        dfFullTable = pd.read_table(strReshaped,sep="\t",header=0)
        
    else:
        # May want to put in the first peak of strPeaks here, and set overlap equal to 0.
        dfFullTable =  pd.DataFrame(0, index=[], columns=['chr','start','end'])
        
    dfFullTable = pd.merge(dfPeaks, dfFullTable, how='left', on=["chr","start","end"])
    dfFullTable.to_csv(strFullTable,sep="\t")
        
    return dfFullTable
    
    
    
    
    
    