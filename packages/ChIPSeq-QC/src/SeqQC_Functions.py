#!/usr/bin/env python

import re
import os
import csv
def GetStem(strFile,strExt):
    reFind = re.compile('(^.*)\.'+ strExt + '$')
    mtchName = re.search(reFind,strFile)    
    return mtchName.group(1)
    
def check_create_dir( strDir ):
    if not os.path.exists( strDir ):
        os.makedirs( strDir )
    return strDir
    
def GetEffectiveGenome (strGenome):
    if strGenome[0:2] == "mm":
        return "mm"
    elif strGenome[0:2] == "hg":
        return "hs"

def AddStrings(strOne,strTwo):
    return str(float(strOne)+float(strTwo))

def ParseFlagStatOutput (strSAMOutput,strFormattedOut,strFileName):
    aaResults = []    
    
    with open(strSAMOutput, 'rb') as fileStats:
        readerStats = csv.reader(fileStats, delimiter=' ')
        for astrLine in readerStats:
            aaResults.append(astrLine)
        zipped = zip(*aaResults)
    
    
    
    """
    Sample SamOutput, 
    At this point, "zipped" only has the two columns of numbers.
    [0] 77724046 + 1192302 in total (QC-passed reads + QC-failed reads)
    [1] 0 + 0 duplicates
    [2] 54400384 + 547488 mapped (69.99%:45.92%)
    [3] 77724046 + 1192302 paired in sequencing
    [4] 38754812 + 599318 read1
    [5] 38969234 + 592984 read2
    [6] 52677948 + 472184 properly paired (67.78%:39.60%)
    [7] 53709634 + 495244 with itself and mate mapped
    [8] 690750 + 52244 singletons (0.89%:4.38%)
    [9] 885166 + 21032 with mate mapped to a different chr
    [10dMa 130100 + 6137 with mate mapped to a different chr (mapQ>=5)
    """
    astrQC = zipped[0]
    astrNotQC = zipped[2]
 
    with open(strFormattedOut,'w') as fileOut:
        fileOut.write('\t'.join([" ",strFileName,"\n"]))
        
        iIndex=0
        fileOut.write('\t'.join(["NumberOfReads",AddStrings(astrQC[iIndex],astrNotQC[iIndex]),"\n"]))
        
        iIndex=2
        fileOut.write('\t'.join(["AlignedReads",AddStrings(astrQC[iIndex],astrNotQC[iIndex]),"\n"]))
        
        #dMapped_PassedQC = int(astrQC[2])/float(astrQC[0])
        #dMapped_FailedQC = int(astrNotQC[2])/float(astrNotQC[0])
        if (float(astrNotQC[0])+float(astrQC[0])) >0:
            dAllMapped = (int(astrNotQC[2])+int(astrQC[2]))/(float(astrNotQC[0])+float(astrQC[0]))
            fileOut.write('\t'.join(["PercentAligned",str(dAllMapped),"\n"]))
        else:
            fileOut.write('\t'.join(["PercentAligned","NA","\n"]))
        
        iIndex=1        
        #fileOut.write('\t'.join(["TotalDuplicates",AddStrings(astrQC[iIndex],astrNotQC[iIndex]),"\n"]))
    return
    
def CheckPrimers (astrPrimerFiles,strFormattedOut):
    iPrimerHits = 0
    iTotalSeqs = 0    
    
    for strPrimerOut in astrPrimerFiles:
        with open(strPrimerOut, 'rb') as filePrimer:
            astrData = filePrimer.readline().split('\t')
            iPrimerHits += int(astrData[1])
            iTotalSeqs += int(astrData[2])
    
    dPctPrimer = float(iPrimerHits/iTotalSeqs)    
    
    with open(strFormattedOut,'a') as fileOut:
        fileOut.write('\t'.join(["PercentPrimerHits",str(dPctPrimer),"\n"]))
        
def GetAlnData (strAlnFile,strFormattedOut,iLines=4,iKeep=3):
    iGetLines = 0
    aaData = []    
    with open(strAlnFile, 'rb') as fileAln:
            for strLine in fileAln:
                mtchName = re.search("METRICS CLASS",strLine)
                if iGetLines >0:
                    aaData.append(strLine.split('\t'))
                    iGetLines = iGetLines-1
                if mtchName:
                    #print strLine
                    iGetLines =iLines

    aaTranspose = zip(*aaData)
    
    with open(strFormattedOut,'a') as fileOut:
        for tupData in aaTranspose:
            fileOut.write('\t'.join([tupData[0].strip(),tupData[iKeep].strip(),"\n"]))    
    
    return
    


def GetComplexity(strDupTxt,strFormattedOut):
    astrData = open(strDupTxt).readline().split('\t')
    if astrData[1]!="-1":
        dComplexity = 1 - float(astrData[1])
        with open(strFormattedOut,"a") as fOut:
            fOut.write("Complexity"+"\t" + str(dComplexity))
        return
    else:
        with open(strFormattedOut,"a") as fOut:
            fOut.write("Complexity"+"\t" + "Error: -1 in dup.txt") 
        return

             
"""
strTest =     "/data/yosef/CD8_effector_diff/tmp/test_prealn/formatted.tab"
ParseFlagStatOutput("/data/yosef/CD8_effector_diff/tmp/test_prealn/flagstat_test.txt",
strTest)

CheckPrimers(["/data/yosef/CD8_effector_diff/tmp/test_prealn/fastqc/primer.1.txt",
"/data/yosef/CD8_effector_diff/tmp/test_prealn/fastqc/primer.2.txt"],
strTest)
GetAlnData ("/data/yosef/CD8_effector_diff/tmp/test_prealn/fastqc/aln_metrics.txt",
strTest)
GetComplexity("/data/yosef/CD8_effector_diff/tmp/test_prealn/fastqc/dup.txt",
strTest)
"""