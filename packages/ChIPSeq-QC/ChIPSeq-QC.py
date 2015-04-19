#!/usr/bin/env python


import argparse
import subprocess as sp
import os
import re
import sys
import src.SeqQC_Functions as SeqQC

c_strStringency = "STRICT"

# Jim Kaminski (editing Nir Yosef's code)
parser = argparse.ArgumentParser(description='')



grpParam = parser.add_argument_group('Input: \n(1) Please enter either: \n * 1 or 2 FASTQ files\
 \n * A BAM file with unaligned reads  \
 \n * A BAM file with aligned reads \
 \n Also, supply (2) a directory containing index files for the reference genome, and (3) add --paired if your reads are paired.' )
grpParam.add_argument('--fastq_1', type=str, dest='strFASTQ_1',default="", help='FASTQ File One')
grpParam.add_argument('--fastq_2', type=str, dest='strFASTQ_2',default="", help='FASTQ File Two (if using paired end data.)')
grpParam.add_argument('--bam', type=str, dest='strBAM', default="",help='BAM file. (Pipeline will convert to FASTQ. Please specify if paired or single end.)')
grpParam.add_argument('--aligned_bam', type=str, dest='strPreAlignedBAM', default="",help='Already aligned BAM file. (Pipeline will skip bowtie/bowtie2 and move to Step3)')
grpParam.add_argument('--paired', dest='fPaired', action='store_true',help='Please add this argument if your data are paired read ends')
grpParam.add_argument('--refgenome_index', type=str, dest='dirRefGenome', help='Please enter the folder where you have stored your reference genome')

grpParam = parser.add_argument_group('Output Directories - Please provide (4) --out, which will tell the pipeline where to output the results.')
grpParam.add_argument('--out', type=str, dest='dirOut', help='Output folder',default=os.getcwd()+os.sep+"out")
grpParam.add_argument('--tmp', type=str, dest='dirTmp', help='Tmp folder - used for Picard tools',default="")

grpParam = parser.add_argument_group('Analysis Options')
grpParam.add_argument('--threads', dest='iThreads', default=1,help='Indicate number of threads for pipeline to use. As of 2-24-2015 this only affects bowtie.')
grpParam.add_argument('--runmacs', dest='fRunMACS', action='store_true',help='Please add this argument if you want to run macs2. This is performed without a control.')
grpParam.add_argument('--remove_dup', dest='fRemoveDup', action='store_true',help='By default, we mark duplicates using Picard. If you wish to remove duplicates, please add this flag.')
parser.set_defaults(fRunMACS=False)
parser.set_defaults(fRemoveDup=False)


grpParam = parser.add_argument_group('If you are not running this on the YosefLab system, please enter full paths to the \
bioinformatics programs and directories below.')
grpParam.add_argument('--picard_dir', type=str, dest='dirPicard', default="/opt/pkg/picard-tools-1.108/lib", help='Please enter the full path for where you have stored Picard tools, default is /opt/pkg/picard-tools-1.108/lib')
grpParam.add_argument('--samtools', type=str, dest='cmdSamtools', default="/opt/genomics/bin/samtools", help='Please enter the full path for where you have stored samtools, default is /opt/genomics/bin/samtools')
grpParam.add_argument('--bowtie', type=str, dest='cmdBowtie', default="/data/yosef/CD8_effector_diff/programs/bowtie-1.1.1/bowtie", help='Please enter the full path for where you have stored bowtie, default is /data/yosef/CD8_effector_diff/programs/bowtie-1.1.1/bowtie')
grpParam.add_argument('--igvtools', type=str, dest='cmdIGV', default="/opt/pkg/IGVTools/igvtools", help='Please enter the full path for where you have stored igvtools, default is /opt/pkg/IGVTools/igvtools')
grpParam.add_argument('--fastqc', type=str, dest='cmdFASTQC', default="/opt/genomics/bin/fastqc", help='Please enter the full path for where you have stored fastqc, default is to call \"fastqc\" ')
grpParam.add_argument('--macs2', type=str, dest='cmdMACS2', default="macs2", help='Please enter the full path for where you have stored macs2, default is to call \"macs2\" ')
parser.set_defaults(fPaired=False)
parser.set_defaults(fRunMACS=False)
args = parser.parse_args()


# Choose index folder for mouse genome.
#c_dirGenomeIndex = "/data/yosef/index_files/mm9/genome/mm9"

dirSrc = os.path.dirname(os.path.realpath(__file__))  + os.sep + "src"
sys.stderr.write("Testing src path: " + dirSrc)

if args.dirTmp =="":
    args.dirTmp = args.dirOut + os.sep + "tmp"
###############################################################################
# Step 0: If needed, convert BAM to FASTQ. We run FASTQC on these files later.
#   Input: Unaligned or aligned bam
#   Output: 1 or 2 FASTQ files in dirOUT, depending on whether or not BAM is flagges as --paired
    
fileLog = open(args.dirOut + os.sep + "ChipSeqQC.log","w")

SeqQC.check_create_dir(args.dirOut)
strGenome = os.path.split(args.dirRefGenome)[1]
print "Preparing to run QC on your reads and align them to the " + strGenome + " genome."
bSkipAln = False

if (args.strBAM!= "" and args.fPaired== False):
    print "Converting BAM to FASTQ for pipeline, treating it as single end data."
    strProjectStem = SeqQC.GetStem(os.path.basename(args.strBAM),"bam")
    strFASTQ_1 = args.dirOut + os.sep + strProjectStem + "_1.fastq"
    fPaired= False
    sp.check_call(["java","-jar","-Xmx2g", args.cmdSamToFASTQ,"I="+args.strBAM, "F="+strFASTQ_1, "VALIDATION_STRINGENCY="+c_strStringency],stderr=fileLog)

elif (args.strBAM!= "" and args.fPaired== True):
    print "Converting BAM to FASTQ for pipeline, treating it as paired end data."
    strProjectStem = SeqQC.GetStem(os.path.basename(args.strBAM),"bam")
    strFASTQ_1 = args.dirOut + os.sep + strProjectStem + "_1.fastq"
    strFASTQ_2 = args.dirOut + os.sep + strProjectStem + "_2.fastq"
    fPaired= True
   
    sp.check_call(["java","-jar","-Xmx2g", args.dirPicard+os.sep+"SamToFastq.jar","I="+args.strBAM, "F="+strFASTQ_1, "F2="+strFASTQ_2, "VALIDATION_STRINGENCY="+c_strStringency],stderr=fileLog)    
 
elif (args.strFASTQ_1!="" and args.strFASTQ_2==""):
    print "Working with one FASTQ file, treating it as single end data."
    fPaired = False
    strFASTQ_1=args.strFASTQ_1
    strProjectStem = SeqQC.GetStem(os.path.basename(args.strFASTQ_1),"fastq")
elif (args.strFASTQ_1!= "" and args.strFASTQ_2!=""):
    print "Working with two FASTQ files, treating them as paired end data."
    fPaired = True
    strFASTQ_1=args.strFASTQ_1
    strFASTQ_2=args.strFASTQ_2
    strProjectStem = SeqQC.GetStem(os.path.basename(args.strFASTQ_1),"fastq")
elif (args.strPreAlignedBAM!=""):
    print "Working with a BAM file of aligned reads."
    bSkipAln=True
    fPaired=args.fPaired
    strProjectStem = SeqQC.GetStem(os.path.basename(args.strPreAlignedBAM),"bam")
    strFASTQ_1 = args.dirOut + os.sep + strProjectStem + "_1.fastq"
    
    strAlignBamStats = args.dirOut + os.sep + "stats_for_initial_bowtie_alignment.txt"
    fAlignBamStats = open(strAlignBamStats,'w')
    sp.call([args.cmdSamtools,"flagstat",args.strPreAlignedBAM],stdout=fAlignBamStats)
    fAlignBamStats.close()    
    
    
    if fPaired == True:
        print "Converting aligned BAM to FASTQ for QC downstream, treating it as paired end data."
        strFASTQ_2 = args.dirOut + os.sep + strProjectStem + "_2.fastq"
        sp.check_call(["java","-jar","-Xmx2g", args.dirPicard+os.sep+"SamToFastq.jar","I="+args.strPreAlignedBAM, 
        "F="+strFASTQ_1, "F2="+strFASTQ_2,"VALIDATION_STRINGENCY="+c_strStringency],stderr=fileLog)    
    else:
        print "Converting aligned BAM to FASTQ for QC downstream, treating it as single end data."
        sp.check_call(["java","-jar","-Xmx2g", args.dirPicard+os.sep+"SamToFastq.jar","I="+args.strPreAlignedBAM,
                       "F="+strFASTQ_1,  "VALIDATION_STRINGENCY"+c_strStringency],stderr=fileLog)    
    
    
else:
    print "\n No FASTQ files or bam files provided. Please load the help by typing \n \
    \'python ATAC-Seq_ChIP-Seq_Pipeline.py --help'."
    
###############################################################################
# Step 1: Align to the reference genome using bowtie. 
#(Skip if starting with a prealigned bam file.)
    # Input: FASTQ file(s) 
    # Output: SAM file with alignment of reads to genome.
if bSkipAln==False:
    print "Step 1: Aligning file(s) to the reference genome."
    strSamInitAlignment =  args.dirOut + os.sep + "initial_bowtie_alignment.sam"
    
    
    fBowtieStdErr = open(args.dirOut+os.sep+"bowtie_stderr.txt",'w')
    
    if fPaired==True:
        sp.check_call([args.cmdBowtie,"--sam",args.dirRefGenome, "-p", str(args.iThreads), "-1", 
                       strFASTQ_1, "-2", strFASTQ_2, "--chunkmbs", "1500", "-X", "2500", "-m", "1", strSamInitAlignment],stderr=fBowtieStdErr)
    else:
        sp.check_call([args.cmdBowtie,"--sam",args.dirRefGenome, "-p", str(args.iThreads), "--chunkmbs", "1500", strFASTQ_1,  
                       strSamInitAlignment],stderr=fBowtieStdErr)
    fBowtieStdErr.close()
    # Note to Jim: Include support for bowtie2 here.
    
    ###############################################################################
    # Step 2:Convert SAM to BAM. 
        # Input: SAM file with alignment of reads to genome.[initial_bowtie_alignment.sam]
        # Output: BAM file with alignment of reads to genome, sorted, without duplicates. [aligned_sorted_noduplicates.bam]
    
    print "Step 2: Convert SAM to BAM."
    strAlignBam = args.dirOut + os.sep + "initial_bowtie_alignment.bam"
    fAlignBam = open(strAlignBam,'w')
    sp.call([args.cmdSamtools,"view","-b", "-S", strSamInitAlignment],stdout=fAlignBam,stderr=fileLog)
    fAlignBam.close()
    
    strAlignBamStats = args.dirOut + os.sep + "stats_for_initial_bowtie_alignment.txt"
    fAlignBamStats = open(strAlignBamStats,'w')
    sp.call([args.cmdSamtools,"flagstat",strAlignBam],stdout=fAlignBamStats,stderr=fileLog)
    fAlignBamStats.close()
    
    
else:
    strAlignBam = args.strPreAlignedBAM

###############################################################################
# Step 3: Sort BAM and Mark or Remove Duplicates.
#   Input: Aligned Bam File
#   Output: (1) a sorted bam file (later removed)
#           (2) a sorted bam file with duplicates marked,
#           (3) an index file for the sorted, duplicates marked file
#           (4) tdf files for IGV

print "Step 3 -Sort and Remove Duplicates"
print "Step 3 - Sorting..."
strSortedBam = args.dirOut + os.sep + "aligned_sorted.bam"
sp.check_call(["java","-Xmx1200m","-jar", args.dirPicard+os.sep+"SortSam.jar","I="+strAlignBam,
"O="+strSortedBam, "SO=coordinate", "TMP_DIR="+args.dirTmp],stderr=fileLog)



if args.fRemoveDup==True:
    print "Step 3 - Removing Duplicates...";
    strSortedNoDupBam = args.dirOut+os.sep + "aligned_sorted_no_duplicates.bam"
    sp.check_call(["java","-Xmx1000m", "-jar",args.dirPicard+os.sep+"MarkDuplicates.jar","REMOVE_DUPLICATES=true",
    "METRICS_FILE="+args.dirOut+os.sep+"duplicate_reads.txt","INPUT="+strSortedBam,
    "OUTPUT="+strSortedNoDupBam, "TMP_DIR="+args.dirTmp,"ASSUME_SORTED=true"],stderr=fileLog)
    strTDFfile =  args.dirOut + os.sep+"aligned_sorted_no_duplicates.tdf"
else:
    print "Step 3 - Marking Duplicates...";
    strSortedNoDupBam = args.dirOut+os.sep + "aligned_sorted_marked_duplicates.bam"
    sp.check_call(["java","-Xmx1000m", "-jar",args.dirPicard+os.sep+"MarkDuplicates.jar","REMOVE_DUPLICATES=false",
    "METRICS_FILE="+args.dirOut+os.sep+"duplicate_reads.txt","INPUT="+strSortedBam,
    "OUTPUT="+strSortedNoDupBam, "TMP_DIR="+args.dirTmp,"ASSUME_SORTED=true"],stderr=fileLog)
    strTDFfile =  args.dirOut + os.sep+"aligned_sorted_marked_duplicates.tdf"
    

print "Step 3 - Indexing...";
sp.check_call([args.cmdSamtools, "index", strSortedNoDupBam])

print "Step 3 - Making TDF files for IGV...";

if fPaired==True:
    sp.check_call([args.cmdIGV, "count","--pairs",strSortedNoDupBam,strTDFfile, strGenome],stderr=fileLog)
else:
    sp.check_call([args.cmdIGV, "count",strSortedNoDupBam,strTDFfile, strGenome],stderr=fileLog)
        

###############################################################################
# Step 4: Run FASTQC on the FASTQ reads
#   Input: FASTQ reads
#   Output: /out/fastqc data

dirFASTQC = args.dirOut + os.sep + "fastqc"
SeqQC.check_create_dir(dirFASTQC)
print "Step 4: Run FASTQC"


sp.call([args.cmdFASTQC, strFASTQ_1, "-o",  dirFASTQC])
if fPaired==True:
    sp.call([args.cmdFASTQC, strFASTQ_2, "-o",  dirFASTQC])


###############################################################################
# Step 5: Check for contamination in FASTQ files with Nir's check_contam.pl.
#   Input: FASTQ files
#   Output: primer.1.txt and primer.2.txt

print "Step 5: Check FASTQ files for contamination."

sp.call(["perl",dirSrc + os.sep + "check_contam.pl",strFASTQ_1,dirFASTQC+os.sep+"primer.1.txt"])
if fPaired==True:
   sp.call(["perl",dirSrc + os.sep + "check_contam.pl",strFASTQ_2,dirFASTQC+os.sep+"primer.2.txt"])

###############################################################################
# Step 6: Report Summary Information on the Sorted BAM File with Picard tools,
# and Nir's count_dup.pl.
#   Input: Sorted BAM file
#   Output: aln_metrics.txt, insert_metrics.txt, DUP.TXT
print "Step 6: Report Summary Information (alignment metrics, insert size metrics, duplicates)."


#system("java -Xmx1000m -jar /opt/pkg/picard-tools-1.108/lib/CollectAlignmentSummaryMetrics.jar INPUT=$OUTPUT_FOLDER/sorted.bam OUTPUT=$OUTPUT_FOLDER/Fastqc/aln_metrics.txt\n");
sp.call(["java","-Xmx1000m","-jar", args.dirPicard + os.sep + "CollectAlignmentSummaryMetrics.jar", "Input="+strSortedBam ,"OUTPUT="+dirFASTQC+os.sep+"aln_metrics.txt", "TMP_DIR="+args.dirTmp],stderr=fileLog)

#system("java -Xmx1000m -jar /opt/pkg/picard-tools-1.108/lib/CollectInsertSizeMetrics.jar INPUT=$OUTPUT_FOLDER/sorted.bam OUTPUT=$OUTPUT_FOLDER/Fastqc/aln_metrics1.txt H=$OUTPUT_FOLDER/Fastqc/ins_metrics.histogram.pdf\n");
sp.call(["java","-Xmx1000m","-jar", args.dirPicard + os.sep + "CollectInsertSizeMetrics.jar", "Input="+strSortedBam,"OUTPUT="+dirFASTQC+os.sep+"insert_metrics.txt", "H="+dirFASTQC+os.sep+"ins_metrics.histogram.pdf", "TMP_DIR="+args.dirTmp],stderr=fileLog)

#system("perl /data/yosef/CD8_effector_diff/scripts/preproc-jim/count_dup.pl $OUTPUT_FOLDER/sorted.bam $OUTPUT_FOLDER/Fastqc/dup.txt\n");
sp.call(["perl",dirSrc + os.sep + "count_dup.pl",strSortedNoDupBam,dirFASTQC+os.sep+"dup.txt"])

# add a Python script here like Nir's
# Have it go through the bam file, line by line, check the flags on each, nd count:
#Number of reads (/opt/genomics/bin/samtools view -c QC/aligned_sorted_marked_duplicates.bam)
#Number of reads that aligned (/opt/genomics/bin/samtools view -c QC/aligned_sorted_marked_duplicates.bam)
#Number of reads that aligned and are duplicates (/opt/genomics/bin/samtools view -c -f 1026 QC/aligned_sorted_marked_duplicates.bam)


strSummaryStats = args.dirOut+os.sep+"SummaryStats.tab"  
SeqQC.ParseFlagStatOutput(strAlignBamStats,strSummaryStats)

###############################################################################
# Step 7: Add the FASTQ check results to SummaryStats.tab


astrPrimerFiles = [dirFASTQC+os.sep+"primer.1.txt"]
if fPaired==True:
    astrPrimerFiles.append(dirFASTQC+os.sep+"primer.2.txt")


SeqQC.CheckPrimers(astrPrimerFiles,strSummaryStats)
SeqQC.GetAlnData (dirFASTQC+os.sep+"aln_metrics.txt",strSummaryStats)







###############################################################################
# Step 8: Can call peaks if desired.
#   Input: Sorted Bam file, with duplicates marked
#   Output - tab delimited file of peaks. (MACS adds xls extension, but it's 
#       really just a tab-delimited file.)

if args.fRunMACS == True:
    strEffectiveGenome = SeqQC.GetEffectiveGenome(strGenome)
    sp.check_call([args.cmdMACS2, "callpeak", "--nomodel", "-t", strSortedNoDupBam, 
    "-f","BAM","-n",args.dirOut+os.sep +"res","-g", strEffectiveGenome],stderr=fileLog)

###############################################################################
# Step 9: Clean up

#Remove sorted_bam file.
os.remove(strSortedBam)

#Remove FASTQ files if we did not start with FASTQ files.
if (args.strFASTQ_1=="" and args.strFASTQ_2==""):
    os.remove(strFASTQ_1)
    if os.path.isfile(strFASTQ_2):
        os.remove(strFASTQ_2)
        
fileLog.close()