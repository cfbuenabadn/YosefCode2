#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015

import subprocess
from string import Template
import argparse
import os;
import doQC;
import sys;
import shutil;
import rnaSeqPipelineUtils;
import cleanup_after_preproc;

parser = argparse.ArgumentParser(description="Process one single cell sample", parents=[rnaSeqPipelineUtils.common_rnaseq_parser])
parser.add_argument('sampleFile1', action='store',
                   help='the first reads file')
parser.add_argument('sampleFile2', action='store', nargs = '?', default='',
                   help='the second reads file (for paired-end. If the --paired_end flag is not specified and this argument is specified then an exception is thrown. If the --paired_end flag is raised but this file is not specified then again an exception is thrown')

args = parser.parse_args()

if(not(args.paired_end) and args.sampleFile2 != ''):
	raise Exception("you cannot specify a second read file unless you raise the --paired_end flag");
elif(args.paired_end and args.sampleFile2 == ''):
	raise Exception("if the --paired_end flag is raised, then you must specify a second reads file");

if(args.reference == "mm10"):
	#settings for mm10 with ERCC spike-ins and other additions required by the BRAIN
	#project (eGFP, tdTomato, CreER)

	KALLISTO_INDEX_FILE = "/data/yosef/index_files/mm10_4brain/index/kallisto_index/kallisto_index_mm10_4brain.idx"

	GENOME_REFERENCE_FILE = "/data/yosef/index_files/mm10_4brain/index/GCF_000001635.23_GRCm38.p3_genomic_4brain.fna";
	BOWTIE2_INDEX = "/data/yosef/index_files/mm10_4brain/index/GRCm38.p3_4brain";
	TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/mm10_4brain/index/tophat_transcriptome_data/GRCm38.p3_4brain_refseq_annot";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/mm10_4brain/index/combinedGFF_4brain.gff";
	RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/mm10_4brain/index/rsem_index/combinedGTF_4brain.gtf";
	RSEM_INDEX = "/data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem";
	#REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat/GRCm38.p3.refFlat";
	#REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat2/refFlat.txt";
	REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat2/refFlat_allonEdited.txt";
	#RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/gencode.vM4.rRNA.interval_list";
	#RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/gencode.vM4.rRNA.interval_list_allonEdited";
	RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/mm10_4BRAIN_rRNA_interval_list_allonEdited.txt";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "null" # Originally I produced this file, but when I updated the rsem dictionary I didn't think it was important enough to update also the ribosomal file. The path in which I had the original ribosomal reference:"/data/yosef/index_files/mm10_4brain/index/gencode/mm10_4BRAIN_rRNA_interval_list_allonEdited_forRSEM.txt";
	RIBOSOMAL_INTERVALS_INDEX_FOR_KALLISTO = "null" #I didn't create it - no time and not that important...
	RSEM_DICTIONARY = "/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";

	#settings for mm10 with ERCC spike-ins
	#BOWTIE2_INDEX = "/data/yosef/index_files/mm10_withERCC/GRCm38.p3_withERCC";
	#TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/mm10_withERCC/tophat_transcriptome_data/GRCm38.p3_refseq_annot";
	#TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/mm10_withERCC/GRCm38.p3.gff";
	#RSEM_INDEX = "/data/yosef/index_files/mm10_withERCC/rsem_index/GRCm38.p3_withERCC_rsem";


elif(args.reference == "z10"):
	#zebrafish
	KALLISTO_INDEX_FILE = ""

	GENOME_REFERENCE_FILE = "/data/yosef/index_files/z10/index/GCF_000002035.5_GRCz10_genomic.fna";
	BOWTIE2_INDEX = "/data/yosef/index_files/z10/index/GCF_000002035.5_GRCz10_genomic";
	TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/z10/index/tophat_transcriptome_data/GRCz10.refseq_annot";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/z10/index/GCF_000002035.5_GRCz10_genomic.gff";

	#RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/z10/index/rsem_index/GCF_000002035.5_GRCz10_genomic_noProblematicGenes.gtf";
	#RSEM_INDEX = "/data/yosef/index_files/z10/index/rsem_index/GCF_000002035.5_GRCz10_genomic_noProblematicGenes.rsem";
	#the  previous files didn't have the control sequences like the eGFP, tdTomatod etc.
	RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/z10/index/rsem_index/combinedGTF_4brain.gtf";
	RSEM_INDEX = "/data/yosef/index_files/z10/index/rsem_index/GRCz10_4brain.rsem";
	REF_FLAT_INDEX = "/data/yosef/index_files/z10/index/refFlat/refFlat_allonEdited.txt";
	RIBOSOMAL_INTERVALS_INDEX = "null" #I didn't create it - no time and not that important...
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "null" #I didn't create it - no time and not that important...
	RIBOSOMAL_INTERVALS_INDEX_FOR_KALLISTO = "null" #I didn't create it - no time and not that important...

	RSEM_DICTIONARY = "/data/yosef/index_files/z10/index/rsem_index/rsemDictionary/z10_rsemGeneMapping.txt";


#hg19 was deprecated - then restored at Michael's request
elif(args.reference == "hg19"):
	KALLISTO_INDEX_FILE = "";
        GENOME_REFERENCE_FILE = "/data/yosef/index_files/hg19_complete/index/GRCh37_p13.fna";
        BOWTIE2_INDEX = "/data/yosef/index_files/hg19_complete/index/GRCh37_p13";
        TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/hg19_complete/index/tophat_transcriptome_data/GRCh37_p13";
        RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg19_complete/index/rsem_index/GRCh37_p13.gtf";
        TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg19_complete/index/GRCh37_p13.gtf";
        RSEM_INDEX ="/data/yosef/index_files/hg19_complete/index/rsem_index/GRCh37_p13_rsem";
        REF_FLAT_INDEX ="/data/yosef/index_files/hg19_complete/index/refFlat/refFlat.txt";
        RIBOSOMAL_INTERVALS_INDEX ="/data/yosef/index_files/hg19_complete/index/gencode/GRCh37_p13.rRNA.interval";
        RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM ="/data/yosef/index_files/hg19_complete/index/gencode/GRCh37_p13.rsem.rRNA.interval";
        RSEM_DICTIONARY = "/data/yosef/index_files/hg19_complete/index/rsem_index/rsemDictionary/rsem_transcript.dict";

	#settings for hg19 with and other additions required by the HIV project (tba)
	#Compiled by Michael, Feb 2015
	#BOWTIE2_INDEX = "/data/yosef/index_files/hg19_HIV/index/GRCh37.p13";
	#TOPHAT2_TRANSCRIPTOME_INDEX = "";
	#TRANSCRIPT_ANNOTATION = "";
	#RSEM_INDEX = "/data/yosef/index_files/hg19_HIV/index/rsem_index/GRCh37.p13_rsem";
	#REF_FLAT_INDEX = "/data/yosef/index_files/hg19_HIV/index/refFlat/refFlat.txt";
	#RIBOSOMAL_INTERVALS_INDEX = "null";
	#RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "null";

elif(args.reference == "hg38"):
	#settings for hg38
	#Compiled by Michael, Feb 2015
	#BOWTIE2_INDEX = "/data/yosef/index_files/hg38/index/GRCh38";
	#TOPHAT2_TRANSCRIPTOME_INDEX = "";
	#TRANSCRIPT_ANNOTATION = "";
	#RSEM_INDEX = "/data/yosef/index_files/hg38/index/rsem_index/GRCh38_rsem";
	#REF_FLAT_INDEX = "/data/yosef/index_files/hg38/index/refFlat/refFlat.txt";
	#RIBOSOMAL_INTERVALS_INDEX = "null";
	#RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "null";

	KALLISTO_INDEX_FILE = ""

	GENOME_REFERENCE_FILE = "/data/yosef/index_files/hg38/index/hs_ref_GRCh38.ref_name.fna";
	BOWTIE2_INDEX = "/data/yosef/index_files/hg38/index/GRCh38";
	TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/hg38/index/tophat_transcriptome_data/GRCh38";
	RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg38/rsem_index/GRCh38.gtf";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg38/index/GRCh38.gtf";
	RSEM_INDEX ="/data/yosef/index_files/hg38/index/rsem_index/GRCh38_rsem";
	REF_FLAT_INDEX ="/data/yosef/index_files/hg38/index/refFlat/refFlat.txt";
	RIBOSOMAL_INTERVALS_INDEX ="/data/yosef/index_files/hg38/index/gencode/rRNA.interval";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM ="/data/yosef/index_files/hg38/index/gencode/rRNA.rsem.interval";
	RSEM_DICTIONARY = "/data/yosef/index_files/hg38/index/rsem_dict.txt";

elif(args.reference == "hg38_HIV"):
	#settings for hg38 with HIV and GFP Sequences
	KALLISTO_INDEX_FILE = "";
	GENOME_REFERENCE_FILE = "/data/yosef/index_files/hg38_HIV/index/GRCh38_p2_HIV_pCmGFP.fna";
	BOWTIE2_INDEX = "/data/yosef/index_files/hg38_HIV/index/GRCh38.p2";
	TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/hg38_HIV/index/tophat_transcriptome_data/GRCh38.p2";
	RSEM_TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg38_HIV/index/rsem_index/GRCh38.p2.gtf";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg38_HIV/index/GRCh38.p2.gtf";
	RSEM_INDEX ="/data/yosef/index_files/hg38_HIV/index/rsem_index/GRCh38.p2_rsem";
	REF_FLAT_INDEX ="/data/yosef/index_files/hg38_HIV/index/refFlat/refFlat.txt";
 	RIBOSOMAL_INTERVALS_INDEX ="/data/yosef/index_files/hg38_HIV/index/gencode/rRNA.interval";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM ="/data/yosef/index_files/hg38_HIV/index/gencode/rRNA.rsem.interval";	
	RSEM_DICTIONARY = "/data/yosef/index_files/hg38_HIV/index/rsem_dict.txt";

else:
	raise Exception("should not happen - unsupported reference genome");



#sample file names should not have the '.gz' extension --> remove it if it is  present
if(args.sampleFile1[-3:] == '.gz'):
	print("Truncating the .gz extension from sample file 1's name");
	args.sampleFile1 = args.sampleFile1[0:-3];
if(args.paired_end and args.sampleFile2[-3:] == '.gz'):
	print("Truncating the .gz extension from sample file 2's name");
	args.sampleFile2 = args.sampleFile2[0:-3];


if(args.sampleFile1[-6:] != ".fastq") or (args.paired_end and args.sampleFile2[-6:] != ".fastq"):
	print "Input file names are" + args.sampleFile1 + (args.sampleFile2 if args.paired_end else "");
	raise Exception("Error: The input file names indicate these are not fastq files, but only this format is supported at present")


args.output_folder = os.path.expanduser(args.output_folder);
args.sampleFile1 = os.path.expanduser(args.sampleFile1);
if(args.paired_end):
	args.sampleFile2 = os.path.expanduser(args.sampleFile2);


#debug code:
#subprocess.call("rm aaa*", shell=True);
#subprocess.call("rm bbb*", shell=True);

print("**********************************************************");
print("**********************************************************");
print("gunzipping files into temporary files");
gunzipCommand = Template("gunzip -c $SAMPLE_FILE1.gz > $SAMPLE_FILE1").substitute(SAMPLE_FILE1=args.sampleFile1);
print(gunzipCommand)
returnCode = subprocess.call(gunzipCommand, shell=True);
if(returnCode != 0):
	raise Exception("gunzip of sampleFile1 failed");

if(args.paired_end):
	#debug code:
	#subprocess.call("cp TFH-H9_S93_L001_R1_001.fastq.gz aaa.fastq.gz", shell=True);
	#subprocess.call("cp TFH-H9_S93_L001_R2_001.fastq.gz bbb.fastq.gz", shell=True);
	#gunzip without deleting the original file
	#subprocess.call("gunzip -c aaa.fastq.gz > aaa.fastq", shell=True);
	#subprocess.call("gunzip -c bbb.fastq.gz > bbb.fastq", shell=True);

	gunzipCommand = Template("gunzip -c $SAMPLE_FILE2.gz > $SAMPLE_FILE2").substitute(SAMPLE_FILE2=args.sampleFile2);
	print(gunzipCommand)
	returnCode = subprocess.call(gunzipCommand, shell=True);
	if(returnCode != 0):
		raise Exception("gunzip of sampleFile2 failed");

#else:
	#debug code:
	#subprocess.call("cp GBCP03_N712_S508_GTAGAGGA-CTAAGCCT_L008_R1_001.fastq.gz aaa.fastq.gz", shell=True);
	#subprocess.call("gunzip -c aaa.fastq.gz > aaa.fastq", shell=True);

if not os.path.exists(args.output_folder):
	os.makedirs(args.output_folder);

#gunzipped files to clear at the end of the run
filesToClear = [args.sampleFile1];
if(args.paired_end):
	filesToClear = filesToClear + [args.sampleFile2]


DO_TRIMMOMATIC = True;
#enter only if the user didn't select --skip_trimmomatic OR he did, but he did not mark the --do_not_rely_on_previous_trimmomatic (which means he does want to rely on previous output)
if(DO_TRIMMOMATIC and not(args.skip_trimmomatic and args.do_not_rely_on_previous_trimmomatic)):
	print("**********************************************************");
	print("**********************************************************");

	trimInput1 = args.sampleFile1;
	if(args.paired_end):
		trimIsPaired = "PE";
		trimInput2 = args.sampleFile2;
		trimOutputParts = [(args.output_folder + "/trimmomatic_output/" + part) for part in ["1.Ptrim.fastq", "1.Utrim.fastq", "2.Ptrim.fastq", "2.Utrim.fastq"]];
		trimOutput = ' '.join(trimOutputParts);
		#trimOutput = Template("$OUTPUT_FOLDER/trimmomatic_output/1.Ptrim.fastq $OUTPUT_FOLDER/trimmomatic_output/1.Utrim.fastq $OUTPUT_FOLDER/trimmomatic_output/2.Ptrim.fastq $OUTPUT_FOLDER/trimmomatic_output/2.Utrim.fastq").substitute(OUTPUT_FOLDER=args.output_folder);

	else:
		trimIsPaired = "SE";
		trimInput2 = "";
		trimOutput = Template("$OUTPUT_FOLDER/trimmomatic_output/1.Utrim.fastq").substitute(OUTPUT_FOLDER=args.output_folder);
		trimOutputParts = [trimOutput];

	if(args.skip_trimmomatic):
		#if we reached this point, the user wants to rely on a previous trimmomatic run, rather than recompute trimmomatic

		print("relying on a previous Trimmomatic run");
		print("guzipping previous Trimmomatic results...");
		for trimFile in trimOutputParts:
			print trimFile;
			print trimOutputParts;
			print args.output_folder;

			gunzipCommand = Template("gunzip -c $INFILE.gz > $INFILE").substitute(INFILE=trimFile);
			print(gunzipCommand)
			returnCode = subprocess.call(gunzipCommand, shell=True);
			if(returnCode != 0):
				raise Exception("gunzip of trimmomatic output file failed");

		filesToClear = filesToClear + trimOutputParts;
	else:
		#if we reached this point, the user wants to run trimmomatic...

		print("doing Trimmomatic");
		sys.stdout.flush();

		if not os.path.exists(args.output_folder + "/trimmomatic_output"):
			os.makedirs(args.output_folder + "/trimmomatic_output");

		#infer the length of reads in the sample by looking at the length of the first read and assuming that all reads are of the same length
		read_length = len(rnaSeqPipelineUtils.GetFirstReadInFastqFile(args.sampleFile1));
		print "Inferred read length in sample is: " + str(read_length);
		if(args.paired_end):
			read_length_otherSide = len(rnaSeqPipelineUtils.GetFirstReadInFastqFile(args.sampleFile2));
			print "Inferred read length of the other end in the paired-end sample is: " + str(read_length_otherSide);
			if read_length != read_length_otherSide:
				print("---> Warning: Trimmomatic step warning: it seems that the two ends of the paired-end experiment have different read lengths.");


		trimMinLen = min(50, int(0.8 * read_length));#16; #36
		trimCrop = ""; #"CROP:50"; --> do not do trimCrop
		trimMinPhred = 15;
		trimWindow = args.trimmomatic_window if args.trimmomatic_window else "4:15"

		trimCmd = Template("java -jar /opt/pkg/Trimmomatic-0.32/trimmomatic-0.32.jar $TRIM_IS_PAIRED -threads $NUM_THREADS -phred33 -trimlog $OUTPUT_FOLDER/trimmomatic_output/trimmomatic_log.txt $TRIM_INPUT1 $TRIM_INPUT2 $TRIM_OUTPUT LEADING:$TRIM_MIN_PHRED TRAILING:$TRIM_MIN_PHRED SLIDINGWINDOW:$TRIM_WINDOW MINLEN:$TRIM_MINLEN $TRIM_CROP").substitute(TRIM_IS_PAIRED=trimIsPaired, OUTPUT_FOLDER=args.output_folder, TRIM_INPUT1=trimInput1, TRIM_INPUT2=trimInput2, TRIM_OUTPUT=trimOutput, TRIM_MINLEN=trimMinLen, TRIM_MIN_PHRED=trimMinPhred, TRIM_CROP=trimCrop, NUM_THREADS=args.num_threads, TRIM_WINDOW=trimWindow);
		print(trimCmd)
		sys.stdout.flush();
		returnCode = subprocess.call(trimCmd, shell=True);
		if(returnCode != 0):
			raise Exception("trimmomatic failed");

		print("**********************************************************");
		print("**********************************************************");
		print("gzipping Trimmomatic results to back them up...");
		for trimFile in trimOutputParts:
			print trimFile;
			print trimOutputParts;
			print args.output_folder;

			gzipCommand = Template("gzip -c $INFILE > $INFILE.gz").substitute(INFILE=trimFile);
			print(gzipCommand)
			returnCode = subprocess.call(gzipCommand, shell=True);
			if(returnCode != 0):
				raise Exception("gzip of trimmomatic output file failed");

		filesToClear = filesToClear + trimOutputParts;


	print("**********************************************************");
	print("**********************************************************");
	print("Trimmomatic succeeded, overriding original input files");

	#replaced the original inputs with the trimmed files
	if(args.paired_end):
		args.sampleFile1 = Template("$OUTPUT_FOLDER/trimmomatic_output/1.Ptrim.fastq").substitute(OUTPUT_FOLDER=args.output_folder);
		args.sampleFile2 = Template("$OUTPUT_FOLDER/trimmomatic_output/2.Ptrim.fastq").substitute(OUTPUT_FOLDER=args.output_folder);
	else:
		args.sampleFile1 = trimOutput;



	print("**********************************************************");
	print("**********************************************************");
	print("checking that Trimmomatic output is not empty...");

	if(os.stat(args.sampleFile1).st_size == 0):
		print "Danger, Will Robinson, sample file 1 was left empty after running trimmomatic!"
	if(args.paired_end and (os.stat(args.sampleFile2).st_size == 0)):
		print "Danger, Will Robinson, sample file 2 was left empty after running trimmomatic!"

	print("**********************************************************");
	print("**********************************************************");
	print("trimmomatic done!");


if(not(KALLISTO_INDEX_FILE)):
	print("Kallisto index has not been defined... skipping Kallisto");
	args.skip_kallisto = True

DO_KALLISTO = True;
if(DO_KALLISTO and not(args.skip_kallisto)):
	print("**********************************************************");
	print("**********************************************************");
	print("Running Kallisto")

	#I decided to do Kallisto before trimmomatic because I was afraid that Kallisto's fragment length estimation will not work correctly
	#if we start trimming reads (and anyhow, once you do not align but rather pseudo-align errors in BPs are less important
	#Allon 15 Jun 2016: cufflinks also does its own estimation in paired-end, so I reconsidered this and moved Kallisto to after trimming

	if not os.path.exists(args.output_folder + "/kallisto_output"):
		os.makedirs(args.output_folder + "/kallisto_output");

	#I do not use multithreaded pseudoalignment because it is not supported with creation of pseudobam
	kallistoDesiredOutputType = "" if args.kallisto_write_hdf5 else "--plaintext" #hdf5 is the default
	kallistoCmd = Template("kallisto quant $KALLISTO_OUTPUT_TYPE --pseudobam -i $KALLISTO_INDEX -o $OUTPUT_FOLDER -b $KALLISTO_BOOTSTRAP_SAMPLES -t 1").\
		substitute(KALLISTO_INDEX=KALLISTO_INDEX_FILE, OUTPUT_FOLDER=args.output_folder + "/kallisto_output",
				   KALLISTO_BOOTSTRAP_SAMPLES=args.kallisto_bootstrap_samples, KALLISTO_OUTPUT_TYPE=kallistoDesiredOutputType);

	if(args.paired_end):
		kallistoCmd += Template(" $SAMPLE_FILE1 $SAMPLE_FILE2").substitute(SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2)
	else:
		kallistoCmd += Template(" --single -l $MEAN_FRAGMENT_LENGTH -s $STD_FRAGMENT_LENGTH $SAMPLE_FILE1").substitute(MEAN_FRAGMENT_LENGTH=args.mean_fragment_length, STD_FRAGMENT_LENGTH=args.std_fragment_length, SAMPLE_FILE1=args.sampleFile1)

	#Kallisto outputs the pseudobam to stdout
	kallistoCmd += Template(" | /opt/pkg/samtools-1.3.1/bin/samtools view -Sb - > $OUTPUT_FILE").substitute(OUTPUT_FILE=os.path.join(args.output_folder, "kallisto_output", "kallisto_out.bam"))

	print(kallistoCmd)
	sys.stdout.flush();
	returnCode = subprocess.call(kallistoCmd, shell=True);
	if(returnCode != 0):
		raise Exception("Kallisto failed");

	print("**********************************************************");
	print("**********************************************************");
	print("Kallisto quantification done")



	print("**********************************************************");
	print("**********************************************************");
	picardDir = Template("$OUTPUT_FOLDER/kallisto_output/picard_output").substitute(OUTPUT_FOLDER=args.output_folder);
	if not(os.path.exists(picardDir)):
		os.makedirs(picardDir);

	#turns out that Kallisto outputs unmapped reads with non-zero MAPQ (mappying qualities) and Picard doesn't like that and throws an error on unmapped reads with non-zero MAPQ
	#At first I tried to run first Picard's CleanSam utility on kallisto's output that should overwrite the MAPQ of any unmappped reads to 0 but then it also threw errors.
	#so instead, I just switched to using samtools sort utility which is actually preferable because it supports multiple threads (not that it will matter much when running through the cluster which sets the number of threads to 1 per process anyhow)

	# print("Clean SAM (required to set MAPQ to 0 in unaligned reads. Turns out Kallisto generated unmapped reads with non-zero MAPQ")
	# cleanSamCommand = Template("picard CleanSam INPUT=$OUTPUT_FOLDER/kallisto_output/kallisto_out.bam OUTPUT=$OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.cleaned.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	# print(cleanSamCommand)
	# sys.stdout.flush();
	# returnCode = subprocess.call(cleanSamCommand, shell=True);
	# if(returnCode != 0):
	# 	raise Exception("CleanSam failed");
    #
    #
	# print("Sort SAM");
	# sortSamCommand = Template("picard SortSam TMP_DIR=$OUTPUT_FOLDER/temp I=$OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.cleaned.bam O=$OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.sorted.bam SO=coordinate").substitute(OUTPUT_FOLDER=args.output_folder);
	# print(sortSamCommand)
	# sys.stdout.flush();
	# returnCode = subprocess.call(sortSamCommand, shell=True);
	# if(returnCode != 0):
	# 	raise Exception("sortSam failed");


	print("**********************************************************");
	print("**********************************************************");
	print("Sort SAM");
	sortSamCommand = Template("/opt/pkg/samtools-1.3.1/bin/samtools sort -@ $NUM_THREADS -T $OUTPUT_FOLDER/temp -o $OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.sorted.bam $OUTPUT_FOLDER/kallisto_output/kallisto_out.bam").substitute(OUTPUT_FOLDER=args.output_folder, NUM_THREADS=args.num_threads);
	print(sortSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(sortSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("sortSam failed");

	print("**********************************************************");
	print("**********************************************************");
	print("Index SAM");
	indexSamCommand = Template("/opt/pkg/samtools-1.3.1/bin/samtools index $OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(indexSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(indexSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("indexSam failed");

	
	print("**********************************************************");
	print("**********************************************************");
	#to conform with the tophat pipeline, I split the reads into aligned and unaligned, and the QC will run only on the aligned ones.

	print("splitting kallisto's output to aligned and unaligned reads...");
	splitCmd1 = Template("/opt/pkg/samtools-1.3.1/bin/samtools view -b -F 4 $OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.sorted.bam > $OUTPUT_FOLDER/kallisto_output/accepted_hits.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd1)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd1, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");


	splitCmd2 = Template("/opt/pkg/samtools-1.3.1/bin/samtools view -b -f 4 $OUTPUT_FOLDER/kallisto_output/picard_output/kallisto_out.sorted.bam > $OUTPUT_FOLDER/kallisto_output/unmapped.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd2)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd2, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");


	print("**********************************************************");
	print("**********************************************************");
	print("WARNING: In this step multiple alignments per read should be collapsed into one, but this is not implemented yet for kallisto (only for rsem)");
	print("As a result, the number of reads and the number of aligned reads in the QC will be wrong");
	shutil.copyfile(os.path.join(args.output_folder, "kallisto_output/accepted_hits.bam"),
					os.path.join(args.output_folder, "kallisto_output/accepted_hits_noMultiple.bam"));





	print("**********************************************************");
	print("**********************************************************");
	print("Kallisto pipeline done!");

#here don't forget to delete temporary files + QC

#I use the terms "cufflinks pipeline" and "tophat pipeline" interchangeably
RUN_CUFFLINKS_PIPELINE = True;
if(RUN_CUFFLINKS_PIPELINE and not(args.skip_tophat)):

	print("**********************************************************");
	print("**********************************************************");
	print("starting cufflinks pipeline");
	print("TODO (Mar. 15): WARNING: the tophat pipeline will output wrong number of reads and ratio of aligned reads in the summary,txt - this is the same problem that was in the rsem pipeline, the accepted_hits.bam will include multiple reads and they have to be filtered out before counting");
	print("TODO (Mar. 15): WARNING: the script count_dup_per_gene.pl relies on some rsem-specific dictionary and needs to be modified...");
	sys.stdout.flush();

	if not os.path.exists(args.output_folder + "/tophat_output"):
		os.makedirs(args.output_folder + "/tophat_output");

	if(args.paired_end):
		tophatCommand = Template("tophat2 --num-threads $NUM_THREADS -o $OUTPUT_FOLDER/tophat_output --transcriptome-index $TOPHAT2_TRANSCRIPTOME_INDEX $BOWTIE2_INDEX $SAMPLE_FILE1 $SAMPLE_FILE2").substitute(OUTPUT_FOLDER=args.output_folder, TOPHAT2_TRANSCRIPTOME_INDEX=TOPHAT2_TRANSCRIPTOME_INDEX, BOWTIE2_INDEX=BOWTIE2_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2, NUM_THREADS=args.num_threads);
	else:
		tophatCommand = Template("tophat2 --num-threads $NUM_THREADS -o $OUTPUT_FOLDER/tophat_output --transcriptome-index $TOPHAT2_TRANSCRIPTOME_INDEX $BOWTIE2_INDEX $SAMPLE_FILE1").substitute(OUTPUT_FOLDER=args.output_folder, TOPHAT2_TRANSCRIPTOME_INDEX=TOPHAT2_TRANSCRIPTOME_INDEX, BOWTIE2_INDEX=BOWTIE2_INDEX, SAMPLE_FILE1=args.sampleFile1, NUM_THREADS=args.num_threads);

	print("**********************************************************");
	print("**********************************************************");
	print("Running TopHat");
	print(tophatCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(tophatCommand, shell=True);
	if(returnCode != 0):
		raise Exception("tophat failed");

	print("**********************************************************");
	print("**********************************************************");
	print("WARNING: In this step multiple alignments per read should be collapsed into one, but this is not implemented yet for tophat (only for rsem)");
	print("As a result, the number of reads and the number of aligned reads in the QC will be wrong");
	shutil.copyfile(os.path.join(args.output_folder, "tophat_output/accepted_hits.bam"),
					os.path.join(args.output_folder, "tophat_output/accepted_hits_noMultiple.bam"));


	print("**********************************************************");
	print("**********************************************************");
	picardDir = Template("$OUTPUT_FOLDER/tophat_output/picard_output").substitute(OUTPUT_FOLDER=args.output_folder);
	if not(os.path.exists(picardDir)):
		os.makedirs(picardDir);

	print("Sort SAM");
	sortSamCommand = Template("picard SortSam TMP_DIR=$OUTPUT_FOLDER/temp I=$OUTPUT_FOLDER/tophat_output/accepted_hits_noMultiple.bam O=$OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam SO=coordinate").substitute(OUTPUT_FOLDER=args.output_folder);
	print(sortSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(sortSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("sortSam failed");



	print("**********************************************************");
	print("**********************************************************");
	print("Index SAM");
	indexSamCommand = Template("/opt/pkg/samtools-1.3.1/bin/samtools index $OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(indexSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(indexSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("indexSam failed");

	print("**********************************************************");
	print("**********************************************************");
	print("Running Cufflinks");

	#in cufflinks, -m/-frag-len-mean is the expected (mean) fragment length - learned automatically for paired-end so I use it only for single end
	if(args.paired_end):
		cufflinksSingleEndArgs = ""
	else:
		#STD_FRAG_LEN should be args.std_fragment_length BUT I found too large values (as in the brain project) lead to core dumps, so I use the cufflinks default which is 80
		cufflinksSingleEndArgs = Template("--frag-len-mean $MEAN_FRAG_LEN --frag-len-std-dev $STD_FRAG_LEN").substitute(MEAN_FRAG_LEN=args.mean_fragment_length, STD_FRAG_LEN=80)


	cufflinksComand = Template("/opt/genomics/bin/cufflinks $SINGLE_END_ARGS --num-threads $NUM_THREADS -G $TRANSCRIPT_ANNOTATION -o $OUTPUT_FOLDER/tophat_output/cuff_output/ $OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder, TRANSCRIPT_ANNOTATION=TRANSCRIPT_ANNOTATION, NUM_THREADS=args.num_threads, SINGLE_END_ARGS=cufflinksSingleEndArgs);
	cufflinksComand = " ".join(cufflinksComand.split()) #replace multiple whitespaces with one
	print(cufflinksComand)
	sys.stdout.flush();
	returnCode = subprocess.call(cufflinksComand, shell=True);
	if(returnCode != 0):
		raise Exception("cufflinks failed");


	print("**********************************************************");
	print("**********************************************************");
	print("Counting features on the Tophat2 alignment");
	#for featureCounts output, see: https://groups.google.com/forum/#!topic/subread/JPPw9lVfgpw

	#IMPORTANT: I give the gtf file that was converted from the GFF3 with cufflink's gffread utilitiy (at least for the mm10_4brain),
	# even though tophat itself used the gff3 annotation. This is because featureCounts wants a gtf file...

	#I don't count reads mapping to multiple meta-features which the manual does not recommend for RNA-Seq for the following reason:
	#A multi-overlap read or fragment is one that overlaps more than one feature, or more than
	#one meta-feature when summarizing at the meta-feature level. featureCounts provides users
	#with the option to either exclude multi-overlap reads or to count them for each feature that
	#is overlapped. The decision whether or not to counting these reads is often determined by the
	#experiment type. We recommend that reads or fragments overlapping more than one gene
	#are not counted for RNA-seq experiments, because any single fragment must originate from
	#only one of the target genes but the identity of the true target gene cannot be confidently
	#determined. On the other hand, we recommend that multi-overlap reads or fragments are
	#counted for most ChIP-seq experiments because epigenetic modifications inferred from these
	#reads may regulate the biological functions of all their overlapping genes.


	#-C: exclude chimeric fragments from the count (those fragments that have their two ends aligned to different chromosomes)
	pairedEndFeatureCountArgs = "-p -C" if(args.paired_end) else ''
	countCmd = Template("/opt/genomics/bin/featureCounts " + pairedEndFeatureCountArgs + " -R -T $NUM_THREADS -t exon -g gene_id -a $ANNOTATION" +
						" -o $OUTPUT_FOLDER/tophat_output/feature_counts.txt $OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam").substitute(NUM_THREADS=args.num_threads, ANNOTATION=RSEM_TRANSCRIPT_ANNOTATION, OUTPUT_FOLDER=args.output_folder)


	print(countCmd)
	sys.stdout.flush();
	returnCode = subprocess.call(countCmd, shell=True);
	if(returnCode != 0):
		raise Exception("tophat2 count features failed");



	print("**********************************************************");
	print("**********************************************************");
	print("cufflinks pipeline done!");


RUN_RSEM_PIPELINE = True;
if(RUN_RSEM_PIPELINE and not(args.skip_rsem)):

	print("**********************************************************");
	print("**********************************************************");
	print("starting rsem pipeline");

	if not os.path.exists(args.output_folder + "/rsem_output/picard_output"):
		os.makedirs(args.output_folder + "/rsem_output/picard_output");

	print("**********************************************************");
	print("**********************************************************");

	RUN_DIRECT_RSEM = True #Run RSEM directly or run manually bowtie2 and then rsem
	if(RUN_DIRECT_RSEM):
		# previous code: since then I have chosen to run bowtie2 manually
		print("Running rsem");
		#Note that in rsem I use the --output-genome-bam flag to generate a genome bam (in addition to the transcriptome bam that is always created) - rsem will also sort this bam. This is necessary for the QC later.
		#the --append-names is useful but I don't use it because it's incompatible with the dictionaries I've produced in the past
		rsemCommand = Template("/opt/pkg/rsem-1.2.31/bin/rsem-calculate-expression --num-threads $NUM_THREADS --bowtie2 $SINGLE_CELL_OR_BULK_PRIOR --sampling-for-bam --output-genome-bam --sort-bam-by-coordinate --sort-bam-memory-per-thread $RSEM_SAMTOOLS_SORT_MEM --estimate-rspd ").substitute(SINGLE_CELL_OR_BULK_PRIOR="--single-cell-prior" if args.single_cell_prior else "", NUM_THREADS=args.num_threads, RSEM_SAMTOOLS_SORT_MEM=args.rsem_samtools_sort_mem);
		if(args.paired_end):
			rsemCommand += Template(" --paired-end $SAMPLE_FILE1 $SAMPLE_FILE2 ").substitute(SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);
		else:
			#the fragment length distribution is used only in single-end
			rsemCommand += Template(" --fragment-length-mean $MEAN_FRAGMENT_LENGTH --fragment-length-sd $STD_FRAGMENT_LENGTH $SAMPLE_FILE1").substitute(SAMPLE_FILE1=args.sampleFile1, MEAN_FRAGMENT_LENGTH=args.mean_fragment_length, STD_FRAGMENT_LENGTH=args.std_fragment_length)

 		rsemCommand += Template(" $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(RSEM_INDEX=RSEM_INDEX, OUTPUT_FOLDER=args.output_folder)
		rsemCommand = ' '.join(rsemCommand.split()) #replace double whitespaces with one whitespace
		print(rsemCommand)
		sys.stdout.flush();
		returnCode = subprocess.call(rsemCommand, shell=True);
		if(returnCode != 0):
			raise Exception("rsem failed");


	else:
		raise Exception("Reached deprecated code. At first I ran bowtie2 directly because rsem didn't allow configuring bowtie's fine paramteres. But since then they've added support, and since rsem configures the bowtie2 params to avoid indels etc. it's best to let rsme run as usual (in case, for example, the documentation on how they set the bowtie2 params is incomplete")

		#I am running bowtie manually and then giving the output to rsem because sometimes the bowtie parameters need to be tweaked and rsem does not support changes of many of the parameters at present
		print("Running bowtie for rsem")
		#this are the params that rsem delivers to bowtie2 by default
		bowtieParams = Template("-q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 " +
								"-x $RSEM_INDEX -p $NUM_THREADS -k 200").substitute(RSEM_INDEX=RSEM_INDEX, NUM_THREADS=args.num_threads)

		if(args.paired_end):
			#additional params that are specific to paired-end
			bowtieParams += Template(" --no-mixed --no-discordant -I 1 -X $RSEM_BOWTIE_MAXINS").substitute(RSEM_BOWTIE_MAXINS=args.rsem_bowtie_maxins)
			bowtieInput = Template("-1 $SAMPLE_FILE1 -2 $SAMPLE_FILE2").substitute(SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2)
		else:
				bowtieInput = Template("-U $SAMPLE_FILE1").substitute(SAMPLE_FILE1=args.sampleFile1)

		bowtieForRsemCmd = Template("/opt/genomics/bin/bowtie2 $BOWTIE_PARAMS $BOWTIE_INPUT | /opt/pkg/samtools-1.3.1/bin/samtools view -S -b -o $OUTPUT_FOLDER/rsem_output/aligned_by_bowtie2.bam -").substitute(BOWTIE_INPUT=bowtieInput, BOWTIE_PARAMS=bowtieParams, OUTPUT_FOLDER=args.output_folder)
		print(bowtieForRsemCmd)
		sys.stdout.flush();
		returnCode = subprocess.call(bowtieForRsemCmd, shell=True);
		if(returnCode != 0):
			raise Exception("bowtie2 for rsem failed");

		#--sampling-for-bam: rsem outputs a bam with multiple possible alignments per read plus their posterior probabilities, which means that the read counts are disrupted
		#this flag tells it to sample one read per the distribution, assign it a MAPQ=100 value (also tag ZW:f:1) and the rest are 0 value (and tag ZW:f:0)

		print("Running rsem");
		#params that are relevant only in paired_end run
		rsemParamsOnlyForPairedEnd = Template("--fragment-length-max $RSEM_BOWTIE_MAXINS").substitute(RSEM_BOWTIE_MAXINS=args.rsem_bowtie_maxins) if args.paired_end else ""
		#Note that in rsem I use the --output-genome-bam flag to generate a genome bam (in addition to the transcriptome bam that is always created) - rsem will also sort this bam. This is necessary for the QC later.
		rsemCommand = Template("/opt/pkg/rsem-1.2.31/bin/rsem-calculate-expression --sort-bam-by-coordinate --num-threads $NUM_THREADS --estimate-rspd " +
								"--samtools-sort-mem $RSEM_SAMTOOLS_SORT_MEM " +
								"$RSEM_PARAMS_ONLY_FOR_PAIRED_END --output-genome-bam --sampling-for-bam --bam $IS_PAIRED_END $OUTPUT_FOLDER/rsem_output/aligned_by_bowtie2.bam " +
								"$RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder,
																										   RSEM_INDEX=RSEM_INDEX,
																										   NUM_THREADS=args.num_threads,
																										   RSEM_PARAMS_ONLY_FOR_PAIRED_END=rsemParamsOnlyForPairedEnd,
																										   IS_PAIRED_END="--paired-end" if args.paired_end else "",
																										   RSEM_SAMTOOLS_SORT_MEM=args.rsem_samtools_sort_mem);

		print(rsemCommand)
		sys.stdout.flush();
		returnCode = subprocess.call(rsemCommand, shell=True);
		if(returnCode != 0):
			raise Exception("rsem failed");


	print("**********************************************************");
	print("**********************************************************");
	print("Counting features on the Bowtie2-RSEM alignment");
	#for featureCounts output, see: https://groups.google.com/forum/#!topic/subread/JPPw9lVfgpw

	#-C: exclude chimeric fragments from the count (those fragments that have their two ends aligned to different chromosomes)
	pairedEndFeatureCountArgs = "-p -C" if(args.paired_end) else ''
	countCmd = Template("/opt/genomics/bin/featureCounts " + pairedEndFeatureCountArgs + " -R -T $NUM_THREADS -t exon -g gene_id -a $ANNOTATION" +
						" -o $OUTPUT_FOLDER/rsem_output/feature_counts.txt $OUTPUT_FOLDER/rsem_output/rsem_output.genome.sorted.bam").substitute(NUM_THREADS=args.num_threads, ANNOTATION=RSEM_TRANSCRIPT_ANNOTATION, OUTPUT_FOLDER=args.output_folder)

	if(False):
		#IMPORTANT NOTE: running featureCounts on RSEM's output never gives Unassigned_MultiMapping. How come no reads are assigned across features?
		#TURNS OUT that featureCounts identifies multiply-mapped reads by the presence of the NH sam tag. This flag is present in the tophat2 output, so
		#featureCounts correctly identifies multiply-mapped reads there, but not in rsem genome output bam that does not have this tag... (no read has this tag, so all reads are considered as non-multiply mapped)

		#THIS IS WHY I DISABLED featureCounts for RSEM - it didn't seem important enough to do that on top of the tophat feature counting. If it becomes important enough, this problem should be solved first

		# ALSO: I checked htseq-count and it seemed to have the same problem - did not find multiply aligned reads in rsem's output, but many of them in tophat's output
		# The ht-seq command I used:
		# htseq-count -f bam rsem_output.genome.sorted.bam /data/yosef/index_files/mm10_4brain/index/rsem_index/combinedGTF_4brain.gtf

		print(countCmd)
		sys.stdout.flush();
		returnCode = subprocess.call(countCmd, shell=True);
		if(returnCode != 0):
			raise Exception("bowtie2-rsem count features failed");
	else:
		print("Skipping featureCounts with RSEM (rsem output bam doesn't have the NH tag that allows correct identification of multiply aligned reads")


	print("**********************************************************");
	print("**********************************************************");
	#to conform with the tophat pipeline, I split the reads into aligned and unaligned, and the QC will run only on the aligned ones.

	print("splitting rsem's output to aligned and unaligned reads...");
	splitCmd1 = Template("/opt/pkg/samtools-1.3.1/bin/samtools view -b -F 4 $OUTPUT_FOLDER/rsem_output/rsem_output.genome.sorted.bam > $OUTPUT_FOLDER/rsem_output/accepted_hits.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd1)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd1, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");


	splitCmd2 = Template("/opt/pkg/samtools-1.3.1/bin/samtools view -b -f 4 $OUTPUT_FOLDER/rsem_output/rsem_output.genome.sorted.bam > $OUTPUT_FOLDER/rsem_output/unmapped.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd2)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd2, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");



	print("removing the multiple-posterior alignments from rsem's aligned output reads...");
	#with the --output-genome-bam flag set, rsem outputs a bam in which the MAPQ are either 1 (chosen) or 0 (not chosen)
	filterCmd = Template("/opt/pkg/samtools-1.3.1/bin/samtools view -b -q 100 $OUTPUT_FOLDER/rsem_output/accepted_hits.bam > $OUTPUT_FOLDER/rsem_output/accepted_hits_noMultiple.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(filterCmd)
	sys.stdout.flush();
	returnCode = subprocess.call(filterCmd, shell=True);
	if(returnCode != 0):
		raise Exception("filtering failed");


	print("Sort SAM");
	sortSamCommand = Template("picard SortSam TMP_DIR=$OUTPUT_FOLDER/temp I=$OUTPUT_FOLDER/rsem_output/accepted_hits_noMultiple.bam O=$OUTPUT_FOLDER/rsem_output/picard_output/sorted.bam SO=coordinate").substitute(OUTPUT_FOLDER=args.output_folder);
	print(sortSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(sortSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("sortSam failed");

	print("**********************************************************");
	print("**********************************************************");
	print("Index SAM");
	indexSamCommand = Template("/opt/pkg/samtools-1.3.1/bin/samtools index $OUTPUT_FOLDER/rsem_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(indexSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(indexSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("indexSam failed");



	print("**********************************************************");
	print("**********************************************************");
	print("rsem pipeline done!");


RUN_QC = True;
if(RUN_QC and not(args.skip_qc)):
	print("**********************************************************");
	print("**********************************************************");
	qcOutputDir = Template("$OUTPUT_FOLDER/fastqc_output").substitute(OUTPUT_FOLDER=args.output_folder);
	if not(os.path.exists(qcOutputDir)):
		os.makedirs(qcOutputDir);

	print("Running QC script");
	#note that the output folder dir argument that I give is without the "/fastqc_output" suffix because the qc script expects it that way
	print("Running on " + args.sampleFile1);
	doQC.RunFastQC(args.sampleFile1, args.output_folder);
	doQC.CheckContam(args.sampleFile1, args.output_folder, "1");

	if(args.paired_end):
		print("Running on " + args.sampleFile2);
		doQC.RunFastQC(args.sampleFile2, args.output_folder);
		doQC.CheckContam(args.sampleFile2, args.output_folder, "2");

	#for consistency with doQC.RunFastQC I give the output folder without the "/fastqc_output" suffix
	commonMetrics = doQC.CollectFastQcData(args.output_folder, args.paired_end);
	#commonMetrics - common to both the rsem and tophat pipelines

	#Kallisto QC is currently disabled because Kallisto produces BAM files that throw exceptions in Picard, even when Picard's stringency is set to low
	RUN_QC_KALLISTO = True;
	if(RUN_QC_KALLISTO and not(args.skip_kallisto_qc)):
		sortedBamFile = args.output_folder + "/kallisto_output/accepted_hits.bam";
		kallistoOutputFolder = args.output_folder + "/kallisto_output";
		print "Warning: the kallisto pipleline still doesn't have the two dictionary required for the collect dup perl... this part will be skipped..."
		#weird, the "picard CollectAlignmentSummaryMetrics" didn't work with the genome reference, but works without it contrary to what I had long ago in cufflinks...
		kallisto_qc_metrics = doQC.CollectData(sortedBamFile, kallistoOutputFolder, REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX_FOR_KALLISTO, args.paired_end, RSEM_TRANSCRIPT_ANNOTATION, RSEM_DICTIONARY, "null"); #patch - use the same dictionaries for tophat and rsem for Nir's count dups
		doQC.WriteQCMetrics(kallistoOutputFolder, commonMetrics, kallisto_qc_metrics);

	RUN_QC_TOPHAT = True;
	if(RUN_QC_TOPHAT and not(args.skip_tophat_qc)):
		sortedBamFile = args.output_folder + "/tophat_output/picard_output/sorted.bam";
		tophatOutputFolder = args.output_folder + "/tophat_output";
		print "Warning: the cufflinks pipleline still doesn't have the two dictionary required for the collect dup perl... this part will be skipped..."
		tophat_qc_metrics = doQC.CollectData(sortedBamFile, tophatOutputFolder, REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX, args.paired_end, RSEM_TRANSCRIPT_ANNOTATION, RSEM_DICTIONARY, GENOME_REFERENCE_FILE); #patch - use the same dictionaries for tophat and rsem for Nir's count dups
		doQC.WriteQCMetrics(tophatOutputFolder, commonMetrics, tophat_qc_metrics);

	RUN_QC_RSEM = True;
	if(RUN_QC_RSEM and not(args.skip_rsem_qc)):
		rsemPicardDir = args.output_folder + "/rsem_output/picard_output/";
		rsemOutputFolder = args.output_folder + "/rsem_output";
		if not(os.path.exists(rsemPicardDir)):
			os.makedirs(rsemPicardDir);

		#sortedBamFile = args.output_folder + "/rsem_output/rsem_output.genome.sorted.bam";
		#use only aligned reads for qc
		sortedBamFile = args.output_folder + "/rsem_output/picard_output/sorted.bam";
		rsem_qc_metrics = doQC.CollectData(sortedBamFile, rsemOutputFolder, REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM, args.paired_end, RSEM_TRANSCRIPT_ANNOTATION, RSEM_DICTIONARY, GENOME_REFERENCE_FILE);
		doQC.WriteQCMetrics(rsemOutputFolder, commonMetrics, rsem_qc_metrics);






	#if(args.paired_end):
	#	qcComand = Template("perl qc.pl $QC_OUTPUT_DIR  $REF_FLAT_FILE $RIBOSOMAL_INTERVALS_INDEX $SAMPLE_FILE1 $SAMPLE_FILE2").substitute(QC_OUTPUT_DIR=args.output_folder, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2, REF_FLAT_FILE=REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX=RIBOSOMAL_INTERVALS_INDEX);
	#else:
	#	qcComand = Template("perl qc.pl $QC_OUTPUT_DIR  $REF_FLAT_FILE $RIBOSOMAL_INTERVALS_INDEX $SAMPLE_FILE1").substitute(QC_OUTPUT_DIR=args.output_folder, SAMPLE_FILE1=args.sampleFile1, REF_FLAT_FILE=REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX=RIBOSOMAL_INTERVALS_INDEX);


	#print(qcComand)
	#returnCode = subprocess.call(qcComand, shell=True);
	#if(returnCode != 0):
	#	raise Exception("qc failed");

	print("**********************************************************");
	print("**********************************************************");
	print("QC done!");



print("**********************************************************");
print("**********************************************************");
print("removing temporary gunzipped file");
for remFile in filesToClear:
	rmCommand = Template("rm -f $REM_FILE").substitute(REM_FILE=remFile);
	print(rmCommand)
	returnCode = subprocess.call(rmCommand, shell=True);
	if(returnCode != 0):
		raise Exception("rm of file to clear failed");

if(not(args.do_not_clean_intermediary_files)):
	cleanup_after_preproc.cleanTemporaryFiles(args.output_folder)


print("**********************************************************");
print("**********************************************************");
print("creating success token file");
tokenFileName = os.path.join(args.output_folder, "success_token.txt")
with open(tokenFileName, "wt") as fout:
	fout.write("I finished successfully!\n")

#print("**********************************************************");
#print("**********************************************************");
#print("NOT CALLED - FIX ME - removing temporary gunzipped file");
#rmCommand = Template("rm -f $SAMPLE_FILE1").substitute(SAMPLE_FILE1=args.sampleFile1);
#print(rmCommand)
##returnCode = subprocess.call(rmCommand, shell=True);
#if(returnCode != 0):
#	raise Exception("rm of sampleFile1 failed");
#	
#if(args.paired_end):
#	rmCommand = Template("rm -f $SAMPLE_FILE2").substitute(SAMPLE_FILE2=args.sampleFile2);
#	print(rmCommand)
#	#returnCode = subprocess.call(rmCommand, shell=True);
#	if(returnCode != 0):
#		raise Exception("rm of sampleFile2 failed");


print("**********************************************************");
print("**********************************************************");
print "processing complete"
