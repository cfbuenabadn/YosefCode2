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


parser = argparse.ArgumentParser(description="Process one single cell sample")
parser.add_argument("--paired_end", action="store_true",
                    help="The sample is paired-end (if this flag is not given, single end is assumed)")
parser.add_argument("-r", "--reference", action="store", required=True,
		    choices=["mm10", "hg38", "hg19"],
                    help="The referernce genome against which to align. Currently supported: mm10 = mm10, with ERCC spike-ins, RefSeq annotations, compiled by Allon.\nhg38 = human, compiled by Michael")                
parser.add_argument("-o", "--output_folder", action="store", required=True,
                    help="The directory to which output is written.")    
parser.add_argument("-p", "--num_threads", action="store", required=False, default=6,
                    help="The number of allocated threads, to be passed to trimmomatic, tophat, cufflinks, and rsem.")    
parser.add_argument('sampleFile1', action='store',
                   help='the first reads file')
parser.add_argument('sampleFile2', action='store', nargs = '?', default='',
                   help='the second reads file (for paired-end. If the --paired_end flag is not specified and this argument is specified then an exception is thrown. If the --paired_end flag is raised but this file is not specified then again an exception is thrown')
parser.add_argument('--skip_trimmomatic', action='store_true',
                   help="don't run trimmomatic before running on the samples")
parser.add_argument('--do_not_rely_on_previous_trimmomatic', action='store_true',
                   help="This flag has effect only if the flag --skip_trimmomatic is also set. The default behavior when --skip_trimmomatic is set is to rely on outputs from a previous trimmomatic run (it is assumed that they already exist, an error is thrown otherwise). However, if the flag --do_not_rely_on_previous_trimmomatic is set, then the program totally skips the trimmomatic phase and feeds the untrimmed reads to tophat/rsem");
parser.add_argument('--skip_tophat', action='store_true',
                   help="skip the tophat pipeline (note that you still have to set the --skip_tophat_qc flag separately if you wish)")
parser.add_argument('--skip_rsem', action='store_true',
                   help="skip the rsem pipeline (note that you still have to set the --skip_rsem_qc flag separately if you wish)")
parser.add_argument('--skip_qc', action='store_true',
                   help="skip the qc part of the pipeline")
parser.add_argument('--skip_tophat_qc', action='store_true',
                   help="skip the qc part of the pipeline only for tophat (ignored if the --skip_qc flag is given, in which case qc is not run in the first place)")                 
parser.add_argument('--skip_rsem_qc', action='store_true',
                   help="skip the qc part of the pipeline only for rsem (ignored if the --skip_qc flag is given, in which case qc is not run in the first place)")                 
                   
                                  
args = parser.parse_args()

if(not(args.paired_end) and args.sampleFile2 != ''):
	raise Exception("you cannot specify a second read file unless you raise the --paired_end flag");
elif(args.paired_end and args.sampleFile2 == ''):
	raise Exception("if the --paired_end flag is raised, then you must specify a second reads file");

if(args.reference == "mm10"):
	#settings for mm10 with ERCC spike-ins and other additions required by the BRAIN
	#project (eGFP, tdTomato, CreER)
	BOWTIE2_INDEX = "/data/yosef/index_files/mm10_4brain/index/GRCm38.p3_4brain";
	TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/mm10_4brain/index/tophat_transcriptome_data/GRCm38.p3_refseq_annot";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/mm10_4brain/index/GRCm38.p3.gff";
	RSEM_INDEX = "/data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem";
	#REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat/GRCm38.p3.refFlat";
	#REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat2/refFlat.txt";
	REF_FLAT_INDEX = "/data/yosef/index_files/mm10_4brain/index/refFlat2/refFlat_allonEdited.txt";
	#RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/gencode.vM4.rRNA.interval_list";
	#RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/gencode.vM4.rRNA.interval_list_allonEdited";
	RIBOSOMAL_INTERVALS_INDEX = "/data/yosef/index_files/mm10_4brain/index/gencode/mm10_4BRAIN_rRNA_interval_list_allonEdited.txt";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "/data/yosef/index_files/mm10_4brain/index/gencode/mm10_4BRAIN_rRNA_interval_list_allonEdited_forRSEM.txt";
		
		
	#settings for mm10 with ERCC spike-ins
	#BOWTIE2_INDEX = "/data/yosef/index_files/mm10_withERCC/GRCm38.p3_withERCC";
	#TOPHAT2_TRANSCRIPTOME_INDEX = "/data/yosef/index_files/mm10_withERCC/tophat_transcriptome_data/GRCm38.p3_refseq_annot";
	#TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/mm10_withERCC/GRCm38.p3.gff";
	#RSEM_INDEX = "/data/yosef/index_files/mm10_withERCC/rsem_index/GRCm38.p3_withERCC_rsem";
	
#hg19 was deprecated - then restored at Michael's request
elif(args.reference == "hg19"):
	#settings for hg19 with and other additions required by the HIV project (tba)
	#Compiled by Michael, Feb 2015
	BOWTIE2_INDEX = "/data/yosef/index_files/hg19_HIV/index/GRCh37.p13";
	TOPHAT2_TRANSCRIPTOME_INDEX = "";
	TRANSCRIPT_ANNOTATION = "";
	RSEM_INDEX = "/data/yosef/index_files/hg19_HIV/index/rsem_index/GRCh37.p13_rsem";
	REF_FLAT_INDEX = "/data/yosef/index_files/hg19_HIV/index/refFlat/refFlat.txt";
	RIBOSOMAL_INTERVALS_INDEX = "null";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM = "null";
	
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
	
	BOWTIE2_INDEX = "/data/yosef/index_files/hg38/index/GRCh38";
	TOPHAT2_TRANSCRIPTOME_INDEX ="/data/yosef/index_files/hg38/index/tophat_transcriptome_data/GRCh38";
	TRANSCRIPT_ANNOTATION = "/data/yosef/index_files/hg38/index/GRCh38.gtf";
	RSEM_INDEX ="/data/yosef/index_files/hg38/index/rsem_index/GRCh38_rsem";
	REF_FLAT_INDEX ="/data/yosef/index_files/hg38/index/refFlat/refFlat.txt";
	RIBOSOMAL_INTERVALS_INDEX ="/data/yosef/index_files/hg38/index/gencode/rRNA.interval";
	RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM ="/data/yosef/index_files/hg38/index/gencode/rRNA.rsem.interval";

		
else:
	raise Exception("should not happen - unsupported reference genome");


#sample file names should not have the '.gz' extension --> remove it if it is  present
if(args.sampleFile1[-3:] == '.gz'):
	print("Truncating the .gz extension from sample file 1's name");
	args.sampleFile1 = args.sampleFile1[0:-3];
if(args.paired_end and args.sampleFile2[-3:] == '.gz'):
	print("Truncating the .gz extension from sample file 2's name");
	args.sampleFile2 = args.sampleFile2[0:-3];
	
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
		
		trimMinLen = 16; #36
		trimCrop = "CROP:50";
		trimMinPhred = 15;


			
		trimCmd = Template("java -jar /opt/pkg/Trimmomatic-0.32/trimmomatic-0.32.jar $TRIM_IS_PAIRED -threads $NUM_THREADS -phred33 -trimlog $OUTPUT_FOLDER/trimmomatic_output/trimmomatic_log.txt $TRIM_INPUT1 $TRIM_INPUT2 $TRIM_OUTPUT LEADING:$TRIM_MIN_PHRED TRAILING:$TRIM_MIN_PHRED SLIDINGWINDOW:4:$TRIM_MIN_PHRED MINLEN:$TRIM_MINLEN $TRIM_CROP").substitute(TRIM_IS_PAIRED=trimIsPaired, OUTPUT_FOLDER=args.output_folder, TRIM_INPUT1=trimInput1, TRIM_INPUT2=trimInput2, TRIM_OUTPUT=trimOutput, TRIM_MINLEN=trimMinLen, TRIM_MIN_PHRED=trimMinPhred, TRIM_CROP=trimCrop, NUM_THREADS=args.num_threads);
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
	if(args.paired_end and (os.stat(args.sampleFile1).st_size == 0)):
		print "Danger, Will Robinson, sample file 2 was left empty after running trimmomatic!"
        	
        print("**********************************************************");
	print("**********************************************************");	
	print("trimmomatic done!");
        
        	
#I use the terms "cufflinks pipeline" and "tophat pipeline" interchangeably
RUN_CUFFLINKS_PIPELINE = True;
if(RUN_CUFFLINKS_PIPELINE and not(args.skip_tophat)):

	print("**********************************************************");
	print("**********************************************************");
	print("starting cufflinks pipeline");
	print("TODO (3/3/15): the tophat pipeline will output wrong number of reads and ratio of aligned reads in the summary,txt - this is the same problem that was in the rsem pipeline, the accepted_hits.bam will include multiple reads and they have to be filtered out before counting");
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
	picardDir = Template("$OUTPUT_FOLDER/tophat_output/picard_output").substitute(OUTPUT_FOLDER=args.output_folder);
	if not(os.path.exists(picardDir)):
		os.makedirs(picardDir);
		
	print("Sort SAM");
	sortSamCommand = Template("picard SortSam TMP_DIR=$OUTPUT_FOLDER/temp I=$OUTPUT_FOLDER/tophat_output/accepted_hits.bam O=$OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam SO=coordinate").substitute(OUTPUT_FOLDER=args.output_folder);
	print(sortSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(sortSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("sortSam failed");
		

		
	print("**********************************************************");
	print("**********************************************************");
	print("Index SAM");
	indexSamCommand = Template("/opt/genomics/bin/samtools index $OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(indexSamCommand)
	sys.stdout.flush();
	returnCode = subprocess.call(indexSamCommand, shell=True);
	if(returnCode != 0):
		raise Exception("indexSam failed");

	print("**********************************************************");
	print("**********************************************************");
	print("Running Cufflinks");
	cufflinksComand = Template("/opt/genomics/bin/cufflinks --num-threads $NUM_THREADS -G $TRANSCRIPT_ANNOTATION -o $OUTPUT_FOLDER/tophat_output/cuff_output/ $OUTPUT_FOLDER/tophat_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder, TRANSCRIPT_ANNOTATION=TRANSCRIPT_ANNOTATION, NUM_THREADS=args.num_threads);
	print(cufflinksComand)
	sys.stdout.flush();
	returnCode = subprocess.call(cufflinksComand, shell=True);
	if(returnCode != 0):
		raise Exception("cufflinks failed");
	
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

	#--sampling-for-bam: rsem outputs a bam with multiple possible alignments per read plus their posterior probabilities, which means that the read counts are disrupted
	#this flag tells it to sample one read per the distribution, assign it a MAPQ=100 value (also tag ZW:f:1) and the rest are 0 value (and tag ZW:f:0)
	
	print("Running rsem");
	#Note that in rsem I use the --output-genome-bam flag to generate a genome bam (in addition to the transcriptome bam that is always created) - rsem will also sort this bam. This is necessary for the QC later.
	if(args.paired_end):
		rsemComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --num-threads $NUM_THREADS --bowtie2 --estimate-rspd --output-genome-bam --sampling-for-bam --paired-end $SAMPLE_FILE1 $SAMPLE_FILE2 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2, NUM_THREADS=args.num_threads);
	else:
		rsemComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --num-threads $NUM_THREADS --bowtie2 --estimate-rspd --output-genome-bam --sampling-for-bam $SAMPLE_FILE1 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, NUM_THREADS=args.num_threads);
		
	print(rsemComand)
	sys.stdout.flush();
	returnCode = subprocess.call(rsemComand, shell=True);
	if(returnCode != 0):
		raise Exception("rsem failed");
		
	
		
	print("**********************************************************");
	print("**********************************************************");
	#to conform with the tophat pipeline, I split the reads into aligned and unaligned, and the QC will run only on the aligned ones.
	
	print("splitting rsem's output to aligned and unaligned reads...");
	splitCmd1 = Template("samtools view -b -F 4 $OUTPUT_FOLDER/rsem_output/rsem_output.genome.sorted.bam > $OUTPUT_FOLDER/rsem_output/accepted_hits.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd1)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd1, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");
		

	splitCmd2 = Template("samtools view -b -f 4 $OUTPUT_FOLDER/rsem_output/rsem_output.genome.sorted.bam > $OUTPUT_FOLDER/rsem_output/unmapped.bam").substitute(OUTPUT_FOLDER=args.output_folder);
	print(splitCmd2)
	sys.stdout.flush();
	returnCode = subprocess.call(splitCmd2, shell=True);
	if(returnCode != 0):
		raise Exception("split failed");
		
		
		
	print("removing the multiple-posterior alignments from rsem's aligned output reads...");
	#with the --output-genome-bam flag set, rsem outputs a bam in which the MAPQ are either 1 (chosen) or 0 (not chosen) 
	filterCmd = Template("samtools view -b -q 100 $OUTPUT_FOLDER/rsem_output/accepted_hits.bam > $OUTPUT_FOLDER/rsem_output/accepted_hits_noMultiple.bam").substitute(OUTPUT_FOLDER=args.output_folder);
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
	indexSamCommand = Template("/opt/genomics/bin/samtools index $OUTPUT_FOLDER/rsem_output/picard_output/sorted.bam").substitute(OUTPUT_FOLDER=args.output_folder);
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
	
	RUN_QC_TOPHAT = True;
	if(RUN_QC_TOPHAT and not(args.skip_tophat_qc)):
		sortedBamFile = args.output_folder + "/tophat_output/picard_output/sorted.bam";
		tophatOutputFolder = args.output_folder + "/tophat_output";
		tophat_qc_metrics = doQC.CollectData(sortedBamFile, tophatOutputFolder, REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX, args.paired_end);
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
		rsem_qc_metrics = doQC.CollectData(sortedBamFile, rsemOutputFolder, REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX_FOR_RSEM, args.paired_end);
		doQC.WriteQCMetrics(rsemOutputFolder, commonMetrics, rsem_qc_metrics);
	
	
	
	
	
	
	#if(args.paired_end):
	#	qcComand = Template("perl /project/eecs/yosef/singleCell/allon_script/preproc/qc.pl $QC_OUTPUT_DIR  $REF_FLAT_FILE $RIBOSOMAL_INTERVALS_INDEX $SAMPLE_FILE1 $SAMPLE_FILE2").substitute(QC_OUTPUT_DIR=args.output_folder, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2, REF_FLAT_FILE=REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX=RIBOSOMAL_INTERVALS_INDEX);
	#else:
	#	qcComand = Template("perl /project/eecs/yosef/singleCell/allon_script/preproc/qc.pl $QC_OUTPUT_DIR  $REF_FLAT_FILE $RIBOSOMAL_INTERVALS_INDEX $SAMPLE_FILE1").substitute(QC_OUTPUT_DIR=args.output_folder, SAMPLE_FILE1=args.sampleFile1, REF_FLAT_FILE=REF_FLAT_INDEX, RIBOSOMAL_INTERVALS_INDEX=RIBOSOMAL_INTERVALS_INDEX);
		
		
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