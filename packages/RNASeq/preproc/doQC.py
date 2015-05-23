#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015

import subprocess
from string import Template
#import argparse
import sys
import os;
import glob;

def RunFastQC(sampleFile, outputFolder):
	cmd = Template("/opt/genomics/bin/fastqc $f -o  $OUTPUT_FOLDER/fastqc_output").substitute(f=sampleFile, OUTPUT_FOLDER=outputFolder);
			
	print(cmd)
	returnCode = subprocess.call(cmd, shell=True);
	if(returnCode != 0):
		raise Exception("fastqc failed");


def CheckContam(sampleFile, outputFolder, sample_end):
	running_folder = os.path.realpath(os.path.dirname(sys.argv[0]))
	scriptToCall = os.path.join(running_folder, "check_contam.pl")
	cmd = Template("perl $SCRIPT_TO_CALL $f  $OUTPUT_FOLDER/fastqc_output/primer.$END.txt").substitute(SCRIPT_TO_CALL=scriptToCall,f=sampleFile, OUTPUT_FOLDER=outputFolder, END=sample_end);

	print(cmd)
	returnCode = subprocess.call(cmd, shell=True);
	if(returnCode != 0):
		raise Exception("CheckContam failed");


def CollectData(bamFile, outputFolder, refFlatAnnotationsFile, ribosomalIntervalsFile, isPairedEnd, transcriptAnnotationFile, transcriptDictionaryFile, genomeReferenceFile):
	rnaMetricsFileName = outputFolder + '/picard_output/rna_metrics.txt';
	cmd = Template("picard CollectRnaSeqMetrics TMP_DIR=$OUTPUT_FOLDER/temp INPUT=$BAM_FILE OUTPUT=$OUTPUT_FILE CHART=$OUTPUT_FOLDER/picard_output/rna_coverage.pdf REF_FLAT=$refFlatFile STRAND=NONE RIBOSOMAL_INTERVALS=$ribosomalIntervalsFile").substitute(BAM_FILE=bamFile, OUTPUT_FOLDER=outputFolder, OUTPUT_FILE=rnaMetricsFileName, refFlatFile=refFlatAnnotationsFile, ribosomalIntervalsFile=ribosomalIntervalsFile);
	print(cmd)
	returnCode = subprocess.call(cmd, shell=True);
	if(returnCode != 0):
		raise Exception("CollectData failed");

	alnMetricsFileName = outputFolder + '/picard_output/aln_metrics.txt';
	#I found in an online forum (http://sourceforge.net/p/samtools/mailman/message/32772099/) that you have to give the genome reference file for this command to work, even though all the data it needs may already be contained in the bam file...
	cmd = Template("picard CollectAlignmentSummaryMetrics INPUT=$BAM_FILE OUTPUT=$OUTPUT_FILE REFERENCE_SEQUENCE=$GENOME_REFERENCE_FILE\n").substitute(BAM_FILE=bamFile, OUTPUT_FILE=alnMetricsFileName, GENOME_REFERENCE_FILE=genomeReferenceFile);
	print(cmd)
	returnCode = subprocess.call(cmd, shell=True);
	if(returnCode != 0):
		raise Exception("CollectData failed");

	insertMetricsFileName = outputFolder + '/picard_output/aln_metrics1.txt';
	cmd = Template("picard CollectInsertSizeMetrics INPUT=$BAM_FILE OUTPUT=$OUTPUT_FILE H=$OUTPUT_FOLDER/picard_output/ins_metrics.histogram.pdf").substitute(BAM_FILE=bamFile, OUTPUT_FILE=insertMetricsFileName, OUTPUT_FOLDER=outputFolder);
	print(cmd)
	returnCode = subprocess.call(cmd, shell=True);
	if(returnCode != 0):
		raise Exception("CollectData failed");

	if(not(transcriptAnnotationFile) or not(transcriptDictionaryFile)):
		print "the index file required to run the count dup logic were not supplied... skipping this part and not creating dup.txt and dup.txt.genes.txt"
	else:
		dupFileName = outputFolder + '/picard_output/dup.txt';

		USE_ONLY_OLD_SCRIPT = True
		if(USE_ONLY_OLD_SCRIPT):
			#use old script until the new one is fixed:
			dupScriptToRun = "count_dup.pl"
			running_folder = os.path.realpath(os.path.dirname(sys.argv[0]))
			dupScriptToRun = os.path.join(running_folder, dupScriptToRun)
			cmd = Template("perl $DUP_SCRIPT_TO_RUN $BAM_FILE $OUTPUT_FILE").substitute(DUP_SCRIPT_TO_RUN=dupScriptToRun, BAM_FILE=bamFile, OUTPUT_FILE=dupFileName);
		else:
			if(isPairedEnd):
				dupScriptToRun = "count_dup_per_gene_nextgen.pl"#"count_dup_per_gene.pl" # "count_dup.pl"
			else:
				dupScriptToRun = "count_dup_per_gene_nextgen_single.pl"#"count_dup_per_gene_single_end.pl"

			running_folder = os.path.realpath(os.path.dirname(sys.argv[0]))
			dupScriptToRun = os.path.join(running_folder, dupScriptToRun)

			#updated script that in addition to total dup counting, also counts them per gene: It creates the dup.txt file as before and also a new file with a postfix of .genes.txt that includes: <gene name> <#dup reads> <tot reads> <ratio> Note that this is a read-level analysis, not fragment.
			cmd = Template("perl $DUP_SCRIPT_TO_RUN $BAM_FILE $OUTPUT_FILE $TRANSCRIPT_ANNOTATION $TRANSCRIPT_DICTIONARY 0").substitute(BAM_FILE=bamFile, OUTPUT_FILE=dupFileName, TRANSCRIPT_ANNOTATION=transcriptAnnotationFile, TRANSCRIPT_DICTIONARY=transcriptDictionaryFile, DUP_SCRIPT_TO_RUN=dupScriptToRun);


		print(cmd)
		returnCode = subprocess.call(cmd, shell=True);
		if(returnCode != 0):
			raise Exception("count dups failed");

		if(False):
			#in addition to the previous call to the script which counted dups, now call it in a way that counts unique dups
			dupUniqueFileName = outputFolder + '/picard_output/dup_unique.txt';
			cmd = Template("perl $DUP_SCRIPT_TO_RUN $BAM_FILE $OUTPUT_FILE $TRANSCRIPT_ANNOTATION $TRANSCRIPT_DICTIONARY 1").substitute(BAM_FILE=bamFile, OUTPUT_FILE=dupUniqueFileName, TRANSCRIPT_ANNOTATION=transcriptAnnotationFile, TRANSCRIPT_DICTIONARY=transcriptDictionaryFile, DUP_SCRIPT_TO_RUN=dupScriptToRun);
			print(cmd)
			returnCode = subprocess.call(cmd, shell=True);
			if(returnCode != 0):
				raise Exception("count dups unique failed");


		#for debug
		if(False):
			#in addition, call the original old script that take LOTS of time to compare the results

			if(isPairedEnd):
				oldDupScriptToRun = "count_dup_per_gene.pl"
			else:
				oldDupScriptToRun = "count_dup_per_gene_single_end.pl"

			oldDupScriptToRun = os.path.join(running_folder, oldDupScriptToRun)


			dupOldFileName = outputFolder + '/picard_output/dup_OLD.txt';
			cmd = Template("perl $OLD_DUP_SCRIPT_TO_RUN $BAM_FILE $OUTPUT_FILE $TRANSCRIPT_ANNOTATION $TRANSCRIPT_DICTIONARY 0").substitute(OLD_DUP_SCRIPT_TO_RUN=oldDupScriptToRun, BAM_FILE=bamFile, OUTPUT_FILE=dupOldFileName, TRANSCRIPT_ANNOTATION=transcriptAnnotationFile, TRANSCRIPT_DICTIONARY=transcriptDictionaryFile, DUP_SCRIPT_TO_RUN=dupScriptToRun);
			print(cmd)
			returnCode = subprocess.call(cmd, shell=True);
			if(returnCode != 0):
				raise Exception("old count dups failed");

			#in addition to the previous call to the script which counted dups, now call it in a way that counts unique dups
			dupOldUniqueFileName = outputFolder + '/picard_output/dup_OLD_unique.txt';
			cmd = Template("perl $OLD_DUP_SCRIPT_TO_RUN $BAM_FILE $OUTPUT_FILE $TRANSCRIPT_ANNOTATION $TRANSCRIPT_DICTIONARY 1").substitute(OLD_DUP_SCRIPT_TO_RUN=oldDupScriptToRun, BAM_FILE=bamFile, OUTPUT_FILE=dupOldUniqueFileName, TRANSCRIPT_ANNOTATION=transcriptAnnotationFile, TRANSCRIPT_DICTIONARY=transcriptDictionaryFile, DUP_SCRIPT_TO_RUN=dupScriptToRun);
			print(cmd)
			returnCode = subprocess.call(cmd, shell=True);
			if(returnCode != 0):
				raise Exception("old count dups unique failed");

		
	#system("/opt/genomics/bin/CollectRnaSeqMetrics TMP_DIR=$OUTPUT_FOLDER/temp INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/rna_metrics.txt CHART=$OUTPUT_FOLDER/picard_output/rna_coverage.pdf REF_FLAT=$refFlatFile STRAND=NONE RIBOSOMAL_INTERVALS=$ribosomalIntervalsFile\n");
    #system("/opt/genomics/bin/CollectAlignmentSummaryMetrics INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/aln_metrics.txt\n");
    #system("/opt/genomics/bin/CollectInsertSizeMetrics INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/aln_metrics1.txt H=$OUTPUT_FOLDER/picard_output/ins_metrics.histogram.pdf\n");
    #system("perl count_dup.pl $OUTPUT_FOLDER/picard_output/sorted.bam $OUTPUT_FOLDER/picard_output/dup.txt\n");

	#  Total number of reads
	#  % of aligned reads	
	with open(alnMetricsFileName) as fin:
		# skip the header lines
		for _ in xrange(6):
			a = next(fin)

		rows = ( line.split('\t') for line in fin );

		firstRow = rows.next();
		secondRow = rows.next();
	    
		if(secondRow[0] != "UNPAIRED"):
			#if this is not an unpaired file, then the 2nd row is the 1st of the pair, the 3rd row is the 2nd of the pair, and we want the 4th rou which is the pair itself
			secondRow = rows.next(); #read (and discard) 3rd row
			secondRow = rows.next(); #read 4th row
		
			if(secondRow[0] != "PAIR"):
				raise Exception("unrecognized structure of alignment metrics file");


	if(len(firstRow) != len(secondRow)):
		raise Exception("length of dictionary rows in alignment metrics file are not the same?!");

	qc_metrics = {};
	#for i in xrange(1, 22): #start from 1 because col 0 is the category column (which is a string: FIRST_OF_PAIR\SECOND_OF_PAIR\PAIR in paired end and UNPAIRED for single end). Similarly columns 22-24 are string values)
	#    qc_metrics[firstRow[i]] = float(secondRow[i]);
	
	#take only the number of aligned reads and the total number of reads, which are of interest from this file	
	nreads = float(secondRow[2]); #NOTE: I take PF reads as the total number of reads (i.e., reads passing Illumina filter) because this seems to be the basis to almost all metrics that are produced by AlignmentSummaryMetrics - there's no need to consider reads that don't pass the Illumina filter). See: http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics
	nalign = float(secondRow[5]);
	r_align = 100 * float(nalign) / nreads;
	
	qc_metrics["nreads"] = nreads;
	qc_metrics["nalign"] = nalign;
	qc_metrics["r_align"] = r_align;
	
	if(True):
		#NOTE: there seems to be a bug in Picard's output that always gives 0 as the number of aligned reads and ratio of aligned reads (maybe because of inconsistencies in the sam format between Picard and Samtools?) So I compute these manually
		try:
			#if this is a paired-end sample, take only the first of the pair (so as not to count all reads twice and get a number twice as large as you should)
			forPairedEnd = ("-f 64" if isPairedEnd else "");
		
			nalign = float(subprocess.check_output([Template("samtools view -c $FOR_PAIRED_END $OUTPUT_FOLDER/accepted_hits_noMultiple.bam").substitute(OUTPUT_FOLDER=outputFolder, FOR_PAIRED_END=forPairedEnd)], shell=True));
			n_unalign = float(subprocess.check_output([Template("samtools view -c $FOR_PAIRED_END $OUTPUT_FOLDER/unmapped.bam").substitute(OUTPUT_FOLDER=outputFolder, FOR_PAIRED_END=forPairedEnd)], shell=True));
			nreads = nalign + n_unalign;
			r_align = 100 * float(nalign) / nreads;
			
			qc_metrics["nreads"] = nreads;
			qc_metrics["nalign"] = nalign;
			qc_metrics["r_align"] = r_align;
		except:
			raise Exception("failed to count alignments in bam files");
			
	#makes sense only for paired end data
	#  MEDIAN_INSERT_SIZE
	#  MEDIAN_ABSOLUTE_DEVIATION		
	if(isPairedEnd):
		with open(insertMetricsFileName) as fin:
			for line in fin:
				if line.startswith("## METRICS CLASS"):
					#skip the headers line
					fin.next();
					dataLine = fin.next().split('\t');
					insert_sz_avg = float(dataLine[0]);
					insert_sz_std = float(dataLine[1]);
					break;
					
	else:
		insert_sz_avg = float('NaN');
		insert_sz_std = float('NaN');
		
	qc_metrics["insert_sz_avg"] = insert_sz_avg;
	qc_metrics["insert_sz_std"] = insert_sz_std;

	if(not(transcriptAnnotationFile) or not(transcriptDictionaryFile)):
		print "the index file required to run the count dup logic were not supplied... skipping this part and not creating or reading data from dup.txt and dup.txt.genes.txt --> these values will be NaN"
		complexity = float('NaN');
		ndupr = float('NaN');
		left_right = float('NaN');
	else:
		#  % of unique alignments (complexity)
		with open(dupFileName) as fin:
			line = fin.next().split('\t');
		left_right = float(line[0]);
		dupf = float(line[1]); #dup ratio for fragments (pair end)

		if(dupf == -1):
			#this is what count_dup.pl outputs when fcounter = 0
			dupf = float('NaN');
		complexity = 1.0 - dupf;
		dupr = float(line[4]); #dup ratio for reads (one end)
		if(dupr == -1):
			dupr = float('NaN');
		ndupr = 1 - dupr; #for compatability with complexity, I output the 1-complement of the dup read ratio, i.e., the non-dup read ratio, or the ratio of unique reads

	qc_metrics["complexity"] = complexity;
	qc_metrics["ndupr"] = ndupr;
	qc_metrics["left_right"] = left_right;
	
	
	#  PCT_RIBOSOMAL_BASES
	#  PCT_CODING_BASES
	#  PCT_UTR_BASES
	#  INTRONIC_BASES
	#  PCT_INTERGENIC_BASES
	#  PCT_MRNA_BASES
	#  MEDIAN_CV_COVERAGE (evenness)
	#  MEDIAN_5PRIME_BIAS
	#  MEDIAN_3PRIME_BIAS
	#  MEDIAN_5PRIME_TO_3PRIME_BIAS
	with open(rnaMetricsFileName) as fin:
		for line in fin:
			if line.startswith("## METRICS CLASS"):
				headerLine = fin.next().split('\t');
				dataLine = fin.next().split('\t');
				break;

	if(len(headerLine) != len(dataLine)):
		raise Exception("data line and header line do not have the same length?!");
	
	for i in (range(10, 16) + range(18, 22)): #in Nir's perl code it is: $i(10..15,18..21)
		if(headerLine[i] in qc_metrics) and \
			(qc_metrics[headerLine[i]] != float(dataLine[i])):
			
			errorMsg = Template("key $KEY already exists in qc_metrics and their values don't match (!$VAL1 != $VAL2)").substitute(KEY=headerLine[i], VAL1=qc_metrics[headerLine[i]], VAL2=dataLine[i]);
			raise Exception(errorMsg);
		
		if((not dataLine[i]) or (dataLine[i] == "?")):
			#the string is empty - the value is not there (for example, if a ribosomal index file was not given to CollectRnaSeqMetrics, the RIBOSOMAL_BASES field will be  empty)
			qc_metrics[headerLine[i]] = float('NaN');
		else:
			qc_metrics[headerLine[i]] = float(dataLine[i]);
		
	return qc_metrics;
		
def CollectFastQcData(output_folder, isPairedEnd):
	#same as in doQC.RunFastQC
	fastQcFolder = output_folder + '/fastqc_output'; 
	fastQcFiles = glob.glob(fastQcFolder + "/*_fastqc/fastqc_data.txt");
	
	print "fastQcFiles are: " 
	print fastQcFiles;
	if(isPairedEnd and len(fastQcFiles) != 2) or \
		(not(isPairedEnd) and len(fastQcFiles) != 1):
			raise Exception("wrong number of FastQC files?!");
	
	#  Sequence duplication level (Total)
	totalDup = 0.0;	
	for fileName in fastQcFiles:
		with open(fileName) as fin:
			for line in fin:
				if line.startswith("#Total Duplicate Percentage"):
					totalDup += float(line.partition('\t')[2]);
	if(isPairedEnd):
		#in paired-end we added total-dup from two fastqc results, so we average them
		totalDup /= 2.0; 
	
	
	# Abundance of primer sequences
	primerFile1 = output_folder + '/fastqc_output/primer.1.txt';
	nhits, nseq = CollectPrimerData(primerFile1);
	if(isPairedEnd):
		primerFile2 = output_folder + '/fastqc_output/primer.2.txt';
		nhits2, nseq2 = CollectPrimerData(primerFile2);
		nhits += nhits2;
		nseq += nseq2;
		
	primer = float(nhits) / nseq;
	
	#these are metrics that are common to the tophat and rsem pipelines
	commonMetrics = {};
	commonMetrics["totalDup"] = totalDup;
	commonMetrics["nhits"] = nhits;
	commonMetrics["nseq"] = nseq;
	commonMetrics["primer"] = primer;
	
	return commonMetrics;


def CollectPrimerData(primerFileName):
	with open(primerFileName) as fin:
		line = next(fin).split('\t')
	nhits = float(line[1]);
	nseq = float(line[2]);
	return nhits, nseq;
	
	
def WriteQCMetrics(outputFolder, commonMetrics, qc_metrics):
	#commonMetrics - ones produced by the fastqc program directly from the reads and are common to both the rsem and the tophat pipelines
	#qc_metrics - produced from the sorted.bam file, and are specific to either the  tophat or the rsem pipelines
	
	outputFileName = outputFolder + '/summary.txt';
	with open(outputFileName, 'wt') as fout:
		fout.write(Template("NREADS\t$nreads\nNALIGNED\t$nalign\nRALIGN\t$r_align\n").substitute(nreads=qc_metrics["nreads"], nalign=qc_metrics["nalign"], r_align=qc_metrics["r_align"]));
		fout.write(Template("TOTAL_DUP\t$total_dup\nPRIMER\t$primer\nINSERT_SZ\t$insert_sz_avg\nINSERT_SZ_STD\t$insert_sz_std\nCOMPLEXITY\t$complexity\nNDUPR\t$ndupr\n").substitute(total_dup=commonMetrics["totalDup"], primer=commonMetrics["primer"], insert_sz_avg=qc_metrics["insert_sz_avg"], insert_sz_std=qc_metrics["insert_sz_std"], complexity=qc_metrics["complexity"], ndupr=qc_metrics["ndupr"]));

		keysToWrite=("PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS");
		for key in keysToWrite:
			if not(key in qc_metrics):
				raise Exception("cannot find qc metric value " + key + "!");
				
			fout.write(Template("$key\t$value\n").substitute(key=key, value=qc_metrics[key]));

#parser = argparse.ArgumentParser(description="Do QC for one single cell sample")
#parser.add_argument("--paired_end", action="store_true",
#                    help="The sample is paired-end (if this flag is not given, single end is assumed)")
#parser.add_argument("--refFlatAnnotationsFile", action="store", required=True,
#		    help="Gene annotations in refFlat format, to be used by Picard's CollectRnaSeqMetrics")                
#parser.add_argument("-o", "--output_folder", action="store", required=True,
#                    help="The directory to which output is written.")    
#parser.add_argument('bamFile', action='store',
#                   help='sorted bam file of aligned reads')
#parser.add_argument('sampleFile2', action='store', nargs = '?', default='',
#                   help='the second reads file (for paired-end. If the --paired_end flag is not specified and this argument is specified then an exception is thrown. If the --paired_end flag is raised but this file is not specified then again an exception is thrown')
#              
#args = parser.parse_args()

