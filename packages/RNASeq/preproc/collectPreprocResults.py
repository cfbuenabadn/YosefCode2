#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015

#pipeline todo:
# 1. collapse gene variants in rsem's collect output (as Bo recommended)
# 2. Better handling of multiple alignments:
# a. RSEM currently randomizes one location to which the read is mapped from the posterior distribution of the reads, it's better to use that probability
# throughout the computation for counting reads, duplicates etc.
# b. In cufflinks the output includes multiple reads per gene and the number of reads and number aligned reads will be wrong (there is a warning about this in the output)
# 3. update the list of protein coding genes (cox2 is not there currently, as Michael pointed out)
# 4. The problem Michael pointed to --- rRNA and VDJ have no exons and so gffread will not put them in the gtf it produces, which means they are not in rsem's dictionary and no reads can map to them
# 5. Fix the bug in countdup - it doesn't work on the Ngai single end for some reason
# 6. debug the ribosomal % index - is it ok? Do we see so few rRNAs because we use RefSeq that has much fewer rRNAs than GenCode, or is there a bug in the index?
# 7. add a capability to count_dup_per_gene for cufflinks (currently it relies on rsem indices...)
# 8. Input unaligned BAM file as Nir requested
# 9. For count_dup_unique: uses too much memory. We can increase the $GRID parameter. It should be a ~linear trade off between memory and running time. Maybe parallelization?
# 10. add a parameter for "declared number of processors" to allow using a different number of processors than was declared... this will allow qsub to put more jobs per machine
# 11. scripts for counting reads: htseq and featureCounts http://genomespot.blogspot.com/2014/09/read-counting-with-featurecounts.html
# 12. support Kallisto
# 13. Broad's qc: http://www.broadinstitute.org/cancer/cga/rnaseqc_run
# 14. Use Prinseq as another way to trim? (supports trimming low complexity reads,and other goodies see http://www.pnas.org/content/suppl/2015/05/14/1507125112.DCSupplemental/pnas.1507125112.sapp.pdf)
# 15. Make sample name meaningful in processFolder: instead of the -N for the project, do -N for the sample name
# 16. add gene coverage, uniformity etc. with bedtools, Shaked suggested the GINI coefficient as a measure of read uniformity dispersion across the transcript
# 17. Cufflinks does not quantify the ERCC and the CreER/eGFP genes because I did not add them to the gff annotation! (I did add them to the gtf annotation that rsem uses...). I added them manually to the cufflinks dictionary as a preparation, but I still have to update the gff annotation that cufflinks uses...

#and make the argparse parent: did this but still replicates command line arguments when writing the PBS script.
#I should have used the: argparse.REMAINDER. All the remaining command-line arguments are gathered into a list. This is commonly useful for command line utilities that dispatch to other command line utilities:
#BUT then I have to give the remainder at the end, and cannot use them upfront....
#maybe convert the parsed namespace object back to a string and omit the arguments unique to processFolder?
#https://docs.python.org/3/library/argparse.html
#https://docs.python.org/3/library/argparse.html#the-parse-args-method

#call Nir's new count dup script


from collections import namedtuple;
import argparse;
import os;
import numpy;
import shutil;

def CollectDupGenesTable(dup_filename, NUM_GENES, NUM_CELLS, subDirectoryList):
	cellsWithDupGenesErrors = [];

	#note that I use object and not string, because numpy string array will by default be of size 1 byte (and they're of constant length anyway) so the string may get truncated
	dupTable = numpy.empty((NUM_GENES, NUM_CELLS), dtype="object");

	print ("starting to collect dup reads (%s) per gene..." % os.path.basename(dup_filename))
	for cell_ind in xrange(NUM_CELLS):
		fullDirPath = subDirectoryList[cell_ind];
		print "Operating on single cell (%d / %d) directory: %s" % (cell_ind+1, NUM_CELLS, fullDirPath);

		dupReadsPerGeneFile = os.path.join(fullDirPath, dup_filename);
		if(os.path.exists(dupReadsPerGeneFile)):
			with open(dupReadsPerGeneFile) as fin:
				rows = ( line.split('\t') for line in fin )
				#I assume that the gene list in the dup.txt.genes.txt (file with dup reads per gene file) is the same as in the rsem dictionary - this is how I built it when Nir and I programmed it
				#take from each line the number of dups  and the total number of reads and put them in a string in the format "%f/%f" - this way one file contains both numbers and the ratio can be easily computed
				dupReadsPerGeneTable = [("%s/%s" % (row[1], row[2])) for row in rows];

			dupTable[:, cell_ind] = dupReadsPerGeneTable;
		else:
			#the rsem dup reads per gene file for this folder does not exist for some reason

			#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
			dupTable[:, cell_ind] = ("%f/%f" % (numpy.NAN, numpy.NAN));

			cellsWithDupGenesErrors.append(dupReadsPerGeneFile);

	return dupTable, cellsWithDupGenesErrors

def CollectQcTable(qc_filename, NUM_CELLS, subDirectoryList):
	cellsWithQCErrors = [];
	QC_FIELDS = ["NREADS", "NALIGNED", "RALIGN", "TOTAL_DUP", "PRIMER", "INSERT_SZ", "INSERT_SZ_STD", "COMPLEXITY", "NDUPR", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS"];
	NUM_QC_FIELDS = len(QC_FIELDS);

	qcTable = numpy.empty((NUM_QC_FIELDS, NUM_CELLS));
	qcTable[:] = numpy.NAN;

	print "starting to collect QC metrics..."
	for cell_ind in xrange(NUM_CELLS):
		fullDirPath = subDirectoryList[cell_ind];
		print "Operating on single cell (%d / %d) directory: %s" % (cell_ind+1, NUM_CELLS, fullDirPath);

		qcSummaryFile = os.path.join(fullDirPath, qc_filename);
		if(os.path.exists(qcSummaryFile)):
			with open(qcSummaryFile) as fin:
				rows = ( line.split('\t') for line in fin )
				cellQCTable = [(row[0], row[1]) for row in rows];

			#I first store the pair and then separate into two lists because using the generator burns it and it cannot be rewound (I could have backed it up first with itertools.tee but that seemed like an overkill)
			cellQCTableEntries = [row[0] for row in cellQCTable];
			cellQCTableValues = [float(row[1]) for row in cellQCTable];


			if(cmp(cellQCTableEntries, QC_FIELDS) != 0):
				print cellQCTableEntries;
				print "******************";
				print QC_FIELDS;
				raise Exception("the QC entries of cell " + d + " do not match the expected QC entries - this should not happen!");

			#if we reached this point, then there's a perfect match between the expected qc fields for the table, and the ones present in the summary file
			qcTable[:, cell_ind] = cellQCTableValues;

		else:
			#the rsem qc results for this folder does not exist for some reason

			#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
			qcTable[:, cell_ind] = numpy.NAN;

			cellsWithQCErrors.append(fullDirPath);

	return qcTable, cellsWithQCErrors

parser = argparse.ArgumentParser(description="Collect output from the preproc pipeline")
parser.add_argument("directoriesToProcess", action="store",
                    help="The folders to read (assumes each folder is a batch, with subdirectories for every cell). Folder are semi-colon (';') separated ")
parser.add_argument("-r", "--reference", action="store", required=True,
		    choices=["mm10", "hg38"],
                    help="The referernce genome against which to align. Currently supported: mm10 = mm10, with ERCC spike-ins, RefSeq annotations, compiled by Allon.\nhg38 = human, compiled by Michael")                
parser.add_argument("-o", "--output_folder", action="store", required=True, default="",
                    help="The directory to which output is written (if not specified: add a '/rsem' folder to the input directory's name)")
parser.add_argument('--skip_collecting_expression', action='store_true',
	help="don't collect and overwrite the gene expression matrices")                 
parser.add_argument('--skip_collecting_qc', action='store_true',
	help="don't collect and overwrite the qc summary files")
parser.add_argument('--skip_collecting_dup_genes', action='store_true',
	help="don't collect and overwrite the dup genes summary files")
parser.add_argument('--skip_collecting_feature_counts', action='store_true',
	help="don't collect and overwrite the feature counts (read counts per gene made with the featureCounts software) files")
parser.add_argument('--skip_collecting_rsem', action='store_true',
	help="don't collect any of the data produced by the rsem pipeline")
parser.add_argument('--skip_collecting_tophat', action='store_true',
	help="don't collect any of the data produced by the tophat pipeline")
args = parser.parse_args();


if(args.reference == "mm10"):
	rsemDictionaryFile = "/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
	cuffDictionaryFile = "/home/eecs/allonwag/data/index_files/mm10_4brain/index/cuffDictionary/mm10_4brain_cuffGeneMapping_noLoci.txt"
elif(args.reference == "hg38"):	
	rsemDictionaryFile="/data/yosef/index_files/hg38/index/rsem_dict.txt";
else:
	raise Exception("should not happen - unsupported reference genome");

############################
### Preparations for reading
############################

args.output_folder = os.path.expanduser(args.output_folder);
if not os.path.exists(args.output_folder):
	os.makedirs(args.output_folder);

if not(args.skip_collecting_rsem):
	if not os.path.exists(os.path.join(args.output_folder, 'rsem')):
		os.makedirs(os.path.join(args.output_folder, 'rsem'))

	#read the mapping I prepared from gene ID to gene name and gene biotype
	rsemGeneRecord = namedtuple('rsemGeneRecord', 'geneID, geneName, geneBiotype')

	with open(rsemDictionaryFile) as fin:
		# skip the 3 header lines
		for i in xrange(3):
			next(fin)

		rows = ( line.strip().split('\t') for line in fin )
		rsemAnnotations = map(rsemGeneRecord._make, rows)

	rsemGeneList = [record.geneID for record in rsemAnnotations];
	NUM_RSEM_GENES = len(rsemGeneList);

	#a dictionary from a geneID to the index in the list (which is also the index in the matrix I'm going  to create) - to be used when reading the individual files
	rsem_mapGeneIDToInd = {rsemGeneList[i]:i for i in xrange(NUM_RSEM_GENES)};



if not(args.skip_collecting_tophat):
	if not os.path.exists(os.path.join(args.output_folder, 'cuff')):
		os.makedirs(os.path.join(args.output_folder, 'cuff'))


	#read the mapping I prepared from gene ID to gene name and gene biotype
	cuffGeneRecord = namedtuple('cuffGeneRecord', 'geneID, geneName, geneBiotype')

	with open(cuffDictionaryFile) as fin:
		# skip the 3 header lines
		for i in xrange(3):
			next(fin)

		rows = ( line.strip().split('\t') for line in fin )
		cuffAnnotations = map(cuffGeneRecord._make, rows)

	cuffGeneList = [record.geneID for record in cuffAnnotations];
	NUM_CUFF_GENES = len(cuffGeneList);


	#a dictionary from a geneID to the index in the list (which is also the index in the matrix I'm going  to create) - to be used when reading the individual files
	cuff_mapGeneIDToInd = {cuffGeneList[i]:i for i in xrange(NUM_CUFF_GENES)};


rsemGeneResultRecord = namedtuple('rsemGeneResultRecord', 'geneID, transcriptIDs, length,  effective_length, expected_count, TPM, FPKM');
cuffGeneResultRecord = namedtuple('cuffGeneResultRecord', 'geneID, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, locus, length, coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status');

#strip leading and trailing semicolons if they're present in the semicolon separated string so as not to confuse the split
args.directoriesToProcess = args.directoriesToProcess.strip(';')
#os.path.expanduser --> if the path begins with a tilde - expand it to the user's homedir
directoriesToProcess = [os.path.expanduser(d) for d in args.directoriesToProcess.split(';')]

subDirectoryList = []
for d in directoriesToProcess:
	#os.listdir lists files and directories, but not the '.' and '..' entries
	#for compatability with Nir's previous structure, I do not allow the cell directory name to be "cuff" or "rsem" (this is where Nir used to save the collect results of the folder)
	curSubDirectoryList = [os.path.join(d, x) for x in os.listdir(d) \
		if (os.path.isdir(os.path.join(d, x))) and (x != "cuff") and (x != "rsem")];

	subDirectoryList += curSubDirectoryList

NUM_CELLS = len(subDirectoryList);


############################
### Collect gene expression
############################


COLLECT_GENE_EXPRESSION = True;
if (COLLECT_GENE_EXPRESSION and not(args.skip_collecting_expression)):
	rsem_cellsWithGeneExpressionErrors = [];
	tophat_cellsWithGeneExpressionErrors = [];

	if not(args.skip_collecting_rsem):
		readCountsTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
		tpmTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
		fpkmTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
		readCountsTable[:] = numpy.NAN;
		tpmTable[:] = numpy.NAN;
		fpkmTable[:] = numpy.NAN;

		print "starting to collect rsem gene expression data..."
		for cell_ind in xrange(NUM_CELLS):
			fullDirPath = subDirectoryList[cell_ind];
			print "Operating on single cell (%d / %d) directory: %s" % (cell_ind+1, NUM_CELLS, fullDirPath);

			rsemOutputGenes = os.path.join(fullDirPath, 'rsem_output/rsem_output.genes.results');

			if(os.path.exists(rsemOutputGenes)):
				with open(rsemOutputGenes) as fin:
					next(fin); #skip the first (header) line
					rows = ( line.split('\t') for line in fin )
					cellGeneExpressionTable = map(rsemGeneResultRecord._make, rows)

				#cellGeneList = [record.geneID for record in cellGeneExpressionTable];

				for record in cellGeneExpressionTable:
					indInMatrix = rsem_mapGeneIDToInd[record.geneID];
					readCountsTable[indInMatrix, cell_ind] = float(record.expected_count)
					tpmTable[indInMatrix, cell_ind] = float(record.TPM);
					fpkmTable[indInMatrix, cell_ind] = float(record.FPKM);
			else:
				#the rsem gene results for this folder does not exist for some reason

				#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
				readCountsTable[:, cell_ind] = numpy.NAN;
				tpmTable[:, cell_ind] = numpy.NAN;
				fpkmTable[:, cell_ind] = numpy.NAN;

				rsem_cellsWithGeneExpressionErrors.append(rsemOutputGenes);


		print "finished collecting rsem data... writing results..."
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/rsem_readCountsTable.txt'), readCountsTable, delimiter="\t", fmt="%g");
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/rsem_tpmTable.txt'), tpmTable, delimiter="\t", fmt="%g");
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/rsem_fpkmTable.txt'), fpkmTable, delimiter="\t", fmt="%g");



	if not(args.skip_collecting_tophat):
		fpkmTable = numpy.empty((NUM_CUFF_GENES, NUM_CELLS));
		fpkmTable[:] = numpy.NAN;

		print "starting to collect cufflinks gene expression data..."
		for cell_ind in xrange(NUM_CELLS):
			fullDirPath = subDirectoryList[cell_ind];
			print "Operating on single cell (%d / %d) directory: %s" % (cell_ind+1, NUM_CELLS, fullDirPath);

			cuffOutputGenes = os.path.join(fullDirPath, 'tophat_output/cuff_output/genes.fpkm_tracking');

			if(os.path.exists(cuffOutputGenes)):
				with open(cuffOutputGenes) as fin:
					next(fin); #skip the first (header) line
					rows = ( line.strip().split('\t') for line in fin )
					cellGeneExpressionTable = map(cuffGeneResultRecord._make, rows)

				for record in cellGeneExpressionTable:
					indInMatrix = cuff_mapGeneIDToInd[record.geneID];
					#make sure that the annotation is consistent
					if record.gene_short_name.upper() != cuffAnnotations[indInMatrix].geneName.upper():
						raise Exception("unexpected mismatch between the cufflinks dictionary and the current cufflinks results file")

					if(record.FPKM_status == "OK"):
						fpkmTable[indInMatrix, cell_ind] = float(record.FPKM);
					else:
						fpkmTable[indInMatrix, cell_ind] = numpy.NAN;
			else:
				#the cufflinks gene results for this folder does not exist for some reason

				#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
				fpkmTable[:, cell_ind] = numpy.NAN;

				tophat_cellsWithGeneExpressionErrors.append(cuffOutputGenes);


		print "finished collecting cufflinks data... writing results..."
		numpy.savetxt(os.path.join(args.output_folder, 'cuff/cuff_fpkmTable.txt'), fpkmTable, delimiter="\t", fmt="%g");


############################
### Write cell list
############################

print "writing cell list..."

#write the ordered cell list file...
#identify each cell by the last two dirs in the directory path, which are the batch name (parent dir) and the cell name (son dir)
with open(os.path.join(args.output_folder, 'rsem/cell_list.txt'), 'wt') as fout:
	for d in subDirectoryList:
		#strip trailing slashes so as to not confuse the os.path.split command
		cur_dir = d.rstrip('/')
		parentDir, sonDir = os.path.split(cur_dir)
		#take only the last element in the parent path, which is the direct parent dir of the cell
		batchName = os.path.basename(parentDir)

		#cell name is composed of the batch ID and the cellID within it
		cellName = batchName + '/' + sonDir

		fout.write(cellName + '\n')

#it is the same list of cells between the cufflinks and the rsem pipelines - for convenience I write it twice...
shutil.copyfile(os.path.join(args.output_folder, 'rsem/cell_list.txt'), os.path.join(args.output_folder, 'cuff/cell_list.txt'),);

############################
### Write gene lists
############################

print "writing gene dictionary..."

#now just copy the gene list that I used, that is actually rsem's dictionary that I prepared in advance from the index
#-- for consistency I simply use the same gene list, in the same order, which is based on the dictionary I built from the annotation file rsem got

shutil.copyfile(rsemDictionaryFile, os.path.join(args.output_folder, 'rsem/gene_list.txt'));

#same deal for the cufflinks pipeline...
shutil.copyfile(cuffDictionaryFile, os.path.join(args.output_folder, 'cuff/gene_list.txt'));


############################
### Collect QC
############################


COLLECT_QC_SUMMARIES = True;
if (COLLECT_QC_SUMMARIES and not(args.skip_collecting_qc)):

	if(not(args.skip_collecting_rsem)):
		rsem_qcTable, rsem_cellsWithQCErrors = CollectQcTable('rsem_output/summary.txt', NUM_CELLS, subDirectoryList)
		print "finished collecting rsem qc summaries... writing rsem qc table..."
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/qc_table.txt'), rsem_qcTable, delimiter="\t", fmt="%g");

	if(not(args.skip_collecting_tophat)):
		tophat_qcTable, tophat_cellsWithQCErrors = CollectQcTable('tophat_output/summary.txt', NUM_CELLS, subDirectoryList)
		print "finished collecting tophat qc summaries... writing tophat qc table..."
		numpy.savetxt(os.path.join(args.output_folder, 'cuff/qc_table.txt'), tophat_qcTable, delimiter="\t", fmt="%g");



############################
### Collect dup genes
############################


#Added 3/9/2015 for the BRAIN project qc assessment:
COLLECT_DUP_GENES = True;
if (COLLECT_DUP_GENES and not(args.skip_collecting_dup_genes)):

	if not(args.skip_collecting_rsem):
		rsem_dupTable, rsem_cellsWithDupGenesErrors = CollectDupGenesTable('rsem_output/picard_output/dup.txt.genes.txt', NUM_RSEM_GENES, NUM_CELLS, subDirectoryList)
		rsem_dupUniqueTable, rsem_cellsWithUniqueDupGenesErrors = CollectDupGenesTable('rsem_output/picard_output/dup_unique.txt.genes.txt', NUM_RSEM_GENES, NUM_CELLS, subDirectoryList)

		#now, write the collected QC results into a file
		print "finished collecting duplicate reads per gene files in rsem pipeline... writing summary table..."
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/dup_reads_per_gene_table.txt'), rsem_dupTable, delimiter="\t", fmt="%s");
		numpy.savetxt(os.path.join(args.output_folder, 'rsem/dup_reads_per_gene_onlyUnique_table.txt'), rsem_dupUniqueTable, delimiter="\t", fmt="%s");
	else:
		rsem_cellsWithDupGenesErrors = []
		rsem_cellsWithUniqueDupGenesErrors = []

	if not(args.skip_collecting_tophat):
		tophat_dupTable, tophat_cellsWithDupGenesErrors = CollectDupGenesTable('tophat_output/picard_output/dup.txt.genes.txt', NUM_CUFF_GENES, NUM_CELLS, subDirectoryList)
		tophat_dupUniqueTable, tophat_cellsWithUniqueDupGenesErrors = CollectDupGenesTable('tophat_output/picard_output/dup_unique.txt.genes.txt', NUM_CUFF_GENES, NUM_CELLS, subDirectoryList)

		#now, write the collected QC results into a file
		print "finished collecting duplicate reads per gene files in tophat pipeline... writing summary table..."
		numpy.savetxt(os.path.join(args.output_folder, 'cuff/dup_reads_per_gene_table.txt'), tophat_dupTable, delimiter="\t", fmt="%s");
		numpy.savetxt(os.path.join(args.output_folder, 'cuff/dup_reads_per_gene_onlyUnique_table.txt'), tophat_dupUniqueTable, delimiter="\t", fmt="%s");
	else:
		tophat_cellsWithDupGenesErrors = []
		tophat_cellsWithUniqueDupGenesErrors = []


COLLECT_FEATURE_COUNTS = True
if(COLLECT_FEATURE_COUNTS and not(args.skip_collecting_feature_counts)):
	#*********************************
	#featureCounts are produced on the tophat2 results (i.e., as part of the tophat2-cufflinks pipeline) and not as part of the rsem pipeline
	#(because I found that rsem's alignments do not have the NH flag in the sam file, which does not allow featureCounts to identify
	#the multiply-aligned reads).
	#*********************************

	FeatureCountsResultRecord = namedtuple('FeatureCountsResultRecord', 'geneID, chr, start, end, strand, length, count');

	cellsWithFeatureCountsErrors = []
	featureCountsTable = numpy.empty((NUM_CUFF_GENES, NUM_CELLS));
	featureCountsTable[:] = numpy.NAN;

	print "starting to collect featureCounts..."
	for cell_ind in xrange(NUM_CELLS):
		fullDirPath = subDirectoryList[cell_ind];
		print "Operating on single cell (%d / %d) directory: %s" % (cell_ind+1, NUM_CELLS, fullDirPath);

		featureCountsForGenes = os.path.join(fullDirPath, 'tophat_output/feature_counts.txt');

		if(os.path.exists(featureCountsForGenes)):
			with open(featureCountsForGenes) as fin:
				next(fin); #skip the two first (header) lines
				next(fin);

				rows = ( line.split('\t') for line in fin )
				cellFeatureCountTable = map(FeatureCountsResultRecord._make, rows)

			##############
			#INTERESTING: at first I programmed the collect feature count while using the rsem dictionary - and everything worked fine no genes appeared int the feature count!
			#but then when replacing the rsem dictionary by the (seemingly more appropriate) cuff dictionary, I found that there are genes that do not appear in the cufflinks dictionary
			#(which is the output format of all the cufflinks quantifications) but appear in featureCount! Probably genes that cufflinks does not think it's worth quantifying  in the first place...
			##############

			for record in cellFeatureCountTable:
				if cuff_mapGeneIDToInd.has_key(record.geneID):
					indInMatrix = cuff_mapGeneIDToInd[record.geneID];
					featureCountsTable[indInMatrix, cell_ind] = float(record.count)
				else:
					#see comment above. This else is never reached when using the rsem dictionary, but can be reached when using the cufflinks dictionary,
					#IMPORTANT: I decided to throw counts data for genes that cufflinks does not quantify
					#THESE are MOSTLY genes with 0 counts but NOT always!
					indInRsemMatrix = rsem_mapGeneIDToInd[record.geneID];
					geneSymbol = rsemAnnotations[indInRsemMatrix].geneName
					print "gene %s = %s appears in featureCounts but not in the cufflinks dictionary (count = %f)\n" % (record.geneID, geneSymbol, float(record.count))

		else:
			#the featureCounts results for this folder do not exist for some reason

			#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
			featureCountsTable[:, cell_ind] = numpy.NAN;
			cellsWithFeatureCountsErrors.append(featureCountsForGenes);


	print "finished collecting feature count data... writing results..."
	numpy.savetxt(os.path.join(args.output_folder, 'cuff/tophat2_featureCountsTable.txt'), featureCountsTable, delimiter="\t", fmt="%g");



print "done collecting all data!"

if (COLLECT_GENE_EXPRESSION and not(args.skip_collecting_expression)):
	if(not(rsem_cellsWithGeneExpressionErrors)):
		print "no cells had errors while collecting their rsem gene expression";
	else:
		print "the following cells had errors while collecting their rsem gene expression:";
		for x in rsem_cellsWithGeneExpressionErrors:
			print x;

	if(not(tophat_cellsWithGeneExpressionErrors)):
		print "no cells had errors while collecting their cufflinks gene expression";
	else:
		print "the following cells had errors while collecting their cufflinks gene expression:";
		for x in tophat_cellsWithGeneExpressionErrors:
			print x;


if (COLLECT_QC_SUMMARIES and not(args.skip_collecting_qc)):

	if(not(args.skip_collecting_rsem)):
		if(not(rsem_cellsWithQCErrors)):
			print "no cells had errors while collecting their rsem qc";
		else:
			print "the following cells had errors while collecting their rsem qc:";
			for x in rsem_cellsWithQCErrors:
				print x;


	if(not(args.skip_collecting_tophat)):
		if(not(tophat_cellsWithQCErrors)):
			print "no cells had errors while collecting their tophat qc";
		else:
			print "the following cells had errors while collecting their tophat qc:";
			for x in tophat_cellsWithQCErrors:
				print x;

if (COLLECT_DUP_GENES and not(args.skip_collecting_dup_genes)):
	if(not(args.skip_collecting_rsem)):
		if(not(rsem_cellsWithDupGenesErrors) and not(rsem_cellsWithUniqueDupGenesErrors)):
			print "no cells had errors while collecting their dup genes (rsem pipeline)";
		else:
			print "the following cells had errors while collecting their dup gene (rsem pipeline)s:";
			for x in rsem_cellsWithDupGenesErrors:
				print x;

			print "the following cells had errors while collecting their dup genes (only unique reads) (rsem pipeline):";
			for x in rsem_cellsWithUniqueDupGenesErrors:
				print x;


	if(not(args.skip_collecting_tophat)):
		if(not(tophat_cellsWithDupGenesErrors) and not(tophat_cellsWithUniqueDupGenesErrors)):
			print "no cells had errors while collecting their dup genes (tophat pipeline)";
		else:
			print "the following cells had errors while collecting their dup genes (tophat pipeline):";
			for x in tophat_cellsWithDupGenesErrors:
				print x;

			print "the following cells had errors while collecting their dup genes (only unique reads)  (tophat pipeline):";
			for x in tophat_cellsWithUniqueDupGenesErrors:
				print x;

	#turns out I cannot use that... There are indeed differences between the lists
	#if(cmp(cellGeneList, rsemGeneList) != 0):
	#	print len(cellGeneList);
	#	print "*****";
	#	print len(rsemGeneList);
	#	print [aa for aa in rsemGeneList if aa not in cellGeneList];
	#	print "*****";
	#	print [aa for aa in cellGeneList if aa not in rsemGeneList];
	#
	#	raise Exception("the gene list of cell " + d + " does not match the general annotation - this should not happen!");
	#
	#readCountsTable[:, cell_ind] = [float(record.expected_count) for record in cellGeneExpressionTable];
	#tpmTable[:, cell_ind] = [float(record.TPM) for record in cellGeneExpressionTable];
	#fpkmTable[:, cell_ind] = [float(record.FPKM) for record in cellGeneExpressionTable];


if (COLLECT_FEATURE_COUNTS and not(args.skip_collecting_feature_counts)):
	if(not(cellsWithFeatureCountsErrors)):
		print "no cells had errors while collecting their featureCounts";
	else:
		print "the following cells had errors while collecting their featureCounts:";
		for x in cellsWithFeatureCountsErrors:
			print x;