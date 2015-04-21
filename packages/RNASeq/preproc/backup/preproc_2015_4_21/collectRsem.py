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

from collections import namedtuple;
import argparse;
import os;
import numpy;
import shutil;

def CollectDupGenesTable(dup_filename, NUM_RSEM_GENES, NUM_CELLS, subDirectoryList):
	cellsWithDupGenesErrors = [];

	#note that I use object and not string, because numpy string array will by default be of size 1 byte (and they're of constant length anyway) so the string may get truncated
	dupTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS), dtype="object");

	print ("starting to collect dup reads (%s) per gene..." % dup_filename)
	for cell_ind in xrange(NUM_CELLS):
		d = subDirectoryList[cell_ind];
		fullDirPath = os.path.join(args.diretoryToProcess, d);
		print "Operating on single cell directory: " + fullDirPath;

		dupReadsPerGeneFile = os.path.join(fullDirPath, 'rsem_output/picard_output/%s.txt.genes.txt' % dup_filename);
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


parser = argparse.ArgumentParser(description="Collect RSEM output")
parser.add_argument("diretoryToProcess", action="store",
                    help="The folder to read (assumes the folder is a batch, with subdirectories for every cell)")
parser.add_argument("-r", "--reference", action="store", required=True,
		    choices=["mm10", "hg38"],
                    help="The referernce genome against which to align. Currently supported: mm10 = mm10, with ERCC spike-ins, RefSeq annotations, compiled by Allon.\nhg38 = human, compiled by Michael")                
parser.add_argument("-o", "--output_folder", action="store", required=True,
                    help="The directory to which output is written.")
parser.add_argument('--skip_collecting_expression', action='store_true',
	help="don't collect and overwrite the gene expression matrices")                 
parser.add_argument('--skip_collecting_qc', action='store_true',
	help="don't collect and overwrite the qc summary files")
parser.add_argument('--skip_collecting_dup_genes', action='store_true',
	help="don't collect and overwrite the dup genes summary files")
          
args = parser.parse_args();

if(args.reference == "mm10"):
		rsemDictionaryFile = "/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
elif(args.reference == "hg38"):	
		rsemDictionaryFile="/data/yosef/index_files/hg38/index/rsem_dict.txt";
else:
	raise Exception("should not happen - unsupported reference genome");

if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder);
        


#read the mapping I prepared from gene ID to gene name and gene Nir type
GeneRecord = namedtuple('GeneRecord', 'geneID, geneName, geneNirType')

with open(rsemDictionaryFile) as fin:
    # skip the 3 header lines
    for i in xrange(3):
        next(fin)

    rows = ( line.split('\t') for line in fin )
    rsemAnnotations = map(GeneRecord._make, rows)

rsemGeneList = [record.geneID for record in rsemAnnotations];
NUM_RSEM_GENES = len(rsemGeneList);

#a dictionary from a geneID to the index in the list (which is also the index in the matrix I'm going  to create) - to be used when reading the individual files
mapGeneIDToInd = {rsemGeneList[i]:i for i in xrange(NUM_RSEM_GENES)};


GeneResultRecord = namedtuple('GeneResultRecord', 'geneID, transcriptIDs, length,  effective_length, expected_count, TPM, FPKM');


#os.listdir lists files and directories, but not the '.' and '..' entries
#for compatability with Nir's previous structure, I do not allow the cell directory name to be "cuff" or "rsem" (this is where Nir used to save the collect results of the folder)
subDirectoryList = [x for x in os.listdir(args.diretoryToProcess) \
	if (os.path.isdir(os.path.join(args.diretoryToProcess, x))) and (x != "cuff") and (x != "rsem")];
NUM_CELLS = len(subDirectoryList);

COLLECT_GENE_EXPRESSION = True;
if (COLLECT_GENE_EXPRESSION and not(args.skip_collecting_expression)):
	cellsWithGeneExpressionErrors = [];

	readCountsTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
	tpmTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
	fpkmTable = numpy.empty((NUM_RSEM_GENES, NUM_CELLS));
	readCountsTable[:] = numpy.NAN;
	tpmTable[:] = numpy.NAN;
	fpkmTable[:] = numpy.NAN;

	print "starting to collect gene expression data..."
	for cell_ind in xrange(NUM_CELLS):
		d = subDirectoryList[cell_ind];
		fullDirPath = os.path.join(args.diretoryToProcess, d);
		print "Operating on single cell directory: " + fullDirPath;
		
		rsemOutputGenes = os.path.join(fullDirPath, 'rsem_output/rsem_output.genes.results');
		
		if(os.path.exists(rsemOutputGenes)):
			with open(rsemOutputGenes) as fin:
				next(fin); #skip the first (header) line
				rows = ( line.split('\t') for line in fin )
				cellGeneExpressionTable = map(GeneResultRecord._make, rows)
				
			#cellGeneList = [record.geneID for record in cellGeneExpressionTable];
			
			for record in cellGeneExpressionTable:
				indInMatrix = mapGeneIDToInd[record.geneID];
				readCountsTable[indInMatrix, cell_ind] = float(record.expected_count)
				tpmTable[indInMatrix, cell_ind] = float(record.TPM);
				fpkmTable[indInMatrix, cell_ind] = float(record.FPKM);
		else:
			#the rsem gene results for this folder does not exist for some reason
			
			#the values were initialized to NaN so seemingly no reason to overwrite again, but I do this just to be on the safe side
			readCountsTable[:, cell_ind] = numpy.NAN;
			tpmTable[:, cell_ind] = numpy.NAN;
			fpkmTable[:, cell_ind] = numpy.NAN;
			
			cellsWithGeneExpressionErrors.append(rsemOutputGenes);
			

	print "finished collecting data... writing results..."
	numpy.savetxt(os.path.join(args.output_folder, 'rsem_readCountsTable.txt'), readCountsTable, delimiter="\t", fmt="%g");
	numpy.savetxt(os.path.join(args.output_folder, 'rsem_tpmTable.txt'), tpmTable, delimiter="\t", fmt="%g");
	numpy.savetxt(os.path.join(args.output_folder, 'rsem_fpkmTable.txt'), fpkmTable, delimiter="\t", fmt="%g");




COLLECT_QC_SUMMARIES = True;
if (COLLECT_QC_SUMMARIES and not(args.skip_collecting_qc)):
	cellsWithQCErrors = [];
	QC_FIELDS = ["NREADS", "NALIGNED", "RALIGN", "TOTAL_DUP", "PRIMER", "INSERT_SZ", "INSERT_SZ_STD", "COMPLEXITY", "NDUPR", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS"];
	NUM_QC_FIELDS = len(QC_FIELDS);
	
	qcTable = numpy.empty((NUM_QC_FIELDS, NUM_CELLS));
	qcTable[:] = numpy.NAN;
	
	print "starting to collect QC metrics..."
	for cell_ind in xrange(NUM_CELLS):
		d = subDirectoryList[cell_ind];
		fullDirPath = os.path.join(args.diretoryToProcess, d);
		print "Operating on single cell directory: " + fullDirPath;
		
		rsemQCSummaryFile = os.path.join(fullDirPath, 'rsem_output/summary.txt');
		if(os.path.exists(rsemQCSummaryFile)):
			with open(rsemQCSummaryFile) as fin:
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
			
			cellsWithQCErrors.append(rsemQCSummaryFile);


	#now, write the collected QC results into a file
	print "finished collecting qc summaries... writing qc table..."
	numpy.savetxt(os.path.join(args.output_folder, 'qc_table.txt'), qcTable, delimiter="\t", fmt="%g");


print "writing cell list..."

#write the ordered cell list file...
with open(os.path.join(args.output_folder, 'cell_list.txt'), 'wt') as fout:
	 fout.writelines(subdir + '\n' for subdir in subDirectoryList);

print "writing gene dictionary..."

#now just copy the gene list that I used, that is actually rsem's dictionary that I prepared in advance from the index
#-- for consistency I simply use the same gene list, in the same order, which is based on the dictionary I built from the annotation file rsem got

shutil.copyfile(rsemDictionaryFile, os.path.join(args.output_folder, 'gene_list.txt'));



#Added 3/9/2015 for the BRAIN project qc assessment:
COLLECT_DUP_GENES = True;
if (COLLECT_DUP_GENES and not(args.skip_collecting_dup_genes)):

	dupTable, cellsWithDupGenesErrors = CollectDupGenesTable('dup', NUM_RSEM_GENES, NUM_CELLS, subDirectoryList)
	dupUniqueTable, cellsWithUniqueDupGenesErrors = CollectDupGenesTable('dup_unique', NUM_RSEM_GENES, NUM_CELLS, subDirectoryList)

	#now, write the collected QC results into a file
	print "finished collecting duplicate reads per gene files... writing summary table..."
	numpy.savetxt(os.path.join(args.output_folder, 'dup_reads_per_gene_table.txt'), dupTable, delimiter="\t", fmt="%s");
	numpy.savetxt(os.path.join(args.output_folder, 'dup_reads_per_gene_onlyUnique_table.txt'), dupTable, delimiter="\t", fmt="%s");




print "done collecting all data!"

if (COLLECT_GENE_EXPRESSION and not(args.skip_collecting_expression)):
	if(not(cellsWithGeneExpressionErrors)):
		print "no cells had errors while collecting their gene expression";
	else:
		print "the following cells had errors while collecting their gene expression:";
		for x in cellsWithGeneExpressionErrors:
			print x;


if (COLLECT_QC_SUMMARIES and not(args.skip_collecting_qc)):
	if(not(cellsWithQCErrors)):
		print "no cells had errors while collecting their qc";
	else:
		print "the following cells had errors while collecting their qc:";
		for x in cellsWithQCErrors:
			print x;


if (COLLECT_DUP_GENES and not(args.skip_collecting_dup_genes)):
	if(not(cellsWithDupGenesErrors)):
		print "no cells had errors while collecting their dup genes";
	else:
		print "the following cells had errors while collecting their dup genes:";
		for x in cellsWithDupGenesErrors:
			print x;

		print "the following cells had errors while collecting their dup genes (only unique reads):";
		for x in cellsWithUniqueDupGenesErrors:
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

    
