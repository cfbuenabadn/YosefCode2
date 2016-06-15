# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Mar 2015

import argparse

#A parser that will have the common arguments of the processSingleSample and the processFolder scripts
common_rnaseq_parser = argparse.ArgumentParser(add_help=False) #from https://docs.python.org/3/library/argparse.html#parents: Note that most parent parsers will specify add_help=False. Otherwise, the ArgumentParser will see two -h/--help options (one in the parent and one in the child) and raise an error.
common_rnaseq_parser.add_argument("--paired_end", action="store_true",
                    help="The sample is paired-end (if this flag is not given, single end is assumed)")
common_rnaseq_parser.add_argument("-r", "--reference", action="store", required=True,
		    choices=["mm10", "hg38", "hg19", "hg38_HIV", "z10"],
                    help="The referernce genome against which to align. Currently supported: mm10 = mm10, with ERCC spike-ins, RefSeq annotations, compiled by Allon.\nhg38 = human, compiled by Michael\nz10 = Zebrafish, GRCz10, compiled by Allon")
common_rnaseq_parser.add_argument("-o", "--output_folder", action="store", required=True,
                    help="The directory to which output is written.")
common_rnaseq_parser.add_argument("-p", "--num_threads", action="store", required=False, default=1,
                    help="The number of allocated threads, to be passed to trimmomatic, tophat, cufflinks, and rsem.")
common_rnaseq_parser.add_argument('--skip_trimmomatic', action='store_true',
                   help="don't run trimmomatic before running on the samples")
common_rnaseq_parser.add_argument('--do_not_rely_on_previous_trimmomatic', action='store_true',
                   help="This flag has effect only if the flag --skip_trimmomatic is also set. The default behavior when --skip_trimmomatic is set is to rely on outputs from a previous trimmomatic run (it is assumed that they already exist, an error is thrown otherwise). However, if the flag --do_not_rely_on_previous_trimmomatic is set, then the program totally skips the trimmomatic phase and feeds the untrimmed reads to tophat/rsem");
common_rnaseq_parser.add_argument('--skip_tophat', action='store_true',
                   help="skip the tophat pipeline (note that you still have to set the --skip_tophat_qc flag separately if you wish)")
common_rnaseq_parser.add_argument('--skip_rsem', action='store_true',
                   help="skip the rsem pipeline (note that you still have to set the --skip_rsem_qc flag separately if you wish)")
common_rnaseq_parser.add_argument('--skip_kallisto', action='store_true',
                   help="skip the Kallisto pipeline (note that you still have to set the --skip_kallisto_qc flag separately if you wish)")
common_rnaseq_parser.add_argument('--skip_qc', action='store_true',
                   help="skip the qc part of the pipeline")
common_rnaseq_parser.add_argument('--skip_tophat_qc', action='store_true',
                   help="skip the qc part of the pipeline only for tophat (ignored if the --skip_qc flag is given, in which case qc is not run in the first place)")
common_rnaseq_parser.add_argument('--skip_rsem_qc', action='store_true',
                   help="skip the qc part of the pipeline only for rsem (ignored if the --skip_qc flag is given, in which case qc is not run in the first place)")
common_rnaseq_parser.add_argument('--skip_kallisto_qc', action='store_true',
                   help="skip the qc part of the pipeline only for Kallisto (ignored if the --skip_qc flag is given, in which case qc is not run in the first place)")
common_rnaseq_parser.add_argument('--do_not_clean_intermediary_files', action='store_true',
                   help="If set, do not clean intermediary files that are produced in the course of running (default: off, i.e. clean the intermediary files)")
common_rnaseq_parser.add_argument('--rsem_bowtie_maxins', action='store', default=1000,
                   help="For paired-end data only (ignored if --paired_end is not set): the maximum fragment length (this is the value of the --fragment-length-max in rsem and -X/--maxins in bowtie2). Defaults to 1000, which is the rsem default")
common_rnaseq_parser.add_argument('--rsem_samtools_sort_mem', action='store', default="1G",
                   help="This is rsem's samtools-sort-mem argument. It controls the memory allocated per thread when samtools sort is run through rsem. Defaults to 1G, which is the rsem default")
common_rnaseq_parser.add_argument('--trimmomatic_window', action='store', default='',
                   help="The trimmomatic sliding window argument. Format: '<windowSize>:<requiredQuality>' ")
common_rnaseq_parser.add_argument('--kallisto_bootstrap_samples', action='store', default='0',
                   help="The number of bootstraps done by Kallisto (default: 0)")
common_rnaseq_parser.add_argument('--mean_fragment_length', action='store', default='200',
                   help="The mean fragment length paramter that Kallisto requires; applies only to single-end; this parameter will be ignored and the fragment length estimated from the data if the --paired_end flag is set (default: 200, as in cufflinks). If given, will be used also in cufflinks (again, used only in SE and ignored in PE).")
common_rnaseq_parser.add_argument('--std_fragment_length', action='store', default='80',
                   help="The standard deviation of fragment length paramter that Kallisto requires; applies only to single-end; this parameter will be ignored and the fragment length estimated from the data if the --paired_end flag is set (default: 80, as in cufflinks). If given, should have been used also in cufflinks (again, used only in SE and ignored in PE), but since it leads to core dumps it is ignored in cufflinks. The std of fragment length in cufflinks is always 80, which is cufflink's default")

#For insert size vs. fragment length see https://www.biostars.org/p/95803/

def GetFirstReadInFastqFile(fastqFilename):
    #a very simple implementation that does not rely on biopython and therefore supports only fastq files and not fasta file
    print "note that the current implementation of this function supports only fastq format... replace it with biopython-based or support fasta explicitly"
    with open(fastqFilename, 'r') as f:
        first_line = f.readline();
        first_read = f.readline();

    if(first_read[-1] != '\n'):
        raise Exception("unexpected - there should be a newline here...");

    #strip the
    first_read = first_read.rstrip('\n');

    return first_read;