#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015

# NOTE: this code assumes that PATH and PYTHONPATH environment variables were set in the .bashrc file


import subprocess
from string import Template
import argparse
import os
import glob
import stat
import socket
import sys
import rnaSeqPipelineUtils


parser = argparse.ArgumentParser(description="Process all the samples in a folder", parents=[rnaSeqPipelineUtils.common_rnaseq_parser])
parser.add_argument('folder', action='store',
                   help='the folder to read')
parser.add_argument("-N", "--job_name", action="store", default='RnaSeq_preproc_pipeline',
                    help="A name for the job submitted to the cluster.")
parser.add_argument('-d', '--do_not_send_to_cluster', action='store_true',
	help="just write the directory structure, but do not call the cluster to execute the jobs (used for debug)")

args = parser.parse_args();

#if the path begins with a tilde - expand it to the user's homedir
args.folder = os.path.expanduser(args.folder);
args.output_folder = os.path.expanduser(args.output_folder);

if(args.folder[-1] == '/'):
	#I assume the folder does not end with the '/' sign - if this is the last char  in the  folder string, then remove it
	args.folder = args.folder[:-1];


if(socket.gethostname() != "zen" and not(args.do_not_send_to_cluster)):
	raise Exception("This script must be run from the zen machine that supports the  queue");

#get the path in which the current script is running
scriptDirectory = os.path.realpath(os.path.dirname(sys.argv[0]))
scriptToRun = os.path.join(scriptDirectory, "processSingleSample.py")
print "Script to execute on each cell is: " + scriptToRun


#get the command line arguments for the processSingleSample script. I prefer not to use argparse.REMAINDER
#(which is normally recommended for command line utilities that dispatch to other command line utilities) because then
#the redundant arguments must be grouped at the end of the command
# argumentsForScriptAsDictionary = vars(args)
# #remove the arguments that are unique to processFolder and are not consumed by processSingleSample:
# keysToExclude = ["folder", "job_name", "do_not_send_to_cluster"]
# for key in keysToExclude:
# 	argumentsForScriptAsDictionary.pop(key, None) #pop is a convenient way to remove from a dictionary only if key exists
# argumentsForScript = ' '.join(["--%s %s" % (key, value) for key,value in vars(args).iteritems()])
# print "Command line arguments for the script to execute on each cell are: " + argumentsForScript
# This way is also problematic because what do you do with the flag arguments? They should be "--flag" and not "--flag True"


print "Running on folder...";
print "Folder is: " + args.folder;

#folder_format = Template("$BASE_FOLDER/*/*_R1_001.fastq.gz").substitute(BASE_FOLDER=args.folder);
#folder_format = Template("$BASE_FOLDER/*/*.R1_L001.fastq.gz").substitute(BASE_FOLDER=args.folder);
folder_format = Template("$BASE_FOLDER/*/*_R1_combined.fastq.gz").substitute(BASE_FOLDER=args.folder);
#to take the other end
#folder_format = Template("$BASE_FOLDER/*/*_R2_combined.fastq.gz").substitute(BASE_FOLDER=args.folder);

sampleList = glob.iglob(folder_format);

for sample1 in sampleList:
	print "processing sample: " + sample1;
	
	if(args.paired_end):
		#sample2 = str.replace(sample1, "_R1_001.fastq.gz", "_R2_001.fastq.gz");
		sample2 = str.replace(sample1, "_R1_combined.fastq.gz", "_R2_combined.fastq.gz");
	else:
		sample2 = "";
		
	sampleName = os.path.basename(sample1).replace("_R1_combined.fastq.gz", "");
	sampleOutputFolder = os.path.join(args.output_folder, sampleName);
	
	#remove this redundant ending if it's in the dir name
	if(sampleOutputFolder.endswith(".fastq.gz")):
		sampleOutputFolder = sampleOutputFolder[:-len(".fastq.gz")];
	
	if not os.path.exists(sampleOutputFolder):
		os.makedirs(sampleOutputFolder);

	cmd = Template("python $SCRIPT_TO_RUN $IS_PAIRED_END \
					-r $REFERENCE -p $NUM_THREADS -o $OUTPUT_FOLDER $SKIP_TRIMMOMATIC $DO_NOT_RELY_ON_PREVIOUS_TRIMMOMATIC \
					$SKIP_TOPHAT $SKIP_RSEM $SKIP_KALLISTO $SKIP_QC $SKIP_TOPHAT_QC $SKIP_RSEM_QC $SKIP_KALLISTO_QC \
					$DO_NOT_CLEAN_INTERMEDIARY_FILES $RSEM_BOWTIE_MAXINS $RSEM_SAMTOOLS_SORT_MEM $TRIMMOMATIC_WINDOW \
					$KALLISTO_BOOTSTRAP_SAMPLES $MEAN_FRAGMENT_LENGTH $STD_FRAGMENT_LENGTH $SAMPLE1 $SAMPLE2").substitute(
		SCRIPT_TO_RUN=scriptToRun,
		IS_PAIRED_END=("--paired_end" if args.paired_end else ""),
		REFERENCE=args.reference,
		NUM_THREADS=args.num_threads,
		OUTPUT_FOLDER=sampleOutputFolder,
		SKIP_TRIMMOMATIC="--skip_trimmomatic" if args.skip_trimmomatic else "",
		DO_NOT_RELY_ON_PREVIOUS_TRIMMOMATIC="--do_not_rely_on_previous_trimmomatic" if args.do_not_rely_on_previous_trimmomatic else "",
		SKIP_TOPHAT="--skip_tophat" if args.skip_tophat else "",
		SKIP_RSEM="--skip_rsem" if args.skip_rsem else "",
		SKIP_KALLISTO="--skip_kallisto" if args.skip_kallisto else "",
		SKIP_QC="--skip_qc" if args.skip_qc else "",
		SKIP_TOPHAT_QC="--skip_tophat_qc" if args.skip_tophat_qc else "",
		SKIP_RSEM_QC="--skip_rsem_qc" if args.skip_rsem_qc else "",
		SKIP_KALLISTO_QC="--skip_kallisto_qc" if args.skip_kallisto_qc else "",
		DO_NOT_CLEAN_INTERMEDIARY_FILES="--do_not_clean_intermediary_files" if args.do_not_clean_intermediary_files else "",
		RSEM_BOWTIE_MAXINS=("--rsem_bowtie_maxins %s" % args.rsem_bowtie_maxins) if args.paired_end else "",
		RSEM_SAMTOOLS_SORT_MEM=("--rsem_samtools_sort_mem %s" % args.rsem_samtools_sort_mem) if args.rsem_samtools_sort_mem else "",
		TRIMMOMATIC_WINDOW=("--trimmomatic_window %s" % args.trimmomatic_window) if args.trimmomatic_window else "",
		KALLISTO_BOOTSTRAP_SAMPLES=("--kallisto_bootstrap_samples %s" % args.kallisto_bootstrap_samples) if args.kallisto_bootstrap_samples else "",
		MEAN_FRAGMENT_LENGTH=("--mean_fragment_length %s" % args.mean_fragment_length) if args.mean_fragment_length else "",
		STD_FRAGMENT_LENGTH=("--std_fragment_length %s" % args.std_fragment_length) if args.std_fragment_length else "",
		SAMPLE1=sample1,
		SAMPLE2=sample2 if args.paired_end else "")


	#a simple way to remove all duplicate whitespaces and replace them with one whitespace. The duplicate whitespaced occur because of the way I implement not transferring optional arguments (it leaves extra double whitespaces where the optional arg could have been)
	cmd = ' '.join(cmd.split())

	print "sending cmd to queue: " + cmd;

	print("**********************************************************");
	print("**********************************************************");
	print("writing PBS script\n");
	

	queueScriptFileName = os.path.join(sampleOutputFolder, "run.sh")
	with open(queueScriptFileName, "wt") as fout:
		#see https://www.millennium.berkeley.edu/wiki/torque_example
		#another useful tutorial: https://wikis.nyu.edu/display/NYUHPC/Tutorial+-+Submitting+a+job+using+qsub
		
		#IMPORTANT NOTE: All of the ### and #PBS lines forming the job header need to be collected at the top of the batch file and *not* interspersed with blank lines or shell command lines. qsub seems to only properly parses the top of the file.
		fout.write("#!/bin/sh\n");
		fout.write("# *** \"#PBS\" lines must come before any non-blank, non-comment lines ***")
		fout.write("### Set the job name\n");
		fout.write(Template("#PBS -N $JOB_NAME\n").substitute(JOB_NAME=args.job_name));
		
		### Redirect stdout and stderr by first telling Torque to redirect do /dev/null
		### and then redirecting yourself via exec. This is the way the IT recommends.
		fout.write("### Redirect stdout and stderr by first telling Torque to redirect do /dev/null and then redirecting yourself via exec. This is the way the IT recommends.\n");
		fout.write("#PBS -e localhost:/dev/null\n");
		fout.write("#PBS -o localhost:/dev/null\n");
		
		fout.write("### Set the queue to which to send\n");
		fout.write("#PBS -q yosef\n");
		
		fout.write("### Limit the resources used\n");
		fout.write(Template("#PBS -l nodes=1:ppn=$NUM_THREADS\n").substitute(NUM_THREADS=args.num_threads));
		
		fout.write("### Change the walltime and cpu time limit from their default (the default is currently an hour)\n");
		fout.write("### The format is:\n");
		fout.write("### #PBS -l walltime=HH:MM:SS\n");
		fout.write("### #PBS -l cput=HH:MM:SS\n");
		fout.write("#PBS -l walltime=36:00:00\n"); #changed from 10h to 36h because of the very long zebrafish runs
		fout.write("#PBS -l cput=36:00:00\n"); #changed from 10h to 36h because of the very long zebrafish runs
		
		fout.write("### Move all your environment variables to the job\n");
		fout.write("#PBS -V\n");
		
		fout.write("### Set the working path to be used by the job\n");
		fout.write(Template("#PBS -d $OUTPUT_DIR\n\n").substitute(OUTPUT_DIR=sampleOutputFolder));
		
		
		#outdated way to redirect stdout and stderr - has security problems in copying the files after writing them
		#fout.write("### Optionally specify destinations for your program's output.\n");
		#fout.write("### Specify localhost and an NFS filesystem to prevent file copy errors.\n");
		#fout.write("### Without these lines or if the paths are invalid, PBS will redirect\n");
		#fout.write("### stdout and stderr to files in your home directory.\n");
		#fout.write(Template("#PBS -e localhost:$OUTPUT_DIR/queueLog.txt -e $OUTPUT_DIR/queueErr.txt\n").substitute(OUTPUT_DIR=sampleOutputFolder));
		#fout.write(Template("#PBS -o localhost:$OUTPUT_DIR/queueLog.txt -e $OUTPUT_DIR/queueErr.txt\n\n").substitute(OUTPUT_DIR=sampleOutputFolder));

		fout.write(Template("mkdir -p $OUTPUT_DIR\n").substitute(OUTPUT_DIR=sampleOutputFolder));
		fout.write(Template("exec 2> $OUTPUT_DIR/queueErr.txt > $OUTPUT_DIR/queueLog.txt\n\n").substitute(OUTPUT_DIR=sampleOutputFolder));
		
		
		fout.write("### Switch to the working directory; by default Torque launches processes from your home directory.\n");
		fout.write("### Jobs should only be run from /work; Torque returns results via NFS.\n");
		fout.write("echo Working directory is $PBS_O_WORKDIR\n");
		fout.write("cd $PBS_O_WORKDIR\n\n");
		 
		fout.write("### Run some informational commands.\n");
		fout.write("echo Running on host `hostname`\n");
		fout.write("echo Time is `date`\n");
		fout.write("echo Directory is `pwd`\n");
		fout.write("echo This jobs runs on the following processors:\n");
		fout.write("echo `cat $PBS_NODEFILE`\n\n");
 
		fout.write("### Define number of processors\n");
		fout.write("NPROCS=`wc -l < $PBS_NODEFILE`\n");
		fout.write("echo This job has allocated $NPROCS cpus\n\n");

		fout.write(cmd + '\n');


	print(Template("PBS script $RUN_FILE written\n").substitute(RUN_FILE=queueScriptFileName));

	st = os.stat(queueScriptFileName);
	os.chmod(queueScriptFileName, st.st_mode | stat.S_IEXEC);
	print("added run permissions to PBS script\n");

	if(not(args.do_not_send_to_cluster)):
		#actually execute the job by sending to cluster
		print("sending PBS script to queue\n");
		queueCmd = Template("qsub $RUN_FILE").substitute(RUN_FILE=queueScriptFileName);
		print(queueCmd);
		returnCode = subprocess.call(queueCmd, shell=True);
		if(returnCode != 0):
			raise Exception("invoking qsub failed")
	
	
		#queueCmd = Template("qsub -o $OUTPUT_DIR/queueLog.txt -e $OUTPUT_DIR/queueErr.txt \
	#	-d $OUTPUT_DIR -N $JOB_NAME \
	#	'$CMD'").substitute(CMD=cmd, OUTPUT_DIR=sampleOutputFolder, JOB_NAME=args.job_name);
	
	#print(queueCmd);
	#returnCode = subprocess.call(queueCmd, shell=True);
	#if(returnCode != 0):
	#	raise Exception("invoking qsub failed");
	