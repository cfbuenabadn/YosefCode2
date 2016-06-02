#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# May 2015


#NOTE that the use of os.walk here is not that crisp because this is a recursive function that iterates on all the subdirectories tree...


import os;
import glob;
import argparse;
from string import Template;
import shutil;
import re;
import subprocess;
import sys


def MakeSureDataInExpectedFormat(dirname):
    #returns the names of the two files in the directory (the two files are side1 and side2 of the paired-end)

    #make sure there are exactly one bz2 file per lane
    globPattern = os.path.join(dirname, "*.bz2");
    print "globPattern is: " + globPattern;
    inFileList = glob.glob(globPattern);

    if(len(inFileList) != 2):
        raise Exception("Expected exactly two bz2 files in the folder %s but found %d".format(dirname, len(inFileList)))

    #make sure one is lane 1 and one is lane 2. For some cells, the L is not always L001 but also L002
    globPatternFile11 = os.path.join(dirname, "*_L001_R1.fastq.bz2");
    print "globPatternFile1 is: " + globPatternFile11;
    inFileList1 = glob.glob(globPatternFile11);

    if(len(inFileList1) == 1):
        #Found L001_R2 file, the second file should be L001_R2
        globPatternFile12 = os.path.join(dirname, "*_L001_R2.fastq.bz2");
    elif(len(inFileList1) > 1):
        raise Exception("should never reach this point (because we already checked there were two bz2 files in the directory)")
    else:
        #no files matched the L001 pattern, make sure that the files are L002 rather than L001
        print "pattern with L001 did not match, looking for pattern with L002"
        globPatternFile11 = os.path.join(dirname, "*_L002_R1.fastq.bz2");
        print "globPatternFile1 is: " + globPatternFile11;
        inFileList1 = glob.glob(globPatternFile11);

        if(len(inFileList1) != 1):
            raise Exception("Weird! The sample is not L001 nor L002")

        globPatternFile12 = os.path.join(dirname, "*_L002_R2.fastq.bz2");
        pattern = "L002"

    print "globPatternFile2 is: " + globPatternFile12;
    inFileList2 = glob.glob(globPatternFile12);
    if(len(inFileList2) != 1):
        raise Exception("unexpected bz2 files in the directory")

    side1_filename = inFileList1[0]
    side2_filename = inFileList2[0]
    return side1_filename, side2_filename





parser = argparse.ArgumentParser(description="Scan a folder of the BRAIN project and concatenate split fasta files into one file in preparation for preprocessing")
parser.add_argument("diretoryToProcess", action="store",
                    help="The directory with the original files to normalize");
parser.add_argument('-d', '--delete_originals', action="store_true",
                    help="delete the original split files that have just been merged");

args = parser.parse_args();

#if the path begins with a tilde - expand it to the user's homedir
args.diretoryToProcess = os.path.expanduser(args.diretoryToProcess);


DATA_ALREADY_SPLIT_TO_DIRECTORIES = False
if(DATA_ALREADY_SPLIT_TO_DIRECTORIES):

    #get only top level directories in pyton:
    #[ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
    #os.walk('.').next()[1]
    #See: http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python

    dirs = os.walk(args.diretoryToProcess).next()[1] #os.walk returns a generator for 3-tuples, with the 2nd element being the directory, that's why I do next() and then take [1]
    dirs = [d for d in dirs if (d.startswith("Sample_LIB") or d.startswith("LIB"))]  #take only into the libraries that contain single cells
    filesToWorkOn = []
    for d in dirs:
        fullDirPath = os.path.join(args.diretoryToProcess, d);
        print "Operating on single cell directory: " + fullDirPath;
        side1_filename, side2_filename = MakeSureDataInExpectedFormat(fullDirPath)
        filesToWorkOn.append(side1_filename)
        filesToWorkOn.append(side2_filename)

    print "Successfully verified the directory's expected format"

else:
    #this code assumes only _001 namely no multiple read files per cell

    filesToWorkOn = []
    filenames = set(os.walk(args.diretoryToProcess).next()[2])
    pattern = re.compile(r"_R[12](_001)?\.fastq\.(indexCount\.tab|bz2|gz)$")
    samplePrefixes =  set([re.sub(pattern, "", f) for f in filenames])
    #make sure the "undetermined reads" are not interpreted as a sample
    samplePrefixes = [x for x in samplePrefixes if not re.search("undetermined", x, re.IGNORECASE)]

    for prefix in samplePrefixes:
        #the  bz2 files:
        #relevantFiles = [prefix + s for s in ["_R1.fastq.bz2", "_R2.fastq.bz2", "_R1.fastq.indexCount.tab"]]
        relevantFiles = [prefix + s for s in ["_R1_001.fastq.gz", "_R2_001.fastq.gz"]]
        if(not(set(relevantFiles).issubset(filenames))):
            #I decided to turn off that check when I added gz support (wheas previously I had a different format with bz2)
            raise Exception("unexpected files in directory - not all required files are present")


        newDirName = os.path.join(args.diretoryToProcess, prefix)
        os.makedirs(newDirName)
        for f in relevantFiles:
            shutil.move(os.path.join(args.diretoryToProcess, f), newDirName)
            newFileName = os.path.join(newDirName, f)

            if(newFileName.endswith("bz2")):
                #add only the bz2 files to the ones that should be transferred to gzip
                filesToWorkOn.append(newFileName)

            if(not(f in filenames)):
                raise Exception("should never happen - expected file name is not in the set of files")

            filenames.remove(f)

    #filenames is not empty at this point!
    if(not(bool(filenames))):
        raise Exception("unexpected extra files in directory")

    print "Successfully verified the directory's expected format"

#This is the original else clause that supported only bz2 from before I got the gz files (15 March 2016)
# else:
#
#     filesToWorkOn = []
#     filenames = set(os.walk(args.diretoryToProcess).next()[2])
#     pattern = re.compile(r"_R[12]\.fastq\.(indexCount\.tab|bz2)$")
#     samplePrefixes =  set([re.sub(pattern, "", f) for f in filenames])
#
#     for prefix in samplePrefixes:
#         relevantFiles = [prefix + s for s in ["_R1.fastq.bz2", "_R2.fastq.bz2", "_R1.fastq.indexCount.tab"]]
#         if(not(set(relevantFiles).issubset(filenames))):
#             raise Exception("unexpected files in directory - not all required files are present")
#
#
#         newDirName = os.path.join(args.diretoryToProcess, prefix)
#         os.makedirs(newDirName)
#         for f in relevantFiles:
#             shutil.move(os.path.join(args.diretoryToProcess, f), newDirName)
#             newFileName = os.path.join(newDirName, f)
#
#             if(newFileName.endswith("bz2")):
#                 #add only the bz2 files to the ones that should be transferred to gzip
#                 filesToWorkOn.append(newFileName)
#
#             if(not(f in filenames)):
#                 raise Exception("should never happen - expected file name is not in the set of files")
#
#             filenames.remove(f)
#
#     #filenames is not empty at this point!
#     if(not(bool(filenames))):
#         raise Exception("unexpected extra files in directory")
#
#     print "Successfully verified the directory's expected format"







#now, transform the bz2 files into gz files with the expected names:
#(if there were no bz2 files, the list filesToWorkOn will be empty)
pattern = re.compile("_R(?P<side>[12]).fastq.bz2$")
for fname in filesToWorkOn:
    #\g: escaping to refer to a named group
    outFile = pattern.sub("_R\g<side>_combined.fastq.gz", fname)
    print "Infile is: %s\nOutfile is: %s"  % (fname, outFile)

    print "decompressing and recompressing..."
    cmd = Template("bunzip2 -c < $INFILE | gzip -c > $OUTFILE").substitute(INFILE=fname, OUTFILE=outFile)
    print("executing: " + cmd);
    returnCode = subprocess.call(cmd, shell=True);
    if(returnCode != 0):
        raise Exception("bz2 or gzip failed");

    if(args.delete_originals):
        print "deleting the original file %s!" % fname
        os.remove(fname);


print "normalization of the directory %s complete" %args.diretoryToProcess

