#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015

import os;
import glob;
import argparse;
from string import Template;
import shutil;
import re;
import subprocess;


def MakeSureOnlyOneCellInDir(dirname):
    for possibleLane in [1, 2]:
        globPattern = os.path.join(dirname, "*_L%03d_R?_???.fastq.gz" % possibleLane);
        print "globPattern is: " + globPattern;

        inFileList = glob.glob(globPattern);
        inFiles = " ".join(inFileList);
        print "inFiles are: " + inFiles;

        if(inFileList): #the list is not empty - there are files in this lane
            #make a new dir only for the current lane and move the relevant files into it
            newDirOnlyForCurrentLane = dirname + "_Lane%03d" % possibleLane
            if not os.path.exists(newDirOnlyForCurrentLane):
			    os.makedirs(newDirOnlyForCurrentLane);

            for file in inFileList:
                shutil.move(file, newDirOnlyForCurrentLane)

            #copy the SampleSheet.csv file into the daughter directory
            shutil.copy(os.path.join(dirname, "SampleSheet.csv"), newDirOnlyForCurrentLane)


    #After going through all the lanes, you can delete the original "SampleSheet.csv" file
    os.remove(os.path.join(dirname, "SampleSheet.csv"))

    for _, _, files in os.walk(dirname):
        #if we indeed covered all the possible lanes (and we deleted the SampleSheet.csv file) - the directory should now be empty
        if(files):
            raise Exception("The directory %s is not empty after the operation?!" % dirname)
    os.rmdir(dirname)


#side is either '1' or '2'
def OperateOnDirectoryOfOneCell(dirname, side):
    globPattern = os.path.join(dirname, "*_R" + side + "_???.fastq.gz");
    print "globPattern is: " + globPattern;
    
    inFileList = glob.glob(globPattern);
    inFiles = " ".join(inFileList);
    print "inFiles are: " + inFiles;
    
    
    if(inFileList): #the list is not empty - there are files to merge!
        
        outFile = os.path.join(dirname, \
                    re.sub("_R" + side + "_(\d\d\d).fastq.gz$", "_R" + side + "_combined.fastq.gz", inFileList[0])
                    );
        
        if(len(inFileList) == 1):
                    #only one file in subdir, no need to merge, just change its name by copying/moving
                    #depending on whether the user wanted the originals deleted or not
                    if(args.delete_originals):
                                        print "moving " + inFiles + " to " + outFile;
                                        shutil.move(inFiles, outFile);
                    else:
                                        print "copying " + inFiles + " to " + outFile;
                                        shutil.copyfile(inFiles, outFile);
                    
        else:
                    #note that OUTFILE already has a ".gz" extension - let us remove it before doing gzip, otherwise gzip does not work...
                    outFile, extension = os.path.splitext(outFile);
                    if(extension != ".gz"):
                                        raise Exception("unexpected extension: " + extension);

                    # IMPORTANT: this may be very expensive and unneeded --> according to the format, gzip files can simply be concatenated (even as they are compressed!...)
                    # HOWEVER I found some online discussions saying that FastQC may have trouble with that format because it's lame in its gzip implementation, so I decided to "play it safe" and concatenate gzips in the expensive way by gunzipping them first
                    #use "gzip -f" to force overwriting even if the file already exists
                    gunzipCommand = Template("gunzip -c $INFILES > $OUTFILE && gzip -f $OUTFILE").substitute(INFILES=inFiles, OUTFILE=outFile);
                    print("executing: " + gunzipCommand);
                    returnCode = subprocess.call(gunzipCommand, shell=True);
                    if(returnCode != 0):
                        raise Exception("gunzip failed");
                    
                    if(args.delete_originals):
                        print "deleting the original split files!"
                        for origFile in inFileList:
                                            print ".....deleting file: " + origFile;
                                            os.remove(origFile);




parser = argparse.ArgumentParser(description="Scan a folder of the BRAIN project and concatenate split fasta files into one file in preparation for preprocessing")
parser.add_argument("diretoryToProcess", action="store",
                    help="The directory with the original files to normalize");
parser.add_argument('-d', '--delete_originals', action="store_true",
                    help="delete the original split files that have just been merged");

args = parser.parse_args();

#In the batch of 3/16/2015 there were cells that were sequenced in two lanes
#--> separate them into two different directories per cell. This operation ALWAYS deletes originals (no matter what the
#command line flag says) to prevent the same cell from appearing twice in the data
MAKE_SURE_ONLY_ONE_CELL_IN_DIR = False
if(MAKE_SURE_ONLY_ONE_CELL_IN_DIR):
    for root, dirs, files in os.walk(args.diretoryToProcess):
        for d in dirs:
            fullDirPath = os.path.join(root, d);
            print "Operating on single cell directory: " + fullDirPath;
            MakeSureOnlyOneCellInDir(fullDirPath)



for root, dirs, files in os.walk(args.diretoryToProcess):
    for d in dirs:       
        fullDirPath = os.path.join(root, d);
        print "Operating on single cell directory: " + fullDirPath;
        OperateOnDirectoryOfOneCell(fullDirPath, '1');
        OperateOnDirectoryOfOneCell(fullDirPath, '2');
        
    