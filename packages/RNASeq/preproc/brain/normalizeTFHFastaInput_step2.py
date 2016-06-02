#!/usr/bin/python

# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Jan 2015


#NOTE that the use of os.walk here is not that crisp because this is a recursive function that iterates on all the subdirectories tree...
#Or maybe it's not a bug but rather a feature... can work recursively on a folder...

#after the first step of normalize TFH, if the original files were gz they were not converted to combined_ in their names...

import os;
import glob;
import argparse;
from string import Template;
import shutil;
import re;
import subprocess;



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




parser = argparse.ArgumentParser(description="TFH normalize step 2")
parser.add_argument("diretoryToProcess", action="store",
                    help="The directory with the original files to normalize");
parser.add_argument('-d', '--delete_originals', action="store_true",
                    help="delete the original split files that have just been merged");

args = parser.parse_args();





for root, dirs, files in os.walk(args.diretoryToProcess):
    for d in dirs:       
        fullDirPath = os.path.join(root, d);
        print "Operating on single cell directory: " + fullDirPath;
        OperateOnDirectoryOfOneCell(fullDirPath, '1');
        OperateOnDirectoryOfOneCell(fullDirPath, '2');
        
    