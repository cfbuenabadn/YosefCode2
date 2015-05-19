# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# May 2015

import argparse
import os
import shutil
import gzip

FILES_TO_REMOVE = [os.path.join('fastqc_output', '1.Ptrim_fastqc'),
                   os.path.join('fastqc_output', '1.Utrim_fastqc'),
                   os.path.join('fastqc_output', '2.Ptrim_fastqc'),
                   os.path.join('fastqc_output', '2.Utrim_fastqc'),
                   os.path.join('rsem_output', 'aligned_by_bowtie2.bam'),
                   os.path.join('rsem_output', 'accepted_hits_noMultiple.bam'),
                   os.path.join('rsem_output', 'accepted_hits.bam'),
                   os.path.join('rsem_output', 'unmapped.bam'),
                   os.path.join('tophat_output', 'accepted_hits.bam'),
                   os.path.join('tophat_output', 'accepted_hits_noMultiple.bam')
                   ]

FILES_TO_COMPRESS = [os.path.join('trimmomatic_output', 'trimmomatic_log.txt')]


def cleanTemporaryFiles(folder):
    print("**********************************************************")
    print("**********************************************************")
    print("Cleaning temporary files in folder {0}".format(folder))

    filesToRemove = [os.path.join(folder, fname) for fname in FILES_TO_REMOVE]
    for fname in filesToRemove:
        if os.path.exists(fname):
            if(os.path.isdir(fname)):
                shutil.rmtree(fname)
            elif(os.path.isfile(fname)):
                os.remove(fname)
            else:
                raise Exception("cleanup: should not reach this point");


    print("**********************************************************")
    print("**********************************************************")
    print("Compressing some files in folder {0}".format(folder))

    filesToCompress = [os.path.join(folder, fname) for fname in FILES_TO_COMPRESS]
    for fname in filesToCompress:
        if os.path.isfile(fname):
            f_in = open(fname, 'rb')
            f_out = gzip.open(fname+'.gz', 'wb')
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()

            os.remove(fname)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process one single cell sample")
    parser.add_argument('folder', action='store',
                   help='the folder to clean')
    parser.add_argument('-r', '--recursive', action='store_true',
                   help="if unset, then the folder had a single sample; if set, then the folder contains multiple subfolders and is scanned recursively")


    args = parser.parse_args()

    #if the path begins with a tilde - expand it to the user's homedir
    args.folder = os.path.expanduser(args.folder)

    if(not(args.folder_has_multiple_samples)):
        #this is a single sample
        cleanTemporaryFiles(args.folder)
    else:
        #clear a folder with multiple samples

        #note that the recursive os.walk will descend into the existing subfolders of the single samples folders, but we don't care because it will catch nothing there
        for root, dirs, files in os.walk(args.folder):
            for d in dirs:
                fullDirPath = os.path.join(root, d)
                print "Operating on directory: " + fullDirPath
                cleanTemporaryFiles(fullDirPath)

    print "done cleaning {0}".format(args.folder)