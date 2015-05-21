# Allon Wagner
# Nir Yosef's lab, UC Berkeley
# Mar 2015

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