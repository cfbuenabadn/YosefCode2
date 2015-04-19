# 'gene40159'    'Prdx5' '1'  rna91078
#   'gene35896'    'Psmd2'           '1'   rna81403,rna81404
import os
import pysam
import numpy
import sys
import subprocess
from string import Template

#a hashable wrapper for the read, so I could use it in sets etc.
class AlignedRead:
    #init with an aligned segment
    def __init__(self, mySegment):
        self.mySegment = mySegment
        self.readName = mySegment.query_name
        self.phred20Score = []
        self.meanPhredScore = []

    def __hash__(self):
        return self.readName.__hash__()

    def __eq__(self, other):
        return self.readName == other.readName

    def GetPhred20Score(self):
        if(not(self.phred20Score)):
            quals = numpy.array(self.mySegment.query_qualities)
            self.phred20Score = numpy.mean(quals >= 20)
            self.meanPhredScore = numpy.mean(quals)

        return self.phred20Score

    def GetMeanPhredScore(self):
        if(not(self.phred20Score)):
            self.GetPhred20Score()

        return self.meanPhredScore

    # An object is hashable if it has a hash value which never changes during its lifetime (it needs a __hash__() method),
    #  and can be compared to other objects (it needs an __eq__() or __cmp__() method). Hashable objects which compare equal must have the same hash value.

transcriptID = 'rna81404'; 'rna81403'; 'rna91078';
seAlignmentFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined/rsem_output/aligned_by_bowtie2.bam'
peAlignmentFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai/OEP02_N712_S508_GTAGAGGA-CTAAGCCT_L002_R1_combined/rsem_output/aligned_by_bowtie2.bam'


seAlignmentFile = os.path.expanduser(seAlignmentFile)
peAlignmentFile = os.path.expanduser(peAlignmentFile)


# for infile in [seAlignmentFile, peAlignmentFile]:
    # print("Sort SAM")
    # sortSamCommand = Template("samtools sort -f -@ 8 $INPUT_FILE $INPUT_FILE.sorted.bam").substitute(INPUT_FILE=infile)
    # print(sortSamCommand)
    # sys.stdout.flush()
    # returnCode = subprocess.call(sortSamCommand, shell=True)
    # if(returnCode != 0):
    #     raise Exception("sortSam failed");

    # print("**********************************************************");
    # print("**********************************************************");
    # print("Index SAM");
    # indexSamCommand = Template("samtools index $INPUT_FILE").substitute(INPUT_FILE=infile+".sorted.bam");
    # print(indexSamCommand)
    # sys.stdout.flush();
    # returnCode = subprocess.call(indexSamCommand, shell=True);
    # if(returnCode != 0):
    #     raise Exception("indexSam failed");

seAlignmentFile += ".sorted.bam"
peAlignmentFile += ".sorted.bam"
peSamFile = pysam.AlignmentFile(peAlignmentFile)
seSamFile = pysam.AlignmentFile(seAlignmentFile)


#take reads that are primary mapped
readsMappedToTranscriptInSE = [AlignedRead(read) for read in seSamFile.fetch(transcriptID) if read.flag in (0, 16)]
#take the primary mapping of the first read. According to SAM's format documentation, exactly one read should have flag & 0x900 == 0 and this is the primary read
readsMappedToTranscriptInPE = [AlignedRead(read) for read in peSamFile.fetch(transcriptID) if (read.flag & 0x3 == 3) and (read.flag & 0xc == 0) and (read.flag & 0x900 == 0)]
#don't take and (read.flag & 0x40)  to insist on first segment because I found that the order may not be preserved and I want to be lenient in PE because this is wat we subtract


print len(readsMappedToTranscriptInSE)
print len(readsMappedToTranscriptInPE)

lostReads = set([read for read in readsMappedToTranscriptInSE]) - set([read for read in readsMappedToTranscriptInPE])
print len(lostReads)
remainingReads = set([read for read in readsMappedToTranscriptInSE]) & set([read for read in readsMappedToTranscriptInPE])
print len(remainingReads)

print "remaining"
for read in remainingReads:
    # print read.mySegment
    # print read.GetPhred20Score()
    print read.mySegment

print "lost"
for read in lostReads:
    # print read.mySegment
    # print read.GetPhred20Score()
    print read.mySegment

print numpy.median([read.GetPhred20Score() for read in remainingReads])
print numpy.median([read.GetPhred20Score() for read in lostReads])
print numpy.median([read.GetMeanPhredScore() for read in remainingReads])
print numpy.median([read.GetMeanPhredScore() for read in lostReads])

#take reads that are primary mapped
readsMappedToTranscriptInSE = [read for read in seSamFile.fetch(transcriptID) if read.flag in (0, 16)]
#take the primary mapping of the first read. According to SAM's format documentation, exactly one read should have flag & 0x900 == 0 and this is the primary read
readsMappedToTranscriptInPE = [read for read in peSamFile.fetch(transcriptID) if (read.flag & 0x3 == 3) and (read.flag & 0xc == 0) and (read.flag & 0x900 == 0)]
#don't take and (read.flag & 0x40)  to insist on first segment because I found that the order may not be preserved and I want to be lenient in PE because this is wat we subtract


print len(readsMappedToTranscriptInSE)
print len(readsMappedToTranscriptInPE)

lostReadNames = set([read.query_name for read in readsMappedToTranscriptInSE]) - set([read.query_name for read in readsMappedToTranscriptInPE])
print len(lostReadNames)
remainingReadsNames = set([read.query_name for read in readsMappedToTranscriptInSE]) & set([read.query_name for read in readsMappedToTranscriptInPE])
print len(remainingReadsNames)




# seSamFile.reset()
# print "remaining"
# for read in seSamFile:
#     if(read.query_name in remainingReadsNames):
#         print read
#         #print read.mapping_quality
#
# seSamFile.reset()
# print "lost"
# for read in seSamFile:
#     if(read.query_name in lostReadNames):
#         print read


        #print read.mapping_quality
# print len(lostReadNames)
# print (readsMappedToTranscriptInSE)
# print (readsMappedToTranscriptInPE)
#
# for read in readsMappedToTranscriptInPE:
#     print read
#
# for read in samfile.fetch(transcriptID):
#     # if(read.flag):
#         print hex(read.flag)

peSamFile.close()
seSamFile.close()