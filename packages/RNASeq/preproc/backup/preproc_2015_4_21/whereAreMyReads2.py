
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


# seAlignmentFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/accepted_hits.bam.sortedByName.bam'
# peAlignmentFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/accepted_hits.bam.sortedByName.bam'

def ReadReadFileToDictionary(readTableFileName):
    table = {}
    count = 0
    with open(readTableFileName) as fin:
        for line in fin:
            count += 1
            parts = line.split()
            if parts[1].startswith('@'):
                continue

            table[parts[1]] = int(parts[0])


        if count % 10000:
            print '.'

    print '\n'
    return table

seReadTableFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/reads.txt'
peReadTableFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/reads.txt'
peUnmappedFile = '~/data/BRAIN/processed3/150202_HS2A/Project_Ngai/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/unmapped.sorted.bam'
seAlignmentFile = os.path.expanduser(seReadTableFile)
peAlignmentFile = os.path.expanduser(peReadTableFile)
peUnmappedFile = os.path.expanduser(peUnmappedFile)

seTable = ReadReadFileToDictionary(seAlignmentFile)
peTable = ReadReadFileToDictionary(peAlignmentFile)

lostReads = set(seTable.keys()) - set(peTable.keys())
print "SE reads: %d" % len(seTable.keys())
print "PE reads: %d" % len(peTable.keys())
print "lost reads: %d" % len(lostReads)

peUnmapped = pysam.AlignmentFile(peUnmappedFile)

peUnmapped1CorrespondingToLost = []
peUnmapped2CorrespondingToLost = []
peUnmappedNotCorrespondingToLost = []
for read in peUnmapped:
    curRead = AlignedRead(read)
    if (read.query_name in lostReads):
            if read.is_read1:
                peUnmapped1CorrespondingToLost.append(curRead)
            else:
                peUnmapped2CorrespondingToLost.append(curRead)


            # print read.query_name
            print read
    else:
        peUnmappedNotCorrespondingToLost.append(curRead)

peUnmapped.close()


[read.GetPhred20Score() for read in peUnmapped1CorrespondingToLost]
[read.GetPhred20Score() for read in peUnmapped2CorrespondingToLost]
[read.GetMeanPhredScore() for read in peUnmapped1CorrespondingToLost]
[read.GetMeanPhredScore() for read in peUnmapped2CorrespondingToLost]


unaccountedLost = lostReads - set(read.mySegment.query_name for read in peUnmapped1CorrespondingToLost)

seSamFile = os.path.expanduser('~/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle/OEP02_N712_S507_GTAGAGGA-AAGGAGTA_L002_R1_combined/tophat_output/accepted_hits.bam.sortedByName.bam')
seSam = pysam.AlignmentFile(seSamFile)
seLostAlignments = [AlignedRead(read) for read in seSam if (read.query_name in lostReads)]
seSam.reset()
seNonLostAlignments = [AlignedRead(read) for read in seSam if not(read.query_name in lostReads)]
seSam.close()

print numpy.median([read.GetPhred20Score() for read in seLostAlignments])
print numpy.median([read.GetPhred20Score() for read in seNonLostAlignments])
print numpy.median([read.GetMeanPhredScore() for read in seLostAlignments])
print numpy.median([read.GetMeanPhredScore() for read in seNonLostAlignments])

#
# for lostReadName in lostReads:
#     print lostReadName + ' ' + str(seTable[lostReadName])


# for readName, Appearances in d.iteritems():

# for infile in [seAlignmentFile, peAlignmentFile]:
#     print("Sort SAM")
#     sortSamCommand = Template("samtools sort -nf -@ 12 $INPUT_FILE $INPUT_FILE.sortedByName.bam").substitute(INPUT_FILE=infile)
#     print(sortSamCommand)
#     sys.stdout.flush()
#     returnCode = subprocess.call(sortSamCommand, shell=True)
#     if(returnCode != 0):
#         raise Exception("sortSam failed");
#
#     print("**********************************************************");
#     print("**********************************************************");
#     print("Index SAM");
#     indexSamCommand = Template("samtools index $INPUT_FILE").substitute(INPUT_FILE=infile+".sorted.bam");
#     print(indexSamCommand)
#     sys.stdout.flush();
#     returnCode = subprocess.call(indexSamCommand, shell=True);
#     if(returnCode != 0):
#         raise Exception("indexSam failed");


# def OutputOnlyUniquelyAlignedReads(inFile, outFile, isPaired):
#     with open(inFile) as fin, open(outFile, 'wt') as out:
#         line = fin.readline()
#         while(line):
#             readName = line.split()[0]
#             linesOfRead = []
#
#             curReadName = readName
#             while(curReadName == readName):
#                 line = fin.readline()
#                 curReadName = line.split()[0]
#
#
#
# seAlignmentFile = os.path.expanduser(seAlignmentFile)
# peAlignmentFile = os.path.expanduser(peAlignmentFile)
# peSamFile = pysam.AlignmentFile(peAlignmentFile)
# seSamFile = pysam.AlignmentFile(seAlignmentFile)
#
# print seSamFile
# seLine = seSamFile.next()
# peLine = peSamFile.next()

# with open(seAlignmentFile) as seInput, open(peAlignmentFile) as peInput:
#     seLine = seInput.next()
#     peLine = peInput.next()