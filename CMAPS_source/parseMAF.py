import gzip
import sys
from collections import defaultdict
from Bio.Seq import Seq

def main():
    if len(sys.argv) < 7:
        print("Usage: python createMAFFeatures.py <alignment file> <species names file> <output file> <chromosome> <reference species> <reference species chromosome length file>")
        exit(1)
    mafFile = sys.argv[1]
    speciesFile = open(sys.argv[2], 'r')
    featFile = sys.argv[3]
    chr = sys.argv[4]
    refSpecies = sys.argv[5] + "." + chr
    chromLengths = sys.argv[6]

    f = gzip.open(mafFile, 'rb')
    o = gzip.open(featFile, 'w')
    chromSizes = open(chromLengths, 'r')

    species = []
    for line in speciesFile:
        species.append(line.strip())
    species = sorted(species)
    print species

    chrLength = {}
    for line in chromSizes:
        splitLine = line.strip().split()
        chrLength[splitLine[0]] = int(splitLine[1])

    idx = 0
    # write header to output file
    o.write("pos")
    for sp in species:
        o.write("," + sp)
    o.write("\n")

    f.readline() # read header line
    scoreLine = f.readline()
    while scoreLine[0] == '#':
        scoreLine = f.readline()

    numBases = 0
    covered = {} # dictionary of bases covered to account for blocks convering duplicate pieces
    while scoreLine != "":
        sequenceLine = f.readline()
        startHuman = -1
        seqHuman = ""
        alignedToHuman = {}

        lastAlignment = []
        reverseComplement = False
        # parse an alignment block
        while sequenceLine != "\n":
            lastAlignment.append(sequenceLine)
            splitSeq = sequenceLine.split()

            type = splitSeq[0]
            if type == 's':
                curSpecies = splitSeq[1]
                
                if curSpecies == refSpecies: # store info for human
                    trueLenHuman = int(splitSeq[3])
                    if splitSeq[4] == '-': # reverse complement
                        startHuman = chrLength[chr] - int(splitSeq[2]) - trueLenHuman
                        seqHuman = Seq(splitSeq[6].upper()).reverse_complement()
                        reverseComplement=True
                    else:
                        startHuman = int(splitSeq[2])
                        seqHuman = splitSeq[6].upper()
                    lenHumanDash = len(seqHuman)
                else:
                    curSeq = ""
                    if reverseComplement:
                        curSeq = Seq(splitSeq[6].upper()).reverse_complement()
                    else:
                        curSeq = splitSeq[6].upper()
                    alignedToHuman[curSpecies[:curSpecies.find(".")]] = curSeq # store the aligned sequence in other species
            sequenceLine = f.readline()

        # when done with an alignment block go through and output whether each species is the same as the reference species
        pos = startHuman
        for i in range(lenHumanDash):
            if seqHuman[i] != '-':
                numBases += 1
                if not (pos in covered):
                    o.write(str(pos) + "," + seqHuman[i])
                    for sp in species:
                        if (sp not in alignedToHuman) or (alignedToHuman[sp][i] == '-'):
                            o.write(",X")
                        else:
                            o.write("," + str(alignedToHuman[sp][i]))
                    o.write("\n")
                    covered[pos] = 1
                pos += 1

        scoreLine = f.readline()
    o.close()
main()


