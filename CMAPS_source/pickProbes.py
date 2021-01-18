# Selects probes based on a minimum number of species they have to cover
import sys
import gzip
from collections import defaultdict

def main():
    if len(sys.argv) != 6:
        print "Usage: python pickProbes.py <probes file> <infinium type> <number of species threshold> <output file> <probes/sites>"
        exit(0)
    probesFile = gzip.open(sys.argv[1], 'r')
    infiniumType = int(sys.argv[2])
    speciesThresh = int(sys.argv[3])
    oFile = gzip.open(sys.argv[4], 'w')
    probesOrSites = sys.argv[5]

    alreadyPicked = defaultdict(bool)
    for line in probesFile:
        splitLine = line.strip().split('\t')
        if infiniumType == 1:
            if len(splitLine) < 10:
                continue
            species = splitLine[9].split(",")
        else:
            if len(splitLine) < 13:
                continue
            species = splitLine[12].split(",")

        if len(species) >= speciesThresh:
            if probesOrSites == "probes":
                oFile.write(line)
            else:
                coord = int(splitLine[3])
                if not alreadyPicked[(splitLine[2], coord)]:
                    oFile.write("chr" + splitLine[2] + "\t" + str(coord) + "\t" + str(coord + 2) + "\n")
                    alreadyPicked[(splitLine[2], coord)] = True

    probesFile.close()
    oFile.close()

main()
