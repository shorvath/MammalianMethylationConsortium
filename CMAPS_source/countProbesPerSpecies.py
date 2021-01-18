# Finds how many CG sites are targeted in each species for a selected set of Infinium 2 probes
import sys
import gzip
from collections import defaultdict

def main():
    if len(sys.argv) != 3:
        print "Usage: python pickInf2Probes.py <probes file> <output file>"
        exit(1)
    probesFile = gzip.open(sys.argv[1], 'r')
    oFile = gzip.open(sys.argv[2], 'w')

    speciesCount = defaultdict(int)
    for line in probesFile:
        splitLine = line.strip().split("\t")
        species = splitLine[4].split(",")
        for curSpecies in species:
            speciesCount[curSpecies] += 1

    for species in speciesCount:
        oFile.write(species + " " + str(speciesCount[species]) + "\n")
    oFile.close()
    probesFile.close()

main()
