# Finds how many CG sites are targeted in each species for a selected set of Infinium 2 probes
import sys
import gzip
from collections import defaultdict

def main():
    if len(sys.argv) != 3:
        print ("Usage: python pickInf2Probes.py <probes file> <output file>")
        exit(1)
    probesFile = open(sys.argv[1], 'r')
    oFile = open(sys.argv[2], 'w')

    speciesCount = defaultdict(int)
    for line in probesFile:
        splitLine = line.strip().split("\t")
        species = splitLine[24].split(",")
        for curSpecies in species:
            if (curSpecies == "A"):
                print(line)
            speciesCount[curSpecies] += 1

    for species in sorted(speciesCount.keys()):
        if (species != "hg19"):
            oFile.write(species + " " + str(speciesCount[species]) + "\n")
    oFile.close()
    probesFile.close()

main()
