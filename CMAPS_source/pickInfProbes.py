# For each CG site this script picks the Infinium 2 probe with the most species covered and outputs how many
# species that probe covers. The output from this can be used to sort by number of species covered and look at the
# coverage of each species.
#
# This script assumes that each CG site has 4 probes in the input file (with Inf 1 and 2 results on the same line)
import sys
import gzip
from collections import defaultdict

def main():
    if len(sys.argv) != 4:
        print "Usage: python pickInfProbes.py <probes file> <output file> <infinium type 1/2>"
        exit(1)
    probesFile = gzip.open(sys.argv[1], 'r')
    oFile = gzip.open(sys.argv[2], 'w')
    infiniumType = int(sys.argv[3])

    CGSpeciesCount = defaultdict(int)
    probesPicked = defaultdict(str)
    for line in probesFile:
        splitLine = line.strip().split("\t")
        if (infiniumType == 2):
            if len(splitLine) > 14:
                species = splitLine[12].split(",")
                if len(species) > 1:
                    chr = "chr" + splitLine[2]
                    coord = int(splitLine[3])
                    SNVlocation = splitLine[13]
                    SNVchoice = splitLine[14]

                    if len(species) - 1 > CGSpeciesCount[(chr, coord)]:
                        CGSpeciesCount[(chr, coord)] = len(species) - 1
                        probesPicked[(chr, coord)] = chr + "\t" + str(coord) + "\t" + str(coord + 2) + "\t" + str(len(species) - 1) + "\t" + ",".join(species) + "\t" + SNVlocation + "\t" + SNVchoice + "\t" + splitLine[4] + "\t" + splitLine[5] + "\t" + splitLine[6] + "\t" + splitLine[7] + "\n"
        else:
            if len(splitLine) > 11:
                species = splitLine[9].split(",")
                if len(species) > 1:
                    chr = "chr" + splitLine[2]
                    coord = int(splitLine[3])
                    SNVlocation = splitLine[10]
                    SNVchoice = splitLine[11]

                    if len(species) - 1 > CGSpeciesCount[(chr, coord)]:
                        CGSpeciesCount[(chr, coord)] = len(species) - 1
                        probesPicked[(chr, coord)] = chr + "\t" + str(coord) + "\t" + str(coord + 2) + "\t" + str(len(species) - 1) + "\t" + ",".join(species) + "\t" + SNVlocation + "\t" + SNVchoice + "\t" + splitLine[4] + "\t" + splitLine[5] + "\t" + splitLine[6] + "\t" + splitLine[7] + "\n"

    for site in probesPicked:
        oFile.write(probesPicked[site])
    oFile.close()
    probesFile.close()

main()
