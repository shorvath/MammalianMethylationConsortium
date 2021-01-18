# Picks 2K Infinium 1 probes

from __future__ import print_function
import sys
import gzip
from collections import defaultdict
import operator

def main():
    if len(sys.argv) != 7:
        print("Usage: python pickProbesInfinium1.py <file of sites already picked> <all probes file> <EPIC manifest file> <num Infinium 1 sites> <output file> <only converted 1/0>")
        exit(1)
    existingSites = gzip.open(sys.argv[1], 'rt')
    allScoresSites = gzip.open(sys.argv[2], 'rt')
    EPICfile = gzip.open(sys.argv[3], 'rt')
    numInfinium1sites = int(sys.argv[4])
    oFile = gzip.open(sys.argv[5], 'wt')
    convertedOnly = int(sys.argv[6])

    print("Parsing file for EPIC probe design . . .")
    EPICdesign = defaultdict(str)
    firstLine = True
    for line in EPICfile:
        splitLine = line.strip().split(",")
        if not firstLine:
            chr = splitLine[1][1:-1]
            coord = int(splitLine[2])
            strand = splitLine[3][1:-1]
            infType = splitLine[9][1:-1]

            EPICdesign[(chr, coord)] = (infType, strand)
        else:
            firstLine = False
    print("Done.")

    pickedSites = defaultdict(bool)
    print("Parsing file for sites already picked . . .")
    firstLine = True
    for line in existingSites:
        if not firstLine:
            splitLine = line.strip().split("\t")
            chr = "chr" + splitLine[3]
            coord = int(splitLine[4])

            pickedSites[(chr, coord)] = line
        else:
            firstLine = False
    print("Done.")

    # Extract Infinium 1 probe numbers
    print("Parsing file for Infinium 1 probes. . .")
    CGsites_Inf1 = defaultdict(int)
    probePicked_Inf1 = defaultdict(str)
    EPICconflict = defaultdict(bool)

    for line in allScoresSites:
        splitLine = line.strip().split("\t")

        if (convertedOnly and splitLine[17] != "C"):
            continue

        if splitLine[15] == "F":
            strand = "+"
        else:
            strand = "-"

        if len(splitLine) > 23: # if we have an infinium 1 for this probe given design score and underlying CG count etc
            species = splitLine[23].split(",")
            if len(species) > 1:
                chr = "chr" + splitLine[3]
                coord = int(splitLine[4])
                SNVlocation = splitLine[27]
                SNVoriginal = splitLine[28]
                SNVchoice = splitLine[29]

                if ((chr, coord) not in pickedSites):
                    if len(species) - 1 > CGsites_Inf1[(chr, coord)]:  # just pick Infinium 1 probe that covers most species
                        CGsites_Inf1[(chr, coord)] = len(species) - 1
                        if ((chr, coord) in EPICdesign): # if this site is on the epic array
                            if (EPICdesign[(chr, coord)][0] == "I") and (splitLine[17] == "C") and (strand == EPICdesign[(chr, coord)][1]): # if we have the same design
                                probePicked_Inf1[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[23:30]) + "\tInf1\t1\t1\n"
                            else:
                                probePicked_Inf1[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[23:30]) + "\tInf1\t1\t0\n"
                        else:
                            probePicked_Inf1[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[23:30]) + "\tInf1\t0\t0\n"
    print("Done.")

    print("Sorting all Infinium 1 sites based on how many species they cover. . .", end = "")
    sorted_CGsites_Inf1 = sorted(CGsites_Inf1.items(), key = operator.itemgetter(1))
    sorted_CGsites_Inf1.reverse()
    print("Done.")

    print("Picking Infinium 1 sites. . .")
    i = 0
    numPicked = 0
    while (numPicked < numInfinium1sites and i < len(sorted_CGsites_Inf1)):
        CGsite = sorted_CGsites_Inf1[i][0]
        if ((CGsite not in pickedSites)):
                pickedSites[CGsite] = True
                oFile.write(probePicked_Inf1[CGsite])
                numPicked += 1
        else:
            print("This shouldn't be happening. We've already picked this site with Infinium 2 but somehow didn't properly check for it.")
        i += 1
    print("Done.")

    oFile.close()
    existingSites.close()
    allScoresSites.close()
    EPICfile.close()

main()