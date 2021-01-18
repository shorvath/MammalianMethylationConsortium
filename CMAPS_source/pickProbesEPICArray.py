# Picks a certain number of probes from the EPIC array ranked in order of number of species we can cover with that
# probe design.

from __future__ import print_function
import sys
import gzip
from collections import defaultdict
import operator

def main():
    if len(sys.argv) != 8:
        print("Usage: python pickProbesEPICArray.py <file of sites already picked> <all probes file> <EPIC manifest file> <num EPIC sites> <only Infinium II 0/1> <output file> <allow >3 CpGs 0/1>")
        exit(1)
    existingSites = gzip.open(sys.argv[1], 'rt')
    allScoresSites = gzip.open(sys.argv[2], 'rt')
    EPICfile = gzip.open(sys.argv[3], 'rt')
    numEPICsites = int(sys.argv[4])
    onlyInf2 = int(sys.argv[5])
    oFile = gzip.open(sys.argv[6], 'wt')
    allow3 = int(sys.argv[7])

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

    print("Parsing file for probes that aren't already included and are on epic array. . .")
    numSpeciesEPIC = defaultdict(int)
    probePicked = defaultdict(str)
    for line in allScoresSites:
        splitLine = line.strip().split("\t")

        if (allow3 or int(splitLine[19]) <= 3):
            chr = "chr" + splitLine[3]
            coord = int(splitLine[4])

            if (splitLine[15] == "F"):
                probeStrand = "+"
            else:
                probeStrand = "-"
            # if it's on converted strand since EPIC is all converted, not already picked and in the EPIC array
            #if (chr, coord) in EPICdesign:
            #    print ("converted: ", splitLine[8] == "C")
            #    print ("already picked: ", (chr, coord) not in pickedSites)
            #    print ("On EPIC:", (chr, coord) in EPICdesign)
            #    print ("Same strand as EPIC", EPICdesign[(chr, coord)][1] == probeStrand)
            if (splitLine[17] == "C") and ((chr, coord) not in pickedSites) and ((chr, coord) in EPICdesign) and (EPICdesign[(chr, coord)][1] == probeStrand):
                if (EPICdesign[(chr, coord)][0] == "I"):
                    if (not onlyInf2):
                        if (len(splitLine) > 23):
                            species = splitLine[23].split(",")
                            numSpeciesEPIC[(chr, coord)] = len(species)
                            SNVlocation = splitLine[27]
                            probePicked[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[23:30]) + "\tInf1\t1\t1\n"
                elif (EPICdesign[(chr, coord)][0] == "II"):
                    if len(splitLine) > 30:  # if we have an infinium 2 for this probe given design score and underlying CG count etc
                        species = splitLine[30].split(",")
                        numSpeciesEPIC[(chr, coord)] = len(species)
                        SNVlocation = splitLine[34]
                        probePicked[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[30:]) + "\tInf2\t1\t1\n"
                else:
                    print("Rogue non I or II Infinium design found.")
    print ("Done.")

    print("Sorting EPIC sites that weren't already picked by how many species they cover. . .")
    sorted_CGsites_EPIC = sorted(numSpeciesEPIC.items(), key = operator.itemgetter(1))
    sorted_CGsites_EPIC.reverse()
    print("Done.")

    print ("Picking EPIC sites. . .", )
    numProbesPicked = 0
    i = 0
    while (numProbesPicked < numEPICsites):
        CGsite = sorted_CGsites_EPIC[i][0]
        if EPICdesign[CGsite][0] == "II":
            numProbesPicked += 1
        elif EPICdesign[CGsite][0] == "I":
            numProbesPicked += 2
        else:
            print("Rogue non I or II Infinium design found")

        if (numProbesPicked <= numEPICsites):
            oFile.write(probePicked[CGsite])
        i += 1
    print ("Done.")

    oFile.close()
    existingSites.close()
    allScoresSites.close()
    EPICfile.close()

main()