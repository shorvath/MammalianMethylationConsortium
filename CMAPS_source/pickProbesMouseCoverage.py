# Picks all Infinium 2 probes that cover mouse when picking the maximum number of species for each probe. Selects
# all probes that cover mouse, then up to the desired number of probes ranked in order of conservation.
from __future__ import print_function
import sys
import gzip
from collections import defaultdict
import operator

def main():
    if len(sys.argv) != 7:
        print("Usage: python pickProbesMouseCoverage.py <probes file> <EPIC file> <output file> <number of sites mouse> <total number of sites> <converted only 0/1>")
        exit(1)
    probesFile = gzip.open(sys.argv[1], 'rt')
    EPICfile = gzip.open(sys.argv[2], 'rt')
    oFile = gzip.open(sys.argv[3], 'wt')
    numSitesMouse = int(sys.argv[4])
    totalSites = int(sys.argv[5])
    convertedOnly = int(sys.argv[6])

    oFile.write("Seq_ID\tForward_Sequence\tGenome_Build\tChromosome\tCoordinate\tDesign_State\tSeq_Length\tForward_CpG_Coord\tTB_Strand\tTop_Sequence\tTop_CpG_Coord\tProbe_Type\tProbeset_ID\tProbeset_Score\tMethyl_Probe_Sequence\tAllele_FR_Strand\tAllele_TB_Strand\tAllele_CO_Strand\tMethyl_Probe_Score\tUnderlying_CpG_Count\tUnMethyl_Probe_Sequence\tUnmethyl_Probe_Score\tInfinium_ScoreSeq_ID\tNum_Species\tSpecies\tProbe_Start_Coord\tProbe_End_Coord\tReference_Probe_Sequence\tSNV_location\tSNV_original\tSNV_change\tInfinium_Type\tIs_EPIC_site\tIs_EPIC_design\n")

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

    # Extract Infinium 2 probe numbers
    print("Parsing file for Infinium 2 probes. . .")
    CGsites_Inf2 = defaultdict(int)
    mouseProbes = defaultdict(int)
    probePicked_Inf2 = defaultdict(str)
    whichSpecies_Inf2 = defaultdict(list)

    for line in probesFile:
        splitLine = line.strip().split("\t")

        if (int(splitLine[19]) <= 3):
            # check if we are picking opposite strand or just converted
            if (convertedOnly and splitLine[17] != "C"):
                continue

            if splitLine[15] == "F":
                strand = "+"
            else:
                strand = "-"

            if len(splitLine) > 30: # if we have an infinium 2 for this probe given design score and underlying CG count etc
                species = splitLine[30].split(",")
                if len(species) > 1:
                    chr = "chr" + splitLine[3]
                    coord = int(splitLine[4])
                    SNVlocation = splitLine[34]
                    SNVoriginal = splitLine[35]
                    SNVchoice = splitLine[36]

                    if len(species) - 1 > CGsites_Inf2[(chr, coord)]:
                        CGsites_Inf2[(chr, coord)] = len(species) - 1
                        whichSpecies_Inf2[(chr, coord)] = species

                        if "mm10" in species:
                            mouseProbes[(chr, coord)] = len(species) - 1
                        if ((chr, coord) in EPICdesign): # if it's on the EPIC array
                            if (EPICdesign[(chr, coord)][0] == "II") and (splitLine[17] == "C") and (strand == EPICdesign[(chr, coord)][1]): # if it's the same design
                                probePicked_Inf2[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[30:]) + "\tInf2\t1\t1\n"
                            else:
                                probePicked_Inf2[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[30:]) + "\tInf2\t1\t0\n"
                        else:
                            probePicked_Inf2[(chr, coord)] = "\t".join(splitLine[:23]) + "\t" + str(len(species)) + "\t" + "\t".join(splitLine[30:]) + "\tInf2\t0\t0\n"
    print("Done.")

    print("Sorting sites that cover mouse by how many other species they cover. . .", end = "")
    sorted_mouseProbes = sorted(mouseProbes.items(), key = operator.itemgetter(1))
    sorted_mouseProbes.reverse()
    print("Done.")

    print("Picking Infinium 2 sites that cover mouse. . .")
    pickedSites = defaultdict(bool)
    numPicked = 0
    i = 0
    while (numPicked < numSitesMouse and i < len(sorted_mouseProbes)):
        CGsite = sorted_mouseProbes[i][0]
        pickedSites[CGsite] = True
        oFile.write(probePicked_Inf2[CGsite])
        numPicked += 1
        i += 1
    print("Done. Picked ", numPicked, " sites.")

    print("Sorting all Infinium 2 sites based on how many species they cover. . .", end = "")
    sorted_CGsites_Inf2 = sorted(CGsites_Inf2.items(), key = operator.itemgetter(1))
    sorted_CGsites_Inf2.reverse()
    print("Done.")

    print("Picking Infinium 2 sites. . .")
    i = 0
    while (numPicked < totalSites and i < len(sorted_CGsites_Inf2)):
        CGsite = sorted_CGsites_Inf2[i][0]
        if ((CGsite not in pickedSites)):
            pickedSites[CGsite] = True
            oFile.write(probePicked_Inf2[CGsite])
            numPicked += 1
        i += 1
    print("Done. Picked ", numPicked, " sites.")

    oFile.close()
    probesFile.close()

main()
