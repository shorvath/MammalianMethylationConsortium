from __future__ import print_function
import sys
import gzip
from collections import defaultdict
from itertools import combinations

def getMaxSpecies(positionComb, maxSpeciesSubset, curChunk, startPos, endPos, refSpecies, pt = False):
    # Returns the maximal subset of species that can be covered with a probe starting at startPos and ending at endPos (end exclusive)
    # with SNVs at positionComb. Also returns the specific nucleotide choice we should design for at every SNV so that
    # the maximum number of species are covered.

    # First keep only the species where the CG is conserved. The CG is always at positions startPos + 49 and startPos + 50
    finalSNVchoice = []
    localMaxSpeciesSubset = set(maxSpeciesSubset)
    for pos in range(startPos, endPos):
        SNVspecies = defaultdict(lambda: defaultdict(list))  # for each possible SNV nucleotide, which species we can cover

        maxSpeciesSubsetCopy = set(localMaxSpeciesSubset)
        for species in maxSpeciesSubsetCopy:
            if (curChunk[pos][species] != curChunk[pos][refSpecies]):  # if we have an SNV in this species
                if (curChunk[pos][species] == "X"):
                    localMaxSpeciesSubset.remove(species) # can't cover if nothing in alignment here
                else:
                    if pos in positionComb:  # check if we have an SNV here
                        SNVspecies[pos][curChunk[pos][species]].append(species) # add to the list of species we can keep if we design for this nucleotide
                    else:
                        localMaxSpeciesSubset.remove(species) # can't cover if no SNP at this position and different
        # find the best SNV choice at this position if we need one
        if len(SNVspecies[pos]) > 0:
            SNVchoices = []
            for i in SNVspecies[pos]:
                SNVchoices.append((len(SNVspecies[pos][i]), i))
            SNVchoices = sorted(SNVchoices, reverse = True)

            count = positionComb.count(pos) # how many SNVs are we picking at this position
            SNVchoices = [x[1] for x in SNVchoices[:count]]
            SNVpurpose = defaultdict(list)
            # remove the species that don't have this SNV at this position
            for i in SNVspecies[pos]:
                if i not in SNVchoices:
                    for species in SNVspecies[pos][i]:
                        localMaxSpeciesSubset.remove(species)
                else:
                    for species in SNVspecies[pos][i]:
                        SNVpurpose[i].append(species)

            for i in SNVchoices:
                finalSNVchoice.append((pos, i, curChunk[pos][refSpecies], SNVpurpose[i]))
        #if pt: print(maxSpeciesSubset)
    return localMaxSpeciesSubset, finalSNVchoice

def getFeasibleSpecies(curChunk, startPos, endPos, refSpecies, numSNV):
    # Return the feasible species for a 100bp chunk centered aroung the CG based on how many species actually
    # have the CG conserved, since we can't design for an SNV in there. Also returns the list of positions where we
    # need an SNV. If multiple SNVs are needed at a position, that position is featured multiple times in the list.
    maxSpeciesSubset = set()
    for species in curChunk[min(curChunk.keys()) + 50].keys():
        if curChunk[min(curChunk.keys()) + 50][species] + curChunk[min(curChunk.keys()) + 51][species] == "CG":
            maxSpeciesSubset.add(species)

    SNVPossibilities = []
    maxSpeciesSubsetCopy = set(maxSpeciesSubset)
    different = defaultdict(set)
    for species in maxSpeciesSubset:
        mismatches = []
        for pos in range(startPos, endPos):
            if curChunk[pos][species] != curChunk[pos][refSpecies]:
                mismatches.append((pos, curChunk[pos][species]))
                if (len(mismatches) > numSNV):
                    break
        if len(mismatches) <= numSNV:
            for i in mismatches:
                different[i[0]].add(i[1])
        else:
            maxSpeciesSubsetCopy.remove(species)

    for pos in range(startPos, endPos):
        SNVPossibilities += ([pos] * len(different[pos]))
    maxSpeciesSubset = set(maxSpeciesSubsetCopy)

    return (maxSpeciesSubset, SNVPossibilities)

def scoreProbe(probeStrand, probeType, numSNV, curChunk, infiniumType):
    # print(probeStrand, probeType, infiniumType)
    # For each probe, go through the 4 options for start and end and get the largest set of species we can
    if (probeStrand == "TOP"):
        if (probeType == "O"):
            if (infiniumType == 1):
                startPos = min(curChunk.keys()) + 1
            elif (infiniumType == 2):
                startPos = min(curChunk.keys())
            else:
                print("Wrong infinium probe type. Must be 1 or 2.")
        elif (probeType == "C"):
            if (infiniumType == 1):
                startPos = min(curChunk.keys()) + 50
            elif (infiniumType == 2):
                startPos = min(curChunk.keys()) + 51
            else:
                print("Wrong infinium probe type. Must be 1 or 2.")
        else:
            print("Wrong probe type found in column 7.")
            exit(1)
    elif (probeStrand == "BOT"):
        if (probeType == "O"):
            if (infiniumType == 1):
                startPos = min(curChunk.keys()) + 51
            elif (infiniumType == 2):
                startPos = min(curChunk.keys()) + 52
            else:
                print("Wrong infinium probe type. Must be 1 or 2.")
        elif (probeType == "C"):
            if (infiniumType == 1):
                startPos = min(curChunk.keys()) + 2
            elif (infiniumType == 2):
                startPos = min(curChunk.keys()) + 1
            else:
                print("Wrong infinium probe type. Must be 1 or 2.")
        else:
            print("Wrong probe type found in column 7.")

    maxSpecies, SNVpositions = getFeasibleSpecies(curChunk, startPos, startPos + 50, "hg19", numSNV)
    probeSequence = ""
    for pos in range(startPos, startPos + 50):
        probeSequence += curChunk[pos]["hg19"]
    bestSpecies = set()
    bestSNVchoice = []

    for positionComb in combinations(SNVpositions, min(numSNV, len(SNVpositions))):
        curMaxSpecies, curSNVchoice = getMaxSpecies(positionComb, maxSpecies, curChunk, startPos, startPos + 50, "hg19", pt=True)
        #if (positionComb == (30864656, 30864659, 30864659)):
        #    print("But out here")
        #    print(curMaxSpecies)
        #    print(len(curMaxSpecies))

        if (len(curMaxSpecies) > len(bestSpecies)):
            bestSpecies = curMaxSpecies
            bestSNVchoice = curSNVchoice
    # returns the largest list of species we can cover, the best combination of SNVs and the best SNV choice
    return bestSpecies, bestSNVchoice, probeSequence, startPos + 1, startPos + 50

def main():
    if len(sys.argv) != 4:
        print("Usage: python getCGConservationScores.py <alignment features file> <probes file> <output file>")
        exit(1)
    alignmentFile = gzip.open(sys.argv[1], 'rt')
    probesFile = gzip.open(sys.argv[2], 'rt')
    ofile = gzip.open(sys.argv[3], 'wt')

    isMammal = {"hg19": 1,"panTro4": 1,"gorGor3": 1,"ponAbe2": 1,"nomLeu3": 1,"rheMac3": 1,"macFas5": 1,"papHam1": 1,"chlSab1": 1,"calJac3": 1,"saiBol1": 1,"otoGar3": 1,"tupChi1": 1,"speTri2": 1,"jacJac1": 1,"micOch1": 1,"criGri1": 1,"mesAur1": 1,"mm10": 1,"rn5": 1,"hetGla2": 1,"cavPor3": 1,"chiLan1": 1,"octDeg1": 1,"oryCun2": 1,"ochPri3": 1,"susScr3": 1,"vicPac2": 1,"camFer1": 1,"turTru2": 1,"orcOrc1": 1,"panHod1": 1,"bosTau7": 1,"oviAri3": 1,"capHir1": 1,"equCab2": 1,"cerSim1": 1,"felCat5": 1,"canFam3": 1,"musFur1": 1,"ailMel1": 1,"odoRosDiv1": 1,"lepWed1": 1,"pteAle1": 1,"pteVam1": 1,"myoDav1": 1,"myoLuc2": 1,"eptFus1": 1,"eriEur2": 1,"sorAra2": 1,"conCri1": 1,"loxAfr3": 1,"eleEdw1": 1,"triMan1": 1,"chrAsi1": 1,"echTel2": 1,"oryAfe1": 1,"dasNov3": 1,"monDom5": 1,"sarHar1": 1,"macEug2": 1,"ornAna1": 1, "falChe1": 0,"falPer1": 0,"ficAlb2": 0,"zonAlb1": 0,"geoFor1": 0,"taeGut2": 0,"pseHum1": 0,"melUnd1": 0,"amaVit1": 0,"araMac1": 0,"colLiv1": 0,"anaPla1": 0,"galGal4": 0,"melGal1": 0,"allMis1": 0,"cheMyd1": 0,"chrPic1": 0,"pelSin1": 0,"apaSpi1": 0,"anoCar2": 0,"xenTro7": 0,"latCha1": 0,"tetNig2": 0,"fr3": 0,"takFla1": 0,"oreNil2": 0,"neoBri1": 0,"hapBur1": 0,"mayZeb1": 0,"punNye1": 0,"oryLat2": 0,"xipMac1": 0,"gasAcu1": 0,"gadMor1": 0,"danRer7": 0,"astMex1": 0,"lepOcu1": 0,"petMar2": 0}

    # Read first 100 bases from alignment file
    refSpecies = "hg19"
    curChunk = defaultdict() # dictionary of dictionaries of the nucleotide at each position for each species

    headerLine = alignmentFile.readline().rstrip().split(",")
    for i in range(102):
        splitLine = alignmentFile.readline().rstrip().split(",")
        curPos = int(splitLine[0])
        curChunk[curPos] = {}

        for j in range(2, len(splitLine)):
            if (isMammal[headerLine[j]]):
                curChunk[curPos][headerLine[j]] = splitLine[j]

    end = False
    count = 0
    for line in probesFile:
        count += 1
        if (count % 10000 == 0):
            print(count)
        splitLine = line.strip().split()
        #print(splitLine)
        coord = int(splitLine[4]) - 1
        designScore = float(splitLine[22])
        tempStrand = splitLine[15].split("_")[-1]
        if (tempStrand == "F"):
            probeStrand = "TOP"
        elif (tempStrand == "R"):
            probeStrand = "BOT"
        else:
            print ("Unknown probe strand found for ", line)
            exit(1)
        probeType = splitLine[17]

        # get to position in in the alignment that is centered on the current CG site
        while (curPos - 51 < coord):
            if (curPos % 10000 == 0):
                print(curPos)
            newLine = alignmentFile.readline()
            if (not newLine):
                end = True
                break # add what to do in this case
            newSplitLine = newLine.rstrip().split(",")
            curPos = int(newSplitLine[0])
            curChunk[curPos] = {}
            curChunk.pop(min(curChunk.keys()))

            for j in range(2, len(newSplitLine)):
                if (isMammal[headerLine[j]]):
                    curChunk[curPos][headerLine[j]] = newSplitLine[j]

        if (max(curChunk.keys()) - min(curChunk.keys()) == 101):
            printEL = False
            if (end):
                break
            # Infinium 1 probe
            if designScore >= 0.6:
                numSNV = 3
            elif designScore >= 0.5:
                numSNV = 2
            elif designScore >= 0.4:
                numSNV = 1
            elif designScore >= 0.3:
                numSNV = 0
            else:
                numSNV = -1

            if numSNV >= 0:
                bestSpecies, bestSNVchoice, probeSequence, probeStart, probeEnd = scoreProbe(probeStrand, probeType, numSNV, curChunk, infiniumType = 1)
                if (len(list(bestSpecies)) > 0):
                    ofile.write(line.strip() + "\t" + ",".join(list(bestSpecies)) + "\t" + str(probeStart) + "\t" + str(probeEnd) + "\t" + probeSequence)
                    if (len(bestSNVchoice) != 0):
                        # write position of SNPs, reference nucleotide and alternate nucleotide
                        ofile.write("\t" + ",".join([str(x[0] + 1) for x in bestSNVchoice]) + "\t" + ",".join([str(x[2]) for x in bestSNVchoice]) + "\t" + ",".join([str(x[1]) for x in bestSNVchoice]))
                        # write which species they each target
                        for x in bestSNVchoice:
                            ofile.write("\t" + ",".join(x[3]))
                    else:
                        ofile.write("\t-\t-\t-")
                    for x in range(len(bestSNVchoice), 3):
                        ofile.write("\t-")
                    printEL = True
            #print(bestSpecies, bestSNVchoice)

            #Infinium 2 probe
            numCpG = int(splitLine[-4])
            if numCpG > 3:
                numSNV = -1
            else:
                if designScore >= 0.6:
                    numSNV = 3 - numCpG
                elif designScore >= 0.5:
                    numSNV = 2 - numCpG
                elif designScore >= 0.4:
                    numSNV = 1 - numCpG
                elif designScore >= 0.3:
                    numSNV = 0
                else:
                    numSNV = -1

            # print("numSNV = ", numSNV)
            if numSNV >= 0:
                bestSpecies, bestSNVchoice, probeSequence, probeStart, probeEnd = scoreProbe(probeStrand, probeType, numSNV, curChunk, infiniumType=2)
                if (len(list(bestSpecies)) > 0):
                    ofile.write("\t" + ",".join(list(bestSpecies)) + "\t" + str(probeStart) + "\t" + str(probeEnd) + "\t" + probeSequence)
                    if (len(bestSNVchoice) != 0):
                        ofile.write("\t" + ",".join([str(x[0] + 1) for x in bestSNVchoice]) + "\t" + ",".join([str(x[2]) for x in bestSNVchoice]) + "\t" + ",".join([str(x[1]) for x in bestSNVchoice]))
                        # write which species they each target
                        for x in bestSNVchoice:
                            ofile.write("\t" + ",".join(x[3]))
                    else:
                        ofile.write("\t-\t-\t-")
                    for x in range(len(bestSNVchoice), 3):
                        ofile.write("\t-")
                    printEL = True

            if (printEL):
                ofile.write("\n")

    alignmentFile.close()
    ofile.close()

main()

