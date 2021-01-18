import gzip
import argparse
from collections import deque

def add_species_coords(main_args):
    probes_file = gzip.open(main_args.probes_file, 'rt')
    maf_seq_file = gzip.open(main_args.maf_seq_file, 'rt')
    output_file = gzip.open(main_args.output_file, 'wt')

    isMammal = {"hg19": 1,"panTro4": 1,"gorGor3": 1,"ponAbe2": 1,"nomLeu3": 1,"rheMac3": 1,"macFas5": 1,"papHam1": 1,"chlSab1": 1,"calJac3": 1,"saiBol1": 1,"otoGar3": 1,"tupChi1": 1,"speTri2": 1,"jacJac1": 1,"micOch1": 1,"criGri1": 1,"mesAur1": 1,"mm10": 1,"rn5": 1,"hetGla2": 1,"cavPor3": 1,"chiLan1": 1,"octDeg1": 1,"oryCun2": 1,"ochPri3": 1,"susScr3": 1,"vicPac2": 1,"camFer1": 1,"turTru2": 1,"orcOrc1": 1,"panHod1": 1,"bosTau7": 1,"oviAri3": 1,"capHir1": 1,"equCab2": 1,"cerSim1": 1,"felCat5": 1,"canFam3": 1,"musFur1": 1,"ailMel1": 1,"odoRosDiv1": 1,"lepWed1": 1,"pteAle1": 1,"pteVam1": 1,"myoDav1": 1,"myoLuc2": 1,"eptFus1": 1,"eriEur2": 1,"sorAra2": 1,"conCri1": 1,"loxAfr3": 1,"eleEdw1": 1,"triMan1": 1,"chrAsi1": 1,"echTel2": 1,"oryAfe1": 1,"dasNov3": 1,"monDom5": 1,"sarHar1": 1,"macEug2": 1,"ornAna1": 1, "falChe1": 0,"falPer1": 0,"ficAlb2": 0,"zonAlb1": 0,"geoFor1": 0,"taeGut2": 0,"pseHum1": 0,"melUnd1": 0,"amaVit1": 0,"araMac1": 0,"colLiv1": 0,"anaPla1": 0,"galGal4": 0,"melGal1": 0,"allMis1": 0,"cheMyd1": 0,"chrPic1": 0,"pelSin1": 0,"apaSpi1": 0,"anoCar2": 0,"xenTro7": 0,"latCha1": 0,"tetNig2": 0,"fr3": 0,"takFla1": 0,"oreNil2": 0,"neoBri1": 0,"hapBur1": 0,"mayZeb1": 0,"punNye1": 0,"oryLat2": 0,"xipMac1": 0,"gasAcu1": 0,"gadMor1": 0,"danRer7": 0,"astMex1": 0,"lepOcu1": 0,"petMar2": 0}

    end_ID = {}
    end_probes = {}
    species_dict = {}
    first_line = True
    for line in probes_file:
        split_line = line.strip().split()
        probe_ID = split_line[0]
        end_probes[probe_ID] = int(split_line[26])
        end_ID[int(split_line[26])] = probe_ID
        species_dict[probe_ID] = split_line[24].split(',')

    header_line = maf_seq_file.readline().strip().split(',')
    output_file.write('probeID\thg19_start\thg19end\t')
    mammal_list = []
    for i in header_line:
        species = i.strip().split('_')[1]
        if species not in mammal_list:
            if isMammal[species] and species != 'hg19':
                output_file.write(species + '_start\t' + species + '_end\t' + species + '_has_insertion\t' + species + \
                                  '_is_split\t')
                mammal_list.append(species)
    output_file.write('\n')

    count_insertion_issues = 0
    count_probes_done = 0

    pos_reference = header_line.index('pos_hg19')
    last_50_lines = deque()
    for line in maf_seq_file:
        split_line = line.strip().split(',')
        last_50_lines.append(split_line)
        pos_hg19 = split_line[pos_reference].split(':')[1]
        #print(pos_hg19)

        if pos_hg19 == 'd': # if deletion at current position in reference continue
            continue
        end_human = int(last_50_lines[-1][pos_reference].split(':')[1])
        start_human = int(last_50_lines[0][pos_reference].split(':')[1])

        while (end_human - start_human >= 50):
            last_50_lines.popleft()
            while last_50_lines[0][pos_reference].split(':')[1] == 'd':
                last_50_lines.popleft()
            end_human = int(last_50_lines[-1][pos_reference].split(':')[1])
            start_human = int(last_50_lines[0][pos_reference].split(':')[1])
        pos_hg19 = int(pos_hg19)

        if pos_hg19 in end_ID:
            #for prev_line in last_50_lines: # have to make sure going over last 50 lines uninterrupted by deletion
            # in human
            #    print(prev_line[pos_reference + 2], sep='', end = '')
            #print('\n')
            cur_probe = end_ID[pos_hg19]
            output_file.write(cur_probe + '\t' + last_50_lines[0][pos_reference] + '\t' + last_50_lines[-1][
                pos_reference] + '\t')

            found_ins = False
            for sp in mammal_list:
                if sp in species_dict[cur_probe]:
                    idx_species = header_line.index('pos_' + sp)
                    strand_species = header_line.index('strand_' + sp)
                    #print(sp)
                    #print(idx_species)
                    has_ins = False
                    for prev_line in last_50_lines: # have to make sure going over last 50 lines uninterrupted by deletion in human
                    #    print(prev_line[idx_species + 2], sep='', end = '')
                        if prev_line[pos_reference] == 'd' and prev_line[idx_species] != 'd':
                            has_ins = True
                            found_ins = True
                    #print('\n')
                    is_split_probe = False # have to check if coordinates in other species are continuous
                    for i in range(1, len(last_50_lines) - 1):
                        last_pos = last_50_lines[i - 1][idx_species].split(':')[1]
                        cur_pos = last_50_lines[i][idx_species].split(':')[1]
                        if last_pos != 'd' and cur_pos != 'd':
                            if last_pos.isdigit() and cur_pos.isdigit():
                                last_pos = int(last_pos)
                                cur_pos = int(cur_pos)
                                if cur_pos != last_pos + 1 and cur_pos != last_pos - 1:
                                    is_split_probe = True
                        else:
                            is_split_probe=True
                    start_sp = last_50_lines[0][idx_species]
                    end_sp = last_50_lines[-1][idx_species]
                    end_pos = end_sp.split(':')[1]
                    start_pos = start_sp.split(':')[1]
                    if end_pos < start_pos:
                        start_sp, end_sp = end_sp, start_sp
                    output_file.write(str(start_sp + 1) + "\t" + str(end_sp + 1) + "\t" + str(has_ins) + '\t' + str(
                        is_split_probe) + '\t')
                else:
                    output_file.write('-\t-\t-\t-\t')
            if found_ins:
                count_insertion_issues += 1
            output_file.write('\n')
            count_probes_done += 1
            print(count_probes_done, '/', len(end_probes.keys()))
            if count_probes_done == len(end_probes.keys()):
                break

    output_file.close()
    print(count_insertion_issues, " probes have insertions in other species not accounted for.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Adds in coordinates of methylation probes in other species, 
    because initial file submitted to Bret did not contain this.''')

    parser.add_argument('-p', '--probesFile', required=True, dest='probes_file',
                        help='File with probes that was submitted to Bret')
    parser.add_argument('-m', '--MAF_seq_file', required=True, dest='maf_seq_file',
                        help='File containing sequence info extracted from MAF files, containing all positions in '
                             'each species.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file to write the probe information to.')

    args = parser.parse_args()

    add_species_coords(args)
