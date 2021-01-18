import gzip

if __name__ == "__main__":
    array_file = open("../mammalianChip_36K_probes.tsv", "r")
    scored_probes_file = gzip.open("../allProbes_scored_GW.tsv.gz", "rt")
    output_file = gzip.open("../36K_probes_SNV_species.tsv.gz", "wt")

    # Construct dictionary of probes submitted to Bret
    inf_type_dict = {}
    forward_type_dict = {}
    line_dict = {}
    chrom_dict = {}
    pos_dict = {}
    for line in array_file:
        if line[:2] == "cg":
            split_line = line.strip().split("\t")
            probe_ID = split_line[0].split("_")[0].split("-")[0]
            inf_type = split_line[31]
            forward_type = split_line[0].split("_")[-1]
            chrom_dict[probe_ID] = split_line[4]
            pos_dict[probe_ID] = split_line[5]

            if probe_ID in line_dict:
                print("WTF")
                print(line_dict[probe_ID])
                print(line)
            inf_type_dict[probe_ID] = inf_type
            forward_type_dict[probe_ID] = forward_type
            line_dict[probe_ID] = split_line
        else:
            header_line = line.strip().split("\t")
    output_file.write("\t".join(header_line[:31]) + "\ttarget_SNV1\ttarget_SNV2\ttarget_SNV3\t" +
                      "\t".join(header_line[31:]) + "\n")

    print("Total number of probes on array: ", len(line_dict.keys()))

    # Go through list of all probes and add SNP info
    probes_done_dict = {}
    probes_done = 0

    for line in scored_probes_file:
        split_line = line.strip().split("\t")
        probe_ID = split_line[0]
        forward_type = split_line[15]
        is_converted = (split_line[17] == "C")

        chrom = split_line[3]
        pos = split_line[4]
        if (probe_ID in line_dict) and (forward_type_dict[probe_ID] == forward_type) and (is_converted) and \
                chrom_dict[probe_ID] == chrom and pos_dict[probe_ID] == pos:
            probes_done += 1
            if probe_ID in probes_done_dict:
                print("WTF")
                print(probe_ID)
                print(line_dict[probe_ID])
                print(line)
                exit(1)
            probes_done_dict[probe_ID] = 1
            if probes_done % 100 == 0:
                print("Probes done so far: ", probes_done)
            inf_type = inf_type_dict[probe_ID]

            output_file.write("\t".join(line_dict[probe_ID][:31]) + "\t")
            if inf_type == "Inf1":
                output_file.write(split_line[30] + "\t" + split_line[31] + "\t" + split_line[32] + "\t")
            else:
                output_file.write(split_line[40] + "\t" + split_line[41] + "\t" + split_line[42] + "\t")
            output_file.write("\t".join(line_dict[probe_ID][31:]) + "\n")

    probes_not_done = open("../probes_not_done.tsv", "w")
    for probe_ID in line_dict:
        if probe_ID not in probes_done_dict:
            probes_not_done.write("\t".join(line_dict[probe_ID]) + "\n")

    output_file.close()
