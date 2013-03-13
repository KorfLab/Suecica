#!/usr/bin/env python3
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import argparse, gzip, re, os
from collections import Counter

def command_line():
    parser = argparse.ArgumentParser(description="DistantDupSearchWithMatePairs.py searches for reads which (1) map within a duplicated region found by DupHMM, and (2) which have mate pairs that map either to another chromosome or to the same chromosome but outside the original read's duplicated region. Statistics are also calculated for the fraction of reads with mate pairs within the same duplicated region vs. outside that region.")
    parser.add_argument('-f', default="Data/Sue/sue1/sliding_window_preprocess/set5/sue_set5_sw_aln.sam.gz", type=str, help='SAM file with forward and reverse reads.', metavar='SAMFile')
    parser.add_argument('-r', default=1000, type=int, help='Range around duplicated region within which mate pairs should not be considered outside of a duplicated region. Default is 1000 bps', metavar='DupRange')
    parser.add_argument('-d', default="Data/As_Random_Regions.txt", type=str, help='File containing a list of all DupHMM duplicated regions', metavar='DuplicatedRegions')
    parser.add_argument('-o', default="Output/Sue1Set5RandomRegionMatePairs/", type=str, help='Output directory for three files: (1) Reads w/ mate pairs that map to a different chromosome (Diff_Chr.txt), (2) Reads w/ mate pairs that map to the same chromosome but outside the first read\'s duplicated region (Same_Chr.txt), and (3) a file containing summary statistics (Summary_Stats.txt)', metavar='OutputDir')
    #parser.add_argument('-d', default="Data/As_Duplicated_Regions.txt", type=str, help='File containing a list of all DupHMM duplicated regions', metavar='DuplicatedRegions')
    #parser.add_argument('-o', default="Output/Sue1Set5DuplicatedRegionMatePairs/", type=str, help='Output directory for three files: (1) Reads w/ mate pairs that map to a different chromosome (Diff_Chr.txt), (2) Reads w/ mate pairs that map to the same chromosome but outside the first read\'s duplicated region (Same_Chr.txt), and (3) a file containing summary statistics (Summary_Stats.txt)', metavar='OutputDir')
    
    args = parser.parse_args()
    reads_filename = args.f
    dup_filename = args.d
    dup_range = args.r
    outdir = args.o
    
    return(reads_filename, dup_filename, dup_range, outdir)

reads_filename, dup_filename, dup_range, outdir = command_line()

if not os.path.exists(outdir): os.makedirs(outdir)

# Grab list of duplicated regions as determined by DupHMM
dup_locations = {}
dup_region_dict = Counter()
pattern = r'(Chr\d{1}):(\d+)-(\d+)'
recomp = re.compile(pattern)
with open(dup_filename) as infile:
    for line in infile:
        match = recomp.search(line)
        if match:
            chromosome, s_pos, e_pos = match.groups()
            s_pos = int(s_pos)
            e_pos = int(e_pos)
            dup_region_dict[(chromosome, s_pos, e_pos)] = 0
            try:
                dup_locations[chromosome].append((s_pos, e_pos))
            except:
                dup_locations[chromosome] = [(s_pos, e_pos)]

# Open SAM file containing reads w/ mate pairs
# Save reads for output which have mate pairs mapping to different chromosomes
# and mate pairs mapping to the same chromosome but outside the original read's
# duplicated region
mate_check = {}
diff_chr_outlines = []
same_chr_outlines = []

diff_chr_reads = 0
same_chr_reads = 0
total_dup_reads = 0

diff_chr_o_fn = os.path.join(outdir, "Diff_Chr.txt")
same_chr_o_fn = os.path.join(outdir, "Same_Chr.txt")
stats_o_fn = os.path.join(outdir, "Summary_Stats.txt")

with gzip.open(reads_filename, 'rb') as infile:
    for line in infile:
        line = line.decode('utf-8')
        line = line.split("\t")
        if len(line[0]) > 3:
            read_name = line[0]
            read_num = int(read_name[-8:-7])
            chromosome = line[2]
            if chromosome == "*" or chromosome[:3] != "Chr": continue
            s_pos = int(line[3])
            
            # Check if this read is within a duplicated region. If so, store it so its mate pair on the next line can be checked
            if read_num == 1:
                e_pos = s_pos + len(line[9])
                for (dup_region_s, dup_region_e) in dup_locations[chromosome]:
                    if (s_pos > dup_region_s and s_pos < dup_region_e) or (e_pos > dup_region_s and e_pos < dup_region_e):
                        mate_check[read_name] = [chromosome, s_pos, e_pos, dup_region_s, dup_region_e]
            
            # Check to see if this read has a mate pair located within a duplicated region
            elif read_num == 2:
                pair_read_name = read_name[:-8] + "1" + read_name[-7:]
                if pair_read_name in mate_check:
                    # Duplicated region mate pair exists. See if they're located on different chromosome or the same chromosome but distant
                    pair_chromosome, pair_s_pos, pair_e_pos, dup_region_s, dup_region_e = mate_check[pair_read_name]
                    e_pos = s_pos + len(line[9])
                    if pair_chromosome != chromosome:
                        # Mate pair is located on a different chromosome
                        outline = "\t".join([pair_read_name, pair_chromosome, str(pair_s_pos), str(pair_e_pos), read_name, chromosome, str(s_pos), str(e_pos)]) + "\n"
                        diff_chr_outlines.append(outline)
                        diff_chr_reads += 2
                        total_dup_reads += 2
                    elif s_pos > dup_region_e + dup_range or e_pos < dup_region_s - dup_range:
                        # Mate pair is located outside of other read's duplicated region
                        distance = 0
                        if s_pos > dup_region_e:
                            distance = s_pos - dup_region_e
                        elif dup_region_s > e_pos:
                            distance = dup_region_s - e_pos
                        outline = "\t".join([pair_read_name, pair_chromosome, str(pair_s_pos), str(pair_e_pos), str(dup_region_s), str(dup_region_e), read_name, chromosome, str(s_pos), str(e_pos), str(distance)]) + "\n"
                        same_chr_outlines.append(outline)
                        same_chr_reads += 2
                        total_dup_reads += 2
                    elif (s_pos >= dup_region_s - dup_range and s_pos <= dup_region_e + dup_range) or (e_pos <= dup_region_e + dup_range and e_pos >= dup_region_s - dup_range):
                        # Mate pair is located within the same duplicated region
                        total_dup_reads += 2
                    else:
                        # Something happened that shouldn't have
                        print_outline = "\t".join([pair_read_name, pair_chromosome, str(pair_s_pos), str(pair_e_pos), str(dup_region_s), str(dup_region_e), read_name, chromosome, str(s_pos), str(e_pos)])
                        print(print_outline)
                    del mate_check[pair_read_name]

with open(diff_chr_o_fn, 'w') as outfile:
    outline = "Read1\tRead1_Chr\tRead1_Start\tRead1_End\tRead2\tRead2_Chr\tRead2_Start\tRead2_End\n"
    outfile.write(outline)
    for outline in diff_chr_outlines:
        outfile.write(outline)

with open(same_chr_o_fn, 'w') as outfile:
    outline = "Read1\tRead1_Chr\tRead1_Start\tRead1_End\tRead1_DupRegion_Start\tRead1_DupRegion_End\tRead2\tRead2_Chr\tRead2_Start\tRead2_End\tRead2DistanceFromDupRegion\n"
    outfile.write(outline)
    for outline in same_chr_outlines:
        outfile.write(outline)

with open(stats_o_fn, 'w') as outfile:
    same_dup_region_reads = total_dup_reads - same_chr_reads - diff_chr_reads
    outline = ""
    if dup_range == 0:
        outline = "Duplicated region reads with mate pairs in the same duplicated region:\t\t" + str(same_dup_region_reads) + "\t(" + str(round(same_dup_region_reads / total_dup_reads * 100, 2)) + "%)\n"
    else:
        outline = "Duplicated region reads with mate pairs in the same duplicated region, +/- " + str(dup_range) + " bps:\t" + str(same_dup_region_reads) + "\t(" + str(round(same_dup_region_reads / total_dup_reads * 100, 2)) + "%)\n"
    outfile.write(outline)
    outline = "Duplicated region reads with mate pairs on same chromosome, outside dup region:\t\t\t" + str(same_chr_reads) + "\t(" + str(round(same_chr_reads / total_dup_reads * 100, 2)) + "%)\n"
    outfile.write(outline)
    outline = "Duplicated region reads with mate pairs on a different chromosome:\t\t\t\t\t\t" + str(diff_chr_reads) + "\t(" + str(round(diff_chr_reads / total_dup_reads * 100, 2)) + "%)\n"
    outfile.write(outline)

# Create file listing which of the duplicated regions contains reads w/ distant mate pairs
dup_dist_filename = dup_filename[:-4] + "_With_Distant_Mate_Pairs.txt"
for line in diff_chr_outlines:
    line = line.split("\t")
    chromosome = line[1]
    s_pos = int(line[2])
    e_pos = int(line[3])
    for entry in dup_locations[chromosome]:
        dup_s_pos, dup_e_pos = entry
        if (s_pos > dup_s_pos and s_pos < dup_e_pos) or (e_pos > dup_s_pos and e_pos < dup_e_pos):
            dup_region_dict[(chromosome, dup_s_pos, dup_e_pos)] += 1

for line in same_chr_outlines:
    line = line.split("\t")
    chromosome = line[1]
    s_pos = int(line[4])
    e_pos = int(line[5])
    for entry in dup_locations[chromosome]:
        dup_s_pos, dup_e_pos = entry
        if s_pos == dup_s_pos and s_pos == dup_e_pos:
            dup_region_dict[(chromosome, dup_s_pos, dup_e_pos)] += 1

with open(dup_dist_filename, 'w') as outfile:
    outline = "Dup_Region\tDistantMatePairs\n"
    outfile.write(outline)
    for entry in sorted(dup_region_dict.items(), key = lambda x: (x[0][0], x[0][1], x[0][2])):
        outline = str(entry[0][0]) + ":" + str(entry[0][1]) + "-" + str(entry[0][2]) + "\t" + str(entry[1]) + "\n"
        outfile.write(outline)