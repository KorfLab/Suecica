#!/usr/bin/env python3
import argparse, random, re
import SAM

def Command_line():
	parser = argparse.ArgumentParser(description="DistantDupSearch_GenerateRandomRegions.py generates a random list of regions which contains the same number of regions as in the provided data file of duplicated region, and the randomly generated regions are also of equal length. The file this script outputs is to be used with 'DistantDupSearchWithMatePairs.py' in order to determine whether or not the number of distant mate pairs found for a list of duplicated regions is comparable to a list containing a random selection of regions.")
	parser.add_argument('-d', default="Data/As_Duplicated_Regions.txt", type=str, help='File containing a list of all DupHMM duplicated regions', metavar='DuplicatedRegions')
	parser.add_argument('-o', default="Data/As_Random_Regions.txt", type=str, help='Output file for list of randomly selected regions', metavar='OutputDir')
	
	args = parser.parse_args()
	dup_file = args.d
	o_file = args.o
	
	return(dup_file, o_file)

dup_file, o_file = Command_line()
chromosome_len = SAM.Chr_Lengths("Data/weigel_col-0/SRR013327,SRR013328_aln.sam.gz")
del chromosome_len["chloroplast"]
del chromosome_len["mitochondria"]

region_sizes = []
with open(dup_file) as infile:
	for line in infile:
		pattern = r"(Chr\d{1}):(\d+)-(\d+)"
		recomp = re.compile(pattern)
		match = recomp.match(line)
		if match:
			(chromosome, s_pos, e_pos) = match.groups()
			dist = int(e_pos) - int(s_pos)
			region_sizes.append(dist)

rand_regions = []
for dist in region_sizes:
	# Randomly generate regions and ensure that no two regions are overlapping
	continue_past_loop = False
	while continue_past_loop == False:
		chromosome = random.sample(chromosome_len.keys(),1)[0]
		max_e_pos = chromosome_len[chromosome] - dist
		s_pos = random.randint(1, max_e_pos)
		e_pos = s_pos + dist
		no_overlap = True
		for region in rand_regions:
			temp_chromosome, temp_s_pos, temp_e_pos = region
			if chromosome == temp_chromosome:
				if (s_pos > temp_s_pos and s_pos < temp_e_pos) or (e_pos > temp_s_pos and e_pos < temp_e_pos):
					no_overlap = False
					break
		if no_overlap == True: continue_past_loop = True
	rand_regions.append((chromosome, s_pos, e_pos))

with open(o_file, 'w') as outfile:
	for region in sorted(rand_regions):
		chromosome, s_pos, e_pos = region
		outline = str(chromosome) + ":" + str(s_pos) + "-" + str(e_pos) + "\n"
		outfile.write(outline)