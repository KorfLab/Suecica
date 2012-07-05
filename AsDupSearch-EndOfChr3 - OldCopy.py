#!/usr/bin/env python3.2
import re, time, decimal; from scipy.stats import fisher_exact, scoreatpercentile; from copy import deepcopy; from math import sqrt;
import sys; sys.path.append('/home/mattchat/Packages/')
from GenomeTools import GFF, Pileup

# This script searches for genes in A. suecica which underwent a single ancestral duplication event (shortly after the A. suecica
# species arose), which was then soon followed by a mutation in one of the duplicates. Specifically, this script narrows down on
# the end of At-Chr.3/Aa-Chr.5, where the end of At-Chr.3 rearranged onto the end of Aa-Chr.5
# Script will focus on At-Chr. 3 genes which map at or beyond AT3G62280. This end of the chromosome was deleted and
# transposed to Aa-Chr. 5, and is known to have undergone duplication. This script's ability to detect SUNs will give
# an idea of how recently the duplication of genes within this chromosomal rearrangement region occurred.
# This is performed by:
# 1)  Looking for positions within genic regions that have a 50/50 consensus/one-nucleotide SNP distribution, which is
#		 indicative of an ancestral duplication->mutation event.
# 1a) If such a position maps to A. thaliana, it is compared to the At parent to ensure that the 50/50 distribution is truly due to
#		 an ancestral duplication->mutation event, rather than simply being inherited heterozygosity. If the position maps to
#		 A. lyrata, then it is compared to the Aa parent.
# 2)  Once a gene is determined to have at least one SUN that is not based on inherited heterozygosity, that gene will have its
#		 read counts in As compared to read counts in At. A 2:1 ratio is expected to be found. As and At must first have their
#		 read counts be normalized, however, as the read coverage is not equal between the two (At mean ~ 37, As mean ~ 145).
#		 Values are normalized by dividing a gene's read counts by the 75th percentile genic read counts for the species being considered.
#		 After normalizing the read count values, if a SUN-containing gene has approximately a doubled read count in As when
#		 compared to At, there is substantial evidence that this gene underwent a single ancestral duplication event, and such a gene
#		 would be added to a list of candidate duplicated genes.
# 2a)  Before being added to a list of candidate duplicated genes, the gene is checked first to ensure it does not better fit a 3:1
#		 normalized read ratio. Presently only interested in genes which have undergone a single duplication event.
# 3)  After the script has run through every gene and generated a list of candidate duplicates, these are output in descending
#		 order of total per-gene SUN counts.
#

def get_genes():
	#Chr1	TAIR9	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
	pattern = r"^Chr(\d{1})\tTAIR9\tgene\t(\d*)\t(\d*)\t\.\t[+-]\t\.\tID=(AT.{1}G\d{5});.*$"
	recomp = re.compile(pattern)
	thal_filename = "../t9_sorted.gff"
	gene_nuc_dict = {}; gene_dict = {}
	ignore_genes = True
	end_ignore = 0
	
	with open(thal_filename) as infile:
		for line in infile:
			match = recomp.match(line)
			if match:
				chr_num = int(match.group(1)); start_pos = int(match.group(2)); end_pos = int(match.group(3)); gene_name = str(match.group(4));
				if ignore_genes == False and chr_num != 3: # 1) No longer searching in Chr. 3 beyond AT3G62280, so resume ignoring genes
					ignore_genes = True
				if gene_name == "AT3G62280" and ignore_genes == True: # 1) Looking at genes at or beyond AT3G62280 on At-Chr. 3: the site of At-Chr. 3 deletion ==> Aa-Chr. 5 concatenation
					ignore_genes = False
					end_ignore = start_pos
				if ignore_genes == False:
					gene_dict["At", gene_name] = [chr_num, start_pos, end_pos, 0]
					for i in range(start_pos, end_pos+1):
						gene_nuc_dict[("At", chr_num, i)] = gene_name
	
	#scaffold_1	JGI	exon	5631	6198		-		name "fgenesh2_kg.1__2__AT1G02190.2"; transcriptId 470048
	pattern2 = r"^scaffold_([1-8])\tJGI\t.+\t(\d+)\t(\d+).*(AT[1-5]G\d{5}).*$"
	recomp2 = re.compile(pattern2)
	lyr_filename = "../a_lyrata_6_igv.gff"
	pos_list = []
	prev_chr = 1
	prev_gene = ""

	with open(lyr_filename) as infile:
		for line in infile:
			match = recomp2.match(line)
			if match:
				chr_num = int(match.group(1)); start = int(match.group(2)); end = int(match.group(3)); gene_name = str(match.group(4))
				if chr_num == 5 and start >= 20000000: # Approx. the At-Chr. 3 concatenated region. Appears to start after 20500000, but allowing extra bases (starting at 20000000) to ensure no genes are missed
					if prev_gene == "":
						prev_gene = gene_name
					if gene_name != prev_gene:
						if ("At", prev_gene) in gene_dict:
							max_start = max_end = 0
							for pair in pos_list:
								(temp_s, temp_e) = pair
								if int(temp_s) > int(max_start):
									max_start = temp_s
								if int(temp_e) > int(max_end):
									max_end = temp_e
							pos_list = []
							for i in range(1,15):
								if ("Al", gene_name, i) not in gene_dict:
									gene_dict[("Al", gene_name, i)] = [prev_chr, max_start, max_end, 0]
									for j in range(max_start, max_end+1):
										gene_nuc_dict[("Al", prev_chr, j)] = (gene_name, i)
									break
						prev_gene = gene_name
						prev_chr = chr_num
					else:
						pos_list.append((start, end))

	#with open("gene_nuc_dict.txt", 'w') as out_test:
	#	for entry in sorted(gene_nuc_dict.items(), key = lambda entry:(entry[0][0],entry[0][1],entry[0][2])):
	#		line = str(entry) + "\n"
	#		out_test.write(str(line))
	#sys.exit(0)
	return (gene_nuc_dict, gene_dict, end_ignore)

def output_pooled(filelist):
	linelist = {}; info_dict = {}
	pooled = []
	ignore_genes = False
	pattern = r"^(A[t|l])Chr(\d{1,2})\t(\d*)\t(\w{1})\t(\d*)\t(.*)\t(.*)\t(.*)$"
	recomp = re.compile(pattern)
	#for filename in filelist:
	#	with open(filename) as file:
	#		linelist[filename] = file.readlines()
	for filename in filelist:
		print("Reading lines from", filename, "\n")
		#for line in linelist[filename]:
		with open(filename) as file:
			for line in file:
				match = recomp.match(line)
				if match:
					species = str(match.group(1)); chr_num = int(match.group(2)); pos = int(match.group(3)); con_nuc = str(match.group(4)); reads = int(match.group(5)); nucs = str(match.group(6)).upper(); qual1 = str(match.group(7)); qual2 = str(match.group(8))
					if (species, chr_num, pos) not in info_dict:
						info_dict[species, chr_num, pos] = [con_nuc, reads, nucs, qual1, qual2]
					else:
						info_dict[species, chr_num,pos][1] += int(reads)
						info_dict[species, chr_num,pos][2] += str(nucs)
						info_dict[species, chr_num,pos][3] += str(qual1)
						info_dict[species, chr_num,pos][4] += str(qual2)

	match = re.search(r"(Thaliana|Thalyrata).*(Ler|Sue|Care)",filelist[0])
	outfilename = ""
	if match:
		align, species = match.groups()
		if species == "Ler":
			outfilename = "Data/SynNatLibraries/ParentalData/Ler/PooledLer_" + str(align) + "Alignment.pileup"
		elif species == "Sue":
			outfilename = "Data/SueAll/Sue-OurRESCAN/PooledSue_" + str(align) + "Alignment.pileup"
		elif species == "Care":
			outfilename = "Data/CareAll/CareOurRESCAN/PooledCare_" + str(align) + "Alignment.pileup"
	with open(outfilename, 'w') as outfile:
		print("Sorting list and outputting pooled pileup\n\n")
		for entry in sorted(info_dict.items(), key = lambda entry: (entry[0][0], entry[0][1], entry[0][2])):
			(species, chr_num, pos), [con_nuc, reads, nucs, qual1, qual2] = entry
			line = str(species) + "Chr" + str(chr_num) + "\t" + str(pos) + "\t" + str(con_nuc) + "\t" + str(reads) + "\t" + str(nucs).upper() + "\t" + str(qual1) + "\t" + str(qual2) + "\n"
			outfile.write(line)

def pool_pileups(filelist, gene_dict, gene_nuc_dict, end_ignore):
	#AtChr5	11729122	C	24	.GGG..G.GGG...G....GGG.G	HGIIIIIIIIIIHHGIB>IIHEGI	F!!!FF!F!!!!FF!FF3F!!!'!
	pattern = r"^(A[t|l])Chr(\d{1,2})\t(\d*)\t(\w{1})\t(\d*)\t(.*)\t.*\t.*$"
	recomp = re.compile(pattern)
	linelist = {}; info_dict = {}
	pooled = []
	ignore_genes = True
	for filename in filelist:
		with open(filename) as file:
			linelist[filename] = file.readlines()
	for filename in filelist:
		for line in linelist[filename]:
			match = recomp.match(line)
			if match:
				species = str(match.group(1)); chr_num = int(match.group(2)); pos = int(match.group(3)); con_nuc = str(match.group(4)); reads = int(match.group(5)); nucs = str(match.group(6)).upper()
				if len(str(chr_num)) == 2: chr_num = int(str(chr_num)[-1])
				if ignore_genes == False and (species == "At" and chr_num > 3) or (species == "Al" and chr_num > 5):
					ignore_genes = True
				elif (species == "At" and chr_num == 3 and pos >= end_ignore) or (species == "Al" and chr_num == 5) and pos >= 20000000:
					ignore_genes = False
				
				if ignore_genes == False:
					if (species, chr_num, pos) not in gene_nuc_dict: continue # skip position b/c not within a gene
					if (species, chr_num, pos) in info_dict: # Dictionary for information to be pooled together
						info_dict[species, chr_num,pos][1] += int(reads)
						info_dict[species, chr_num,pos][2] += str(nucs)
					else:
						info_dict[species, chr_num, pos] = [con_nuc, reads, nucs]
	return info_dict

def reads_and_percentile_calc(genic_reads_filelist, gene_reads, end_ignore, species, percentile_75):
	genic_dict = {}
	pattern4 = r'^(AT[1-5]G\d{5}|scaffold_\d{1}_ID_\d{6}_(?:n|AT[1-5]G\d{5}))\t\d*\t(\d{1})\t(\d*)\t(\d*)\t.*$'
	recomp4 = re.compile(pattern4)
	pattern5 = r'(AT[1-5]G\d{5})'
	recomp5 = re.compile(pattern5)
	for file in genic_reads_filelist:
		with open(file) as infile:
			linelist = infile.readlines()
			for line in linelist:
				match = recomp4.match(line)
				if match:
					gene_name = str(match.group(1)); chr_num = int(match.group(2)); start_pos = int(match.group(3)); counts = int(match.group(4))
					if len(gene_name) > 9:
						match = recomp5.search(gene_name)
						if match:
							thal_gene_name = str(match.group(1))
							if chr_num == 5 and start_pos >= 20840000:
								if (species, thal_gene_name) in gene_reads:
									gene_reads[species, thal_gene_name] += 1
								else:
									gene_reads[species, thal_gene_name] = 1
					else:
						if chr_num == 3 and start_pos >= end_ignore:
							if (species, gene_name) in gene_reads:
								gene_reads[species, gene_name] += 1
							else:
								gene_reads[species, gene_name] = 1
					if gene_name in genic_dict:
						genic_dict[gene_name] += counts
					else:
						genic_dict[gene_name] = counts
	
	reads_array = sorted((reads for reads in genic_dict.values()))
	with open("Temp.txt", 'w') as outfile:#
		for entry in sorted(genic_dict):#
			outline = str(entry) + "\n"#
			outfile.write(outline)#
	sys.exit(0)#
	percentile_75[species] = scoreatpercentile(reads_array,75, (1,max(reads_array)))
	return(gene_reads, percentile_75)

def SNP_count(nucs, recomp2, recomp3):
	SNPs = {x: 0 for x in ('A','T','C','G')}
	maxSNP = ""
	maxSNP_value = 0
	countahead = 0
	indel = ""
	consensus_num = len(recomp3.findall(nucs))
	for m in recomp2.finditer(nucs):
		prev = str(m.group(1))
		SNP = str(m.group(2))
		if prev != "None": # nucleotide belongs to an indel
			countahead = int(prev[1]) - 1
			indel = str(prev) + str(SNP)
		elif countahead > 0: # continue adding nucleotides 1 at a time to indel
			indel += str(SNP)
			countahead -= 1
			if countahead == 0:
				if indel in SNPs:
					SNPs[indel] += 1
				else:
					SNPs[indel] = 1
		else: # SNP not part of an indel
			SNPs[SNP] += 1
	SNPcount = sum([v for v in SNPs.values()])
	maxSNP = max([k for k, v in SNPs.items()],key=lambda k: SNPs[k])
	maxSNP_value = SNPs[maxSNP]
	totalusedreads = consensus_num + maxSNP_value
	totalreads = consensus_num + SNPcount
	skip = 0
	if totalusedreads < 12: # skip position b/c < 12 reads to be used in Fisher's Exact Test. Example: ...^F........-2AT  Number of reads = 12 (11 consensus, 1 indel)
		skip = 1
	elif SNPcount == 0: # skip position b/c no SNPs are present (note: when this is looking at A. thaliana and not A. suecica, script will not skip position but instead use for SUN evidence)
		skip = 2
	return (consensus_num, totalreads, maxSNP, maxSNP_value, skip)

start = time.time()

gene_nuc_dict, gene_dict, end_ignore = get_genes()

#-1A, +2C, or just a nucleotide without a preceding [+-]\d{1}
pattern2 = r"([+-]\d{1})?([ATCG]{1})"
recomp2 = re.compile(pattern2)
pattern3 = r"[,.]"
recomp3 = re.compile(pattern3)
percentile_75 = {}

gene_reads = {}
species_list = ["As","At","Aa"]
for entry in gene_dict.keys():
	if len(entry) == 2:
		(species, gene_name) = entry
		for species in species_list:
			gene_reads[species, gene_name] = 0

# Set filelist to A. suecica pileups to pool them
#Following two lines of code used to pool together Suecica pileups into "PooledSue_Thaliana/ThalyrataAlignment.pileup"
#filelist = ["Data/SueAll/Sue-OurRESCAN/lib1Sue/ThalyrataAlignment/Pileup/libSue1F_pileup.txt","Data/SueAll/Sue-OurRESCAN/lib1Sue/ThalyrataAlignment/Pileup/libSue1R_pileup.txt","Data/SueAll/Sue-OurRESCAN/lib2Sue/ThalyrataAlignment/Pileup/libSue1_pileup.txt"]
#output_pooled(filelist)
filelist = ["Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/RESCANlibsSue_sxn_pileup.txt"]
sue_pool_dict = pool_pileups(filelist, gene_dict, gene_nuc_dict, end_ignore)
genic_reads_filelist = ["Data/Tilman_counts/As/sue_set1_gene_count.txt","Data/Tilman_counts/As/sue_set2.5_gene_count.txt","Data/Tilman_counts/As/sue_set2_F_side_gene_count.txt","Data/Tilman_counts/As/sue_set4_F_side_gene_count.txt"]
gene_reads, percentile_75 = reads_and_percentile_calc(genic_reads_filelist, gene_reads, end_ignore, "As", percentile_75)

# Set filelist to A. thaliana pileups to pool them
#Following two lines of code used to pool together Ler pileups into "PooledLer_Thaliana/ThalyrataAlignment.pileup"
#filelist = ["Data/SynNatLibraries/ParentalData/Ler/lib1Ler/ThalyrataAlignment/Pileup/libLer1F_pileup.txt","Data/SynNatLibraries/ParentalData/Ler/lib1Ler/ThalyrataAlignment/Pileup/libLer1R_pileup.txt","Data/SynNatLibraries/ParentalData/Ler/lib2Ler/ThalyrataAlignment/Pileup/libLer2_pileup.txt"]
#output_pooled(filelist)
filelist = ["Data/SynNatLibraries/ParentalData/Ler/PooledLer_ThalyrataAlignment.pileup"]
ler_pool_dict = pool_pileups(filelist, gene_dict, gene_nuc_dict, end_ignore)
genic_reads_filelist = ["Data/Tilman_counts/At/tsu-1/tsu-1_SRR013335_gene_count.txt", "Data/Tilman_counts/At/tsu-1/tsu-1_SRR013337_gene_count.txt", "Data/Tilman_counts/At/col-0/col-0_SRR013327_gene_count.txt","Data/Tilman_counts/At/col-0/col-0_SRR013328_gene_count.txt","Data/Tilman_counts/At/bur-0/bur-0_SRR013331_gene_count.txt","Data/Tilman_counts/At/bur-0/bur-0_SRR013333_gene_count.txt"]
gene_reads, percentile_75 = reads_and_percentile_calc(genic_reads_filelist, gene_reads, end_ignore, "At", percentile_75)

# Set filelist to A. arenosa pileups to pool them
#Following two lines of code used to pool together Ler pileups into "PooledLer_Thaliana/ThalyrataAlignment.pileup"
#filelist = ["Data/CareAll/CareOurRESCAN/lib2Care1/ThalyrataAlignment/Pileup/lib2Care1_pileup.txt","Data/CareAll/CareOurRESCAN/lib1Care1/ThalyrataAlignment/Pileup/libCare1F_pileup.txt", "Data/CareAll/CareOurRESCAN/lib1Care1/ThalyrataAlignment/Pileup/libCare1R_pileup.txt"]
#output_pooled(filelist)
filelist = ["Data/CareAll/CareOurRESCAN/PooledCare_ThalyrataAlignment.pileup"]
care_pool_dict = pool_pileups(filelist, gene_dict, gene_nuc_dict, end_ignore)
genic_reads_filelist = ["Data/Tilman_counts/Aa/black_forest/black_forest_1_F_side_gene_count.txt","Data/Tilman_counts/Aa/black_forest/black_forest_2_F_side_gene_count.txt","Data/Tilman_counts/Aa/borbasii/borbasii_1_F_side_gene_count.txt","Data/Tilman_counts/Aa/borbasii/borbasii_2_F_side_gene_count.txt","Data/Tilman_counts/Aa/care/care_PE_F_gene_count.txt","Data/Tilman_counts/Aa/care/care_SE_gene_count.txt","Data/Tilman_counts/Aa/strecno/sd_strecno_gene_count.txt"]
gene_reads, percentile_75 = reads_and_percentile_calc(genic_reads_filelist, gene_reads, end_ignore, "Aa", percentile_75)

#sys.exit(0)  # Used for when the four lines of code pooling individual Sue and Ler pileups are uncommented

sue_sun_filename = "Output/IsSUN_Sue.txt" # Debugging file
sue_notsun_filename = "Output/NotSUN_Sue.txt" # Debugging file
ler_sun_filename = "Output/IsSUN_Ler.txt" # Debugging file
ler_notsun_filename = "Output/NotSUN_Ler.txt" # Debugging file
care_sun_filename = "Output/IsSUN_Care.txt" # Debugging file
care_notsun_filename = "Output/NotSUN_Care.txt" # Debugging file

# SUN-finding algorithm that runs when finding SUNs in As. Looks for 50/50 Con/SNP in As, 100/0 in At, and ~2X reads in As compared to At
with open(sue_sun_filename, 'w') as sue_sun_file: # Debugging file
	with open(sue_notsun_filename, 'w') as sue_notsun_file: # Debugging file
		with open(ler_sun_filename, 'w') as ler_sun_file: # Debugging file
			with open(ler_notsun_filename, 'w') as ler_notsun_file: # Debugging file
				with open(care_sun_filename, 'w') as care_sun_file: # Debugging file
					with open(care_notsun_filename, 'w') as care_notsun_file: # Debugging file
						for entry in sorted(sue_pool_dict.items(), key = lambda entry: (entry[0][0], entry[0][1], entry[0][2])):
							(species, chr_num, pos), [con_nuc, reads, nucs] = entry
							consensus_num, totalreads, maxSNP, maxSNP_value, skip = SNP_count(nucs, recomp2, recomp3)
							if skip == 1 or skip == 2: continue # Either the total usable reads (consensus + maxSNP)  is < 12 (1), or no SNPs are present (2)
							is_sun = False
							#				Actual		Ideal (for 50/50)
							#	Consensus	  a			  b
							#		SNPs	  c			  b
							oddsratio, pvalue = fisher_exact([[consensus_num, totalreads*0.5], [maxSNP_value, totalreads*0.5]])
							if pvalue < 0: pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
							fisher_line = "Consensus: " + str(consensus_num) + "\t" + str(maxSNP) + ": " + str(maxSNP_value) + "\tP-value=" + str(pvalue) + "\n" # Debugging file
							if pvalue > 0.10:
								# To a 90% confidence level, cannot reject the hypothesis that the observed SNP frequency is the same as 50/50 consensus/SNP
								is_sun = True
							
							if is_sun == True:
								# 1a) If we believe position is a SUN, must check the SNP distribution for parent species to ensure that observed candidate SUN is not inherited heterozygosity
								if (species == "At" and (species, chr_num, pos) in ler_pool_dict) or (species == "Aa" and (species, chr_num, pos) in care_pool_dict):
									if species == "At":
										[parent_con_nuc, parent_reads, parent_nucs] = ler_pool_dict[(species, chr_num, pos)]
									else:
										[parent_con_nuc, parent_reads, parent_nucs] = care_pool_dict[(species, chr_num, pos)]
									parent_consensus_num, parent_totalreads, parent_maxSNP, parent_maxSNP_value, parent_skip = SNP_count(parent_nucs, recomp2, recomp3)
									if parent_skip == 1: continue # Total usable reads (consensus + maxSNP) in At/Aa is < 12, so skip
									
									#				Actual		Ideal (for 100/0)
									#	Consensus	  a			  b
									#		SNPs	  c			  d
									#print ("parent_consensus_num: ", str(parent_consensus_num)," parent_consensus_ideal: ", str(parent_totalreads), " parent_maxSNP_value: ", str(parent_maxSNP_value),sep="") #
									parent_pvalue = 1
									if parent_totalreads != parent_consensus_num: # If equal, no SNPs at all, so p-value=1. Don't bother with Fisher Exact Test (time-saving feature)
										oddsratio, parent_pvalue = fisher_exact([[parent_consensus_num, parent_totalreads], [parent_maxSNP_value, 0]])
									if parent_pvalue < 0: parent_pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
									
									#AtChr5	11729122	C	24	.GGG..G.GGG...G....GGG.G	HGIIIIIIIIIIHHGIB>IIHEGI	F!!!FF!F!!!!FF!FF3F!!!'!
									parent_line = "AtChr" + str(chr_num) + "\t" + str(pos) + "\t" + str(con_nuc) + "\t" + str(parent_reads) + "\t" + str(parent_nucs).upper() + "\n"
									parent_fisher_line = "Consensus: " + str(parent_consensus_num) + "\t" + str(parent_maxSNP) + ": " + str(parent_maxSNP_value) + "\tP-value=" + str(parent_pvalue) + "\n" # Debugging file
									if parent_pvalue > 0.10:
										# To a 90% confidence level, cannot reject the hypothesis that the observed SNP frequency is the same as 100/0 consensus/SNP
										if species == "At":
											ler_sun_file.write(parent_line)#Debugging file
											ler_sun_file.write(parent_fisher_line)#Debugging file
										else:
											care_sun_file.write(parent_line)#Debugging file
											care_sun_file.write(parent_fisher_line)#Debugging file
									else:
										# To a 90% confidence level, can reject the hypothesis that the observed SNP frequency is the same as 100/0 consensus/SNP
										is_sun = False
										if species == "At":
											ler_notsun_file.write(parent_line)#Debugging file
											ler_notsun_file.write(parent_fisher_line)#Debugging file
										else:
											care_notsun_file.write(parent_line)#Debugging file
											care_notsun_file.write(parent_fisher_line)#Debugging file
								else:
									# Nucleotide position is not within a gene
									is_sun = False
							
							#AtChr5	11729122	C	24	.GGG..G.GGG...G....GGG.G	HGIIIIIIIIIIHHGIB>IIHEGI	F!!!FF!F!!!!FF!FF3F!!!'!
							line = "AtChr" + str(chr_num) + "\t" + str(pos) + "\t" + str(con_nuc) + "\t" + str(reads) + "\t" + str(nucs).upper() + "\n" # Debugging file
							if is_sun == False:
								sue_notsun_file.write(line)#Debugging file
								sue_notsun_file.write(fisher_line)#Debugging file
							else:
								# 2) Testing for ~2X read count in A. suecica when compared to A. arenosa
								#    Normalization of read counts for a given species is accomplished by dividing a gene's read counts by the 75th percentile genic read count for that species
								#    As is being normalized by At: Divide As genic read counts by As's 75th percentile genic read count, then multiply As genic read counts by At's 75th percentile genic read count
								gene_name = gene_nuc_dict[(species, chr_num, pos)]
								sue_totalreads_norm = gene_reads["As", gene_name] / percentile_75["As"]
								parent_totalreads = 0
								if species == "At":
									sue_totalreads_norm *= percentile_75["At"]
									parent_totalreads = gene_reads["At", gene_name]
								else:
									sue_totalreads_norm *= percentile_75["Aa"]
									parent_totalreads = gene_reads["Aa", gene_name]
								#						Actual		Ideal (for 2X:1X)
								#	A. suecica			  a			  2b
								#	A. thaliana/arenosa	  b			  b
								oddsratio, pvalue_2X = fisher_exact([[sue_totalreads_norm, parent_totalreads*2], [parent_totalreads, parent_totalreads]])
								extra_line = "Sue norm reads: " + str(sue_totalreads_norm) + "\tAre norm reads: " + str(parent_totalreads) + "\tP-value=" + str(pvalue_2X) + "\n\n" # Debugging file
								fisher_line += extra_line
								if pvalue_2X > 0.10:
									gene_dict[species, gene_name][3] += 1
									sue_sun_file.write(line) # Debugging file
									sue_sun_file.write(fisher_line) #Debugging file
									care_sun_file.write(extra_line) # Debugging file
									continue
								sue_notsun_file.write(line) # Debugging file
								sue_notsun_file.write(fisher_line) # Debugging file
								care_notsun_file.write(extra_line) # Debugging file

# 3) Outputting a descending order list of top per-gene SUN counts in At-mapping genes in A. suecica
outfilename = "Output/Sue_SUN_Genes.txt"
gene_list = sorted(gene_dict.items(), key = lambda x: (x[1][3], x[0][0], x[1][0], x[1][1]), reverse=True)
with open(outfilename, 'w') as outfile:
	line = "Gene\tChr\tStart_pos\tEnd_pos\tSUN_Count\n"
	outfile.write(line)
	for entry in gene_list:
		try:
			(species, gene_name), [chr_num, start_pos, end_pos, SUN_count] = entry
			line = str(species) + ", " + str(gene_name) + "\t" + str(chr_num) + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(SUN_count) + "\n"
		except:
			(species, gene_name, i), [chr_num, start_pos, end_pos, SUN_count] = entry
			line = str(species) + ", " + str(gene_name) + "\t" + str(chr_num) + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(SUN_count) + "\n"
		outfile.write(line)

total = round((time.time() - start) / 60,2)
print("\nSUNs per gene counted (" + str(total) + " minutes)")