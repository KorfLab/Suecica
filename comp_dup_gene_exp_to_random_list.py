#!/usr/bin/env python3
import argparse, random
from collections import Counter
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF, SAM

def command_line():
	parser = argparse.ArgumentParser(description="comp_dup_gene_exp_to_random_list.py establishes whether or not a list of genes believed to be duplicated, on average, have higher expression in Thalyrata-mapped RNA-Seq data when compared to an equally sized list of randomly selected genes.")
	#parser.add_argument('-i', default="Data/RNA-Seq/AsLeaf_vs_ColLeaf/mRNA_Expression_0_hr_Technicals_Combined.tsv", type=str, help='Input file containing mRNA expression data', metavar='ExpressionFile')
	parser.add_argument('-i', default="Data/RNA-Seq/AsLeaf_vs_ColLeaf/mRNA_Expression_0_hr_Technicals_Combined_Sue_Col_Ratio.tsv", type=str, help='Input file containing mRNA expression data', metavar='ExpressionFile')
	parser.add_argument('-c', default=7, type=int, help='Column number of gene expression count table that contains the A. suecica data to be tested.', metavar='ColumnNumber')
	parser.add_argument('-sam', default="Data/RNA-Seq/sue1_leaf_mRNA-Seq/set3/sue_mRNA_set3_sw_aln.sam.gz", type=str, help='SAM file which read counts originated from. This will be read to get the total number of reads within the library for RPKM calculations.', metavar='')
	parser.add_argument('-d', default="Data/Duplicated_Genes/2500bp_window/Sue1Set5/As_Duplications.txt", type=str, help='Input file containing the list of genes believed to be duplicated in A. suecica', metavar='DuplicatedGenesFile')
	parser.add_argument('-gff', default="../Genomes/Thalyrata.gff", type=str, help='GFF file. To be used in RPKM calculation', metavar='GFF_File')
	parser.add_argument('-sample', default="", type=str, help='Input file containing the list of genes which you wish you randomly sample from. By default, if this parameter is not included, all genes found within the GFF file will be sampled from.', metavar='SamplingGeneList')
	parser.add_argument('-l', default=100000, type=int, help='Number of times the script should loop over randomly generating gene lists and comparing this sum to the duplicated gene list\'s sum', metavar='LoopCount')
	
	args = parser.parse_args()
	exp_file = args.i
	column_num = args.c
	sam_file = args.sam
	dup_file = args.d
	gff_file = args.gff
	sample_file = args.sample
	loop_count = args.l
	
	return(exp_file, column_num, sam_file, dup_file, gff_file, sample_file, loop_count)

exp_file, column_num, sam_file, dup_file, gff_file, sample_file, loop_count = command_line()

# Open GFF file and make dictionary containing gene names & their positions
gff_genes_dict = GFF.Parse(gff_file, create_nuc_dict=False)

# Store A. suecica list of DupHMM duplicated genes
dup_list = []
with open(dup_file) as infile:
	for gene_name in infile:
		gene_name = gene_name.strip()
		dup_list.append(gene_name)

# Get total number of reads w/ quality >= 20 for RPKM calculations
#parse_sam = SAM.Parse(sam_file)
#parse_sam.total_reads()
#parse_sam.start()
#total_reads = parse_sam.get_total_reads()

# Dictionary which stores whether or not a gene has <10 reads for any given sample
unusable_genes = {}

# Store A. suecica gene counts
gene_count_dict = Counter()
with open(exp_file) as infile:
	for line in infile:
		line = line.split()
		gene_name = line[0]
		if gene_name != "Gene": # Not looking at first line
			#gene_counts = int(line[column_num-1])
			gene_counts = float(line[column_num-1])
			
			# Raw read counts
			gene_count_dict[gene_name] = round(gene_counts,3)
			
			# Determine whether or not a gene has <10 reads for all combined technical replicate expression values
			if ((int(line[1]) + int(line[2])) < 10) or ((int(line[3]) + int(line[4])) < 10):
				unusable_genes[gene_name] = False
			
			# RPKM read counts
			#(chromosome, s_pos, e_pos) = gff_genes_dict.loc(gene_name)
			#dist_kb = (e_pos - s_pos) / 1000
			#RPK = gene_counts / dist_kb
			#RPKM = RPK / (total_reads / 1000000)
			#gene_count_dict[gene_name] = round(RPKM,3)

# If sampling file was provided, read in list of genes to support (in addition to dups)
if sample_file != "":
	sample_list = []
	with open(sample_file) as infile:
		for gene_name in infile:
			gene_name = gene_name.strip()
			sample_list.append(gene_name)
	for gene_name in dup_list:
		if gene_name not in sample_list:
			sample_list.append(gene_name)
	gene_count_dict = {key: value for key, value in gene_count_dict.items() if key in sample_list}

# Only make dup vs. non-dup comparison possible for genes w/ >=10 counts in both Col and Sue combined technical replicates
gene_count_dict = {gene: c for gene, c in gene_count_dict.items() if gene not in unusable_genes.keys()}
dup_list = [gene for gene in dup_list if gene not in unusable_genes.keys()]

# Sum up A. suecica duplicated gene counts
sue_dup_sum = sum([c for gene, c in gene_count_dict.items() if gene in dup_list])

# Create sampling gene dictionary that does not contain duplicated genes
sample_dict = {gene: c for gene, c in gene_count_dict.items() if gene not in dup_list}

# Generate a random list of genes for loop_count times, sum the counts for this random list, and compare
# this sum to the duplicated gene's list sum. Calculate the number of times out of loop_count randomly
# generated lists the duplicated gene list was greater than the randomly generated list.
sue_greater = 0
dup_list_len = len(dup_list)
for i in range(0,loop_count):
	rand_list = random.sample(sample_dict.keys(), dup_list_len)
	rand_list_sum = sum([c for gene, c in sample_dict.items() if gene in rand_list])
	if rand_list_sum < sue_dup_sum:
		sue_greater += 1

pct_greater = round(sue_greater / loop_count * 100, 1)
print("The duplicated list of genes' RPKM value was greater than an equally sized, random list of genes for "
	  + str(sue_greater) + " out of " + str(loop_count) + " randomly generated lists. (", str(pct_greater), "%)", sep='')