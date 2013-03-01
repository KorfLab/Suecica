#!/usr/bin/env python3
import argparse, random
from collections import Counter
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF, SAM

def command_line():
	parser = argparse.ArgumentParser(description="CompDupGeneExpToRandomList.py establishes whether or not a list of genes believed to be duplicated, on average, have higher expression in Thalyrata-mapped RNA-Seq data when compared to an equally sized list of randomly selected genes.")
	parser.add_argument('-i', default="Data/RNA-Seq/AsLeaf_vs_ColLeaf/mRNA_Expression_Technicals_Combined.tsv", type=str, help='Input file containing mRNA expression data for A. suecica', metavar='ExpressionFile')
	parser.add_argument('-c', default=4, type=int, help='Column number of gene expression count table that contains the A. suecica data to be tested.', metavar='ColumnNumber')
	parser.add_argument('-sam', default="Data/RNA-Seq/sue1_leaf_mRNA-Seq/set3/sue_mRNA_set3_sw_aln.sam.gz", type=str, help='SAM file which read counts originated from. This will be read simply to get an idea of the total number of reads within the library', metavar='')
	parser.add_argument('-d', default="Data/As_Duplications_(Aa_Filtered).txt", type=str, help='Input file containing the list of genes believed to be duplicated in A. suecica', metavar='DuplicatedGenesFile')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file. To be used in RPKM calculation', metavar='GFF_File')
	parser.add_argument('-l', default=10000, type=int, help='Number of times the script should loop over randomly generating gene lists and comparing this sum to the duplicated gene list\'s sum', metavar='LoopCount')
	
	args = parser.parse_args()
	exp_file = args.i
	column_num = args.c
	sam_file = args.sam
	dup_file = args.d
	gff_file = args.gff
	loop_count = args.l
	
	return(exp_file, column_num, sam_file, dup_file, gff_file, loop_count)

exp_file, column_num, sam_file, dup_file, gff_file, loop_count = command_line()

# Open GFF file and make dictionary containing gene names & their positions
gff_genes_dict = GFF.Parse_genes(gff_file, create_nuc_dict=True)

# Store A. suecica list of DupHMM duplicated genes
dup_list = []
with open(dup_file) as infile:
	for line in infile:
		line = line.strip()
		dup_list.append(line)
dup_list_len = len(dup_list)

# Get total number of reads w/ quality >= 20
total_reads = SAM.total_reads(sam_file)

# Store A. suecica gene counts
gene_count_dict = Counter()
with open(exp_file) as infile:
	for line in infile:
		line = line.split()
		gene_name = line[0]
		if gene_name != "Gene":
			gene_counts = int(line[column_num-1])
			(chromosome, s_pos, e_pos) = gff_genes_dict.loc(gene)
			dist_kb = (e_pos - s_pos) / 1000
			RPK = count / dist_kb
			RPKM = RPK / (total_reads / 1000000)
			gene_count_dict[gene_name] = gene_counts

# Sum up A. suecica duplicated gene counts
sue_dup_sum = sum([c for gene, c in gene_count_dict.items() if gene in dup_list])

# Generate a random list of genes for loop_count times, sum the counts for this random list, and compare
# this sum to the duplicated gene's list sum. Calculate the number of times out of loop_count randomly
# generated lists the duplicated gene list was NOT greater than the randomly generated list.
sue_not_greater = 0
for i in range(0,loop_count):
	rand_list = random.sample(gene_count_dict.keys(), dup_list_len)
	rand_list_sum = sum([c for gene, c in gene_count_dict.items() if gene in rand_list])
	if rand_list_sum > sue_dup_sum:
		sue_not_greater += 1

print("The duplicated list of genes was not greater than an equally sized, randomly generated list of genes' expression sum for " + str(sue_not_greater) + " out of " + str(loop_count) + " randomly generated lists.")