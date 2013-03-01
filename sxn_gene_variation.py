#!/usr/bin/env python3
import argparse, os, gzip
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF, SAM

def command_line():
	parser = argparse.ArgumentParser(description="sxn_dup_gene_variation.py totals the number of reads mapping to each gene for all provided SxN F2 SAM file libraries. It then outputs the mean and variance of each gene into a .tsv file format.")
	parser.add_argument('-i', default="Data/Sue/SxN_F2s/VarianceExp", type=str, help='Input folder containing all the SxN F2 SAM libraries.', metavar='Input_Folder')
	parser.add_argument('-o', default="Output/SxN_Variance.tsv", type=str, help='Tab-separated value table containing the mean and variance for each gene.', metavar='Output_File')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file containing start and end positions for each gene', metavar='GFF_File')
	
	args = parser.parse_args()
	in_folder = args.i
	o_file = args.o
	gff_file = args.gff
	
	return(in_folder, o_file, gff_file)

in_folder, o_file, gff_file = command_line()

# Open GFF file and make dictionary containing gene names & their positions
# Then add the per-library gene read counts to a dictionary
gff_genes_dict = GFF.Parse_genes(gff_file, create_nuc_dict=True)
library_gene_counts = {}
for (gene, [chromosome, spos, epos]) in gff_genes_dict:
	library_gene_counts[gene] = []
for root, subfolders, files in os.walk(in_folder):
	for f_name in files:
		print("Working on", str(f_name))
		f_in = os.path.join(root,f_name)
		temp_gene_counts = SAM.reads_per_gene(f_in, gff_genes_dict)
		
		# Count up the total number of reads in the library
		gzipped = False
		total_reads = 0
		if f_in.endswith(".gz"):
			infile = gzip.open(f_in, 'rb')
			gzipped = True
		else:
			infile = open(f_in)
		for line in infile:
			if gzipped == True: line = str(line, encoding='utf8')
			line = line.split()
			if line[0][:1] != "@" and int(line[4]) >= 20: total_reads += 1
		
		# Change read counts to RPKM value
		for gene, count in temp_gene_counts.items():
			(chromosome, s_pos, e_pos) = gff_genes_dict.loc(gene)
			dist_kb = (e_pos - s_pos) / 1000
			RPK = count / dist_kb
			RPKM = RPK / (total_reads / 1000000)
			library_gene_counts[gene].append(RPKM)

# Calculate gene read count means and variance
gene_means = {}
gene_variances = {}
for gene, RPKM_list in library_gene_counts.items():
	gene_mean = sum(RPKM_list) / len(RPKM_list)
	gene_means[gene] = gene_mean
	
	gene_var = sum( [(x-gene_mean)**2 for x in RPKM_list] ) / (len(RPKM_list) - 1)
	gene_variances[gene] = gene_var

with open(o_file, 'w') as outfile:
	outline = "Gene\tMean\tVariance\n"
	outfile.write(outline)
	for gene in sorted(gff_genes_dict.gene_dict.keys()):
		outline = "\t".join([str(gene), str(gene_means[gene]), str(gene_variances[gene])]) + "\n"
		outfile.write(outline)