#!/usr/bin/env python3
from collections import Counter
import gzip
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF

#dup_filename = ""
#dup_filename = "Data/As_Duplications_Meiosis_Genes_(Aa_Filtered).txt"
dup_filename = "Data/As_Duplications_(Aa_Filtered).txt"
high_filename = "Data/Sue/SxN_F2s/HighPerformers/Al_filtered_SAM/HighPerformers.sam.gz"
low_filename = "Data/Sue/SxN_F2s/LowPerformers/Al_filtered_SAM/LowPerformers.sam.gz"
#outfilename = "Output/High_vs_Low_Meiosis_Reads.txt"
outfilename = "Output/High_vs_Low_AllDupGene_Reads.txt"
#outfilename = "Output/High_vs_Low_AllGene_Reads.txt" # For this: dup_filename = ""
#high_filename = "Data/Sue/SxN_F2s/RandomExp/HighPerformers/Al_filtered_SAM/HighPerformers.sam.gz"
#low_filename = "Data/Sue/SxN_F2s/RandomExp/LowPerformers/Al_filtered_SAM/LowPerformers.sam.gz"
#outfilename = "Output/High_vs_Low_RandomGroups_AllDupGene_Reads.txt"
gff_file = "../Thalyrata.gff"

high_total_cov = 0
low_total_cov = 0
high_gene_reads = Counter()
low_gene_reads = Counter()

# Open GFF file and make dictionary containing gene names & their positions
gff_genes_dict = GFF.Parse_genes(gff_file, create_nuc_dict=True)

# Open list of duplicated meiosis genes and store them
if dup_filename != "":
	dup_genes = []
	with open(dup_filename) as infile:
		for line in infile:
			line = line.strip()
			dup_genes.append(line)
	for entry in dup_genes:
		high_gene_reads[entry] = 1
		low_gene_reads[entry] = 1

# Examine all duplicated genes, not just meiosis duplicated genes
if dup_filename == "":
	dup_genes = gff_genes_dict.gene_dict.keys()
	for entry in gff_genes_dict.gene_dict.keys():
		high_gene_reads[entry] = 1
		low_gene_reads[entry] = 1

# Count up the number of reads per gene for the High Performers file
with gzip.open(high_filename, 'rb') as infile:
	for line in infile:
		line = str(line, encoding='utf8')
		if line[:1] != "@":
			line = line.split("\t")
			high_total_cov += 1
			chromosome = str(line[2])
			s_pos = int(line[3])
			gene_name = gff_genes_dict.get_gene_name(chromosome, s_pos)
			if gene_name != None:
				high_gene_reads[gene_name] += 1

# Count up the number of reads per gene for the Low Performers file
with gzip.open(low_filename, 'rb') as infile:
	for line in infile:
		line = str(line, encoding='utf8')
		if line[:1] != "@":
			line = line.split("\t")
			low_total_cov += 1
			chromosome = str(line[2])
			s_pos = int(line[3])
			gene_name = gff_genes_dict.get_gene_name(chromosome, s_pos)
			if gene_name != None:
				low_gene_reads[gene_name] += 1

# Normalize high performing individuals' read counts
high_dif = low_total_cov / high_total_cov
if high_total_cov > low_total_cov:
	print("High Performer group reads are scaled down by ", str(high_dif), sep='')
else:
	print("High Performer group reads are scaled up by ", str(high_dif), sep='')
high_gene_reads = {gene_name: x * high_dif for gene_name, x in high_gene_reads.items()}

# Output ratio of high performing vs. low performing reads
with open(outfilename, 'w') as outfile:
	outline = "Gene\tHighPerformanceNormalizedCount\tLowPerformanceNormalizedCount\tRatio\n"
	outfile.write(outline)
	for gene in sorted(dup_genes):
		outline = "\t".join([gene, str(high_gene_reads[gene]), str(low_gene_reads[gene]), str(round(high_gene_reads[gene]/low_gene_reads[gene],2)), "\n"])
		outfile.write(outline)