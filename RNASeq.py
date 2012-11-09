#!/usr/bin/env python3.2
import argparse
from collections import Counter
from random import sample
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF, SAM

def Command_line():
	parser = argparse.ArgumentParser(description="RNASeq.py compares the normalized expression levels per gene between two given input files. \
									 A list of duplicated genes of interest can be provided, and their overall variation in differential \
									 expression compared to variation in differential expression among the remaining genes.")
	parser.add_argument('-d', default="Data/RNA-Seq/sue1-mRNA-Seq/sue_set3_Matt/Thalyrata/sue_mRNA_set3_sw_aln.sam.gz", type=str, help='Experimental SAM file containing duplications.', metavar='ExpSAM')
	parser.add_argument('-c', default="Data/RNA-Seq/Col_RNA-Seq/Thalyrata/SRR493036_aln.sam.gz", type=str, help='Control SAM file.', metavar='ControlSAM')
	parser.add_argument('-l', default="Data/Suecica_Duplications_Aa_Filtered.txt", type=str, help='List of duplicated genes for differential expression analysis in experimental vs. control data sets.', metavar='GeneList')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file for the reference genome. Used to identify start and end positions for genes.', metavar='GFF_File')
	parser.add_argument('-s', default="N", choices='YN', help='Choose \'Y\' if you want the expression of your list of genes of interest to be compared to an equally sized random sample of non-interest genes. Leave this to its default \'N\' if you wish for the list to be compared to the remaining list of genes.', metavar='Sample')
	parser.add_argument('-o', default="Output/Sue1Set3_vs_Col_RNA-Seq.txt", type=str, help='Output filepath for expression comparison results', metavar='OutputPath')
	
	args = parser.parse_args()
	dup_SAM_file = args.d
	c_SAM_file = args.c
	gene_file = args.l
	gff_file = args.gff
	use_sample = args.s
	outpath = args.o
	
	return(dup_SAM_file, c_SAM_file, gene_file, gff_file, use_sample, outpath)

dup_SAM_file, c_SAM_file, gene_file, gff_file, use_sample, outpath = Command_line()
dup_list = []
with open(gene_file) as infile:
	for line in infile:
		line = line.strip()
		dup_list.append(line)
gff_genes_dict = GFF.Parse_Genes(gff_file)
comp_dict = Counter({key: 0 for key in gff_genes_dict.gene_dict.keys() if key not in dup_list})
if use_sample == "Y":
	# Construct a list of randomly selected genes to compare to the list of duplicated genes
	comp_list = sample(comp_dict.keys(),len(dup_list))
	comp_dict = Counter({gene: 0 for gene in comp_list})
dup_sam_counter = SAM.Reads_Per_Gene(dup_SAM_file, gff_genes_dict)
control_sam_counter = SAM.Reads_Per_Gene(c_SAM_file, gff_genes_dict)
diff_exp = {gene_name: round(exp_count / control_sam_counter[gene_name], 5) for gene_name, exp_count in dup_sam_counter.items()}
with open(outpath, 'w') as outfile:
	for entry in sorted(diff_exp):
		outline = str(gene_name) + "\n"
		outfile.write(outline)