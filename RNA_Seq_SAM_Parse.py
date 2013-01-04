#!/usr/bin/env python3
import argparse, gzip, sys
from collections import Counter
sys.path.append("../Packages")
from GenomeTools import SAM, GFF

def Command_line():
	parser = argparse.ArgumentParser(description="RNA_Seq_SAM_Parse.py parses mapped reads from a SAM file, and given gene locations from a GFF file, will output the number of reads mapped to each gene.")
	parser.add_argument('-i', default="Data/RNA-Seq/sue_mRNA-Seq_leaf/SE_processing/sue_mRNA_leaf_sw_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="Data/RNA-Seq/Col_RNA-Seq/Thalyrata/SRR493036_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="Data/RNA-Seq/sue1-mRNA-Seq/set3/Thalyrata/sue_mRNA_set3_sw_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file')
	parser.add_argument('-o', default="Data/RNA-Seq/sue_mRNA-Seq_leaf/SE_processing/sue_mRNA_leaf_sw_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="Data/RNA-Seq/Col_RNA-Seq/Thalyrata/SRR493036_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="Data/RNA-Seq/sue1-mRNA-Seq/set3/Thalyrata/sue_mRNA_set3_sw_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	
	args = parser.parse_args()
	i_file = args.i
	gff_file = args.gff
	o_file = args.o

	return(i_file, gff_file, o_file)

i_file, gff_file, o_file = Command_line()
gff_obj = GFF.Parse_Genes(gff_file, "AtChr1;AtChr2;AtChr3;AtChr4;AtChr5")
gene_counts = SAM.Reads_Per_Gene(i_file, gff_obj)
	
with open(o_file, 'w') as outfile:
	outline = "Gene\tReads_mapped\n"
	outfile.write(outline)
	for gene, count in sorted(gene_counts.items(), key = lambda gene_counts: gene_counts[1], reverse=True):
		outline = str(gene) + "\t" + str(count) + "\n"
		outfile.write(outline)