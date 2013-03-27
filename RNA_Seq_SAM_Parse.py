#!/usr/bin/env python3
import argparse, gzip, sys
from collections import Counter
sys.path.append("../Packages")
from GenomeTools import SAM, GFF

def command_line():
	parser = argparse.ArgumentParser(description="RNA_Seq_SAM_Parse.py parses mapped reads from a SAM file, and given gene locations from a GFF file, will output the number of reads mapped to each gene.")
	#parser.add_argument('-i', default="Data/RNA-Seq/sue_mRNA-Seq_leaf/SE_processing/sue_mRNA_leaf_sw_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="Data/RNA-Seq/sue1-mRNA-Seq/set3/Thalyrata/sue_mRNA_set3_sw_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="Data/RNA-Seq/Col_RNA-Seq/Thalyrata/SRR493036_aln.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_01nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_02nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_03nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_04nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_05nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_06nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_07nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_08nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_09nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	parser.add_argument('-i', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/SAM/libWhan_10nc.sam.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	#parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file')
	parser.add_argument('-gff', default="../TAIR9_GFF3_genes_transposons.gff", type=str, help='GFF file')
	#parser.add_argument('-o', default="Data/RNA-Seq/sue_mRNA-Seq_leaf/SE_processing/sue_mRNA_leaf_sw_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="Data/RNA-Seq/sue1-mRNA-Seq/set3/Thalyrata/sue_mRNA_set3_sw_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="Data/RNA-Seq/Col_RNA-Seq/Thalyrata/SRR493036_aln_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_01nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_02nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_03nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_04nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_05nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_06nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_07nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_08nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	#parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_09nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	parser.add_argument('-o', default="/home/mattchat/SuecicaDupSearch/Data/RNA-Seq/Col0_leaf_RNA-seq/libWhan_10nc_gene_count.tsv", type=str, help='Output file for gene-mapped read counts', metavar='OutputFile')
	
	args = parser.parse_args()
	i_file = args.i
	gff_file = args.gff
	o_file = args.o

	return(i_file, gff_file, o_file)

i_file, gff_file, o_file = command_line()
#gff_obj = GFF.Parse_Genes(gff_file, "AtChr1;AtChr2;AtChr3;AtChr4;AtChr5")
gff_obj = GFF.parse_genes(gff_file, "Chr1;Chr2;Chr3;Chr4;Chr5")
parse_sam = SAM.Parse(i_file)
parse_sam.reads_per_gene(gff_obj)
parse_sam.start()
gene_counts = parse_sam.get_reads_per_gene()
	
with open(o_file, 'w') as outfile:
	outline = "Gene\tReads_mapped\n"
	outfile.write(outline)
	for gene, count in sorted(gene_counts.items(), key = lambda gene_counts: gene_counts[1], reverse=True):
		outline = str(gene) + "\t" + str(count) + "\n"
		outfile.write(outline)