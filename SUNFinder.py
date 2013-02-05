#!/usr/bin/env python3.2
import argparse, re
from collections import Counter
from scipy.stats import fisher_exact
from random import sample
import sys; sys.path.append('../Packages/')
from GenomeTools import Pileup, GFF

# This script aids in the search for genes in A. suecica which underwent a single ancestral duplication event (shortly after the A. suecica
# species arose), which was then soon followed by at least one mutation in one of the duplicates. This script performs the following:
# 1) Grabs a list of duplicated genes and determines their start and end positions based on a GFF file. Optionally may search for SUNs +/-
#    some value relative to the start and end positions.
# 2) Finds nucleotide positions within/around supposedly duplicated regions that have a 50/50 reference/one-nucleotide SNP
#	 distribution, which is indicative of an ancestral duplication->mutation event.
#
# 3) After the script has run through every gene in duplicated regions and generated a list of SUN-containing genes, genes are output
#	 in descending order of total per-gene SUN counts.

def Command_line():
	parser = argparse.ArgumentParser(description="SUNFinder.py aids in the search for genes in A. suecica which underwent a single ancestral \
									 duplication event (shortly after the A. suecica species arose), which was then soon followed by at least \
									 one mutation in one of the duplicates. This script performs the following: \
									 (1) Grabs a list of duplicated genes and determines their start and end positions based on a GFF file. \
									 (2) Finds genic nucleotide positions within supposedly duplicated regions that have a 50/50 \
									 reference/one-nucleotide SNP distribution, which is indicative of an ancestral duplication->mutation event.\
									 (3) After the script has run through every gene in duplicated regions and generated a list of SUN-containing \
									 genes, genes are output in descending order of total per-gene SUN counts.")
	#parser.add_argument('-p', default="Data/Sue/sue1/single_bp_preprocess/set5/libSUE1_set5_aln_Aa_Filtered.pileup.gz", type=str, help='Pileup filename containing duplications with SNPs.', metavar='Pileup')
	parser.add_argument('-p', default="Data/Sue/SxN_F2s/HighPerformers/Al_filtered_SAM/HighPerformers.pileup.gz", type=str, help='Pileup filename containing duplications with SNPs.', metavar='Pileup')
	#parser.add_argument('-p', default="Data/Sue/SxN_F2s/LowPerformers/Al_filtered_SAM/LowPerformers.pileup.gz", type=str, help='Pileup filename containing duplications with SNPs.', metavar='Pileup')
	#parser.add_argument('-d', default="Data/As_Duplications_(Aa_Filtered).txt", type=str, help='File containing list of duplicated genes.', metavar='Duplications')
	parser.add_argument('-d', default="Data/As_Duplications_Meiosis_Genes_(Aa_Filtered).txt", type=str, help='File containing list of duplicated genes.', metavar='Duplications')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file for the reference genome. Used to identify start and end positions for genes (Default="../Thalyrata.gff").', metavar='GFF_File')
	parser.add_argument('-adj', default=0, help="Adjust search range for SUNs to +/- some value beyond/before the 3\'/5\' ends", metavar='+/- adjusted search space')
	parser.add_argument('-r', default='N', choices='YN', help='Set this option to \'Y\' to instead output the SUN count for a list of randomly selected genes equal in size to the list of duplicated genes.', metavar='RandomGenes')
	#parser.add_argument('-o', default="Output/Sue1Set5_AaFiltered_SUNGenes.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	#parser.add_argument('-o', default="Output/Sue1Set5_AaFiltered_SUNGenes_Random.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	#parser.add_argument('-o', default="Output/HighPerformers_AaFiltered_DupSUNGenes_67-33.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	#parser.add_argument('-o', default="Output/LowPerformers_AaFiltered_DupSUNGenes_67-33.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	parser.add_argument('-o', default="Output/HighPerformers_AaFiltered_AllGenes_50-50.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	#parser.add_argument('-o', default="Output/LowPerformers_AaFiltered_AllGenes_50-50.txt", type=str, help='Output filepath for SUN results', metavar='OutFilepath')
	
	args = parser.parse_args()
	pileup_file = args.p
	dup_file = args.d
	gff_file = args.gff
	adj_value = args.adj
	use_random = args.r
	sun_genes_filename = args.o
	
	return(pileup_file, dup_file, gff_file, adj_value, use_random, sun_genes_filename)
	
pileup_file, dup_file, gff_file, adj_value, use_random, sun_genes_filename = Command_line()

dup_list = []
with open(dup_file) as infile:
	for line in infile:
		line = line.strip()
		dup_list.append(line)
gff_genes_dict = GFF.Parse_Genes(gff_file)
if adj_value > 0:
	for gene_name, [chr, start_pos, end_pos] in gff_genes_dict.gene_dict.items():
		start_pos = start_pos - adj_value
		end_pos = end_pos + adj_value
		gff_genes_dict.gene_dict[gene_name] = [chr, start_pos, end_pos]
if use_random == "Y":
	# Construct a list of randomly selected genes to replace the list of duplicated genes
	#dup_list = sample(gff_genes_dict.gene_dict.keys(),len(dup_list))
	dup_list = gff_genes_dict.gene_dict.keys() # Delete me soon
dup_positions = ";".join([chr + ":" + str(s_pos) + "-" + str(e_pos) for gene_name, [chr,s_pos,e_pos] in sorted(gff_genes_dict.gene_dict.items(), key = lambda x: (x[1][0], x[1][1]) ) if gene_name in dup_list ] )
dup_pileup = Pileup.Parse(pileup_file, False, dup_positions)

suns_in_gene = Counter()
suns_in_gene_prob = {}

# SUN-finding algorithm that runs when finding SUNs in As. Looks for 50/50 Con/SNP in As
with open(sun_genes_filename[:-4] + "_SUNs.txt", 'w') as sue_sun_file: # Debugging file
	fisher_line = "Pos\tRef\tRefValue\tMaxSNP\tMaxSNPValue\tP-value\n"
	sue_sun_file.write(fisher_line)
	for entry in dup_pileup:
		(chr, pos), [ref, nuclist] = entry
		if ref not in ["A","C","G","T"]: continue # Skip ambiguous nucleotide positions
		gene_name = gff_genes_dict.get_gene_name(chr, pos)
		if gene_name in dup_list:
			totalreads, consensus_num, maxSNPs = Pileup.Fisher_SNP_Info(ref,nuclist)
			if totalreads < 12 or sum([int(v) for v in maxSNPs.values()]) == 0: continue # Either the total usable reads (consensus + maxSNP)  is < 12, or no SNPs are present
			is_sun = True
			#				Actual		Ideal (for 50/50)
			#	Consensus	  a			  b
			#		SNPs	  c			  d
			maxSNP, maxSNPvalue = maxSNPs.popitem(False)
			oddsratio, pvalue = fisher_exact([[consensus_num, totalreads*0.5], [maxSNPvalue, totalreads*0.5]])
			print(consensus_num,totalreads*0.5,sep='\t')
			print(maxSNPvalue,totalreads*0.5,sep='\t')
			print("\n")
			#oddsratio, pvalue = fisher_exact([[consensus_num, totalreads*(2/3)], [maxSNPvalue, totalreads*(1/3)]])
			if pvalue < 0: pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
			fisher_line = '\t'.join([str(chr) + ":" + str(pos),str(ref),str(consensus_num),str(maxSNP),str(maxSNPvalue),str(pvalue)]) + "\n"
			if pvalue < 0.10: is_sun = False # To a 90% confidence level, can reject the hypothesis that the observed SNP frequency is the same as 50/50 consensus/SNP
			
			if is_sun == True:
				if gene_name not in suns_in_gene:
					suns_in_gene[gene_name] = 1
					suns_in_gene_prob[gene_name] = [pvalue]
				else:
					suns_in_gene[gene_name] += 1
					suns_in_gene_prob[gene_name].append(pvalue)
				sue_sun_file.write(fisher_line)

gene_dict = {}
for gene_info in gff_genes_dict:
	gene_name, [chr, spos, epos] = gene_info
	if gene_name in suns_in_gene:
		avg_prob = sum(suns_in_gene_prob[gene_name]) / suns_in_gene[gene_name]
		gene_dict[gene_name] = [chr, spos, epos, suns_in_gene[gene_name], avg_prob]
	else:
		gene_dict[gene_name] = [chr, spos, epos, 0, 0]

# 3) Outputting a descending order list of top per-gene SUN counts for genes in the supposedly duplicated region of At-Chr.3 in A. suecica
with open(sun_genes_filename, 'w') as outfile:
	line = "Gene\tChr\tStart_pos\tEnd_pos\tSUN_Count\tSUN_Avg_Prob\n"
	outfile.write(line)
	for entry in sorted(gene_dict.items(), key = lambda x: (x[1][3], x[1][4]), reverse=True):
		gene_name, [chr_num, spos, epos, SUN_count, SUN_count_avg_prob] = entry
		if gene_name in dup_list:
			line = "\t".join([gene_name,chr_num,str(spos),str(epos),str(SUN_count), str(SUN_count_avg_prob)]) + "\n"
			outfile.write(line)