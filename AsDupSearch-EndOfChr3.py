#!/usr/bin/env python3.2
from collections import Counter
from scipy.stats import fisher_exact
import sys; sys.path.append('/home/mattchat/Packages/')
from GenomeTools import Pileup, GFF

# This script searches for genes in A. suecica which underwent a single ancestral duplication event (shortly after the A. suecica
# species arose), which was then soon followed by a mutation in one of the duplicates. Specifically, this script narrows down on
# the end of At-Chr.3/Aa-Chr.5, where the end of At-Chr.3 rearranged onto the end of Aa-Chr.5
# This is performed by:
# 1)  The script looks at the supposedly duplicated region and finds the average read count value over this region, which is then compared to the average
#		 read count outside of this supposedly duplication region (as well as outside the end of At-Chr.3 which has undergone a translocation/duplication event).
#	     Should a Fisher's Exact test find the ratio of reads in the supposedly segmentally duplicated region to non-duplicated region to be 2:1, the program continues
#		 looking for SUNs that might identify unique positions.
#
# 2)  Looking for intragenic nucleotide positions within the supposedly duplicated region that have a 50/50 consensus/one-nucleotide SNP
#		 distribution, which is indicative of an ancestral duplication->mutation event.
#
# 3)  After the script has run through every gene in the duplicated region and generated a list of SUN-containing genes, these are output in descending
#		 order of total per-gene SUN counts.

def dup_read_check(sue_dup_pileup):
	# Checks for 2:1 dup region:non-dup region read count ratio
	sue_pileup = Pileup.Parse(sue_filename,True,["AlChr11","AlChr12","AlChr13","AlChr14","AlChr15","AlChr16","AlChr17","AlChr18"])
	#sue_pileup = Pileup.Parse(sue_filename,True,"AlChr15")
	non_dup_avg = sue_pileup.Coverage()
	rand = 1000
	dup_avg, rand_overlap = sue_pileup.Coverage("AlChr15:20840605-21219165", rand)
	
	## 1) If we believe the region to be segmentally duplicated, must check the read count level to see that it is approximately doubled compared to the
	##		 average read count level among genes not in the supposed duplicated region
	##
	##								Actual		Ideal (for 2:1)
	##	At-Chr.3 Dup Region			  a			  2b
	##	At-Chr.3 Non-dup Region		  b			  b
	oddsratio, pvalue = fisher_exact([[dup_avg, 2*non_dup_avg], [non_dup_avg, non_dup_avg]])
	if pvalue < 0: pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
	print("Duplicated region avg:\t",round(dup_avg,3))
	print("Non-dup region avg:\t",round(non_dup_avg,3))
	print("Probability of duplicated region having 2:1 ratio with non-dup region: ",round(pvalue,3))
	if pvalue < 0.10:  # To a 90% confidence level, can reject the hypothesis that the observed dup:non-dup read count ratio is not 2:1
		print("Supposedly duplicated region does not have a 2:1 read count ratio with the non-duplicated region of At-Chr.3")
		sys.exit()
	else:
		print("Supposedly duplicated region has a 2:1 read count ratio with the non-duplicated region of At-Chr.3")
		print("Number of times dup. region read coverage was less than the read coverage from ", rand, " lists of randomly selected positions' coverage: ", rand_overlap, sep="")

dup_genes = GFF.Parse_Genes("../thalyrata.gff","AlChr15:20840605-21219165")

sue_filename = "Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/cleaned-RESCANlibsSue_sxn_pileup_simplified.txt"
sue_dup_pileup = Pileup.Parse(sue_filename,True,"AlChr15:20840605-21219165")

#dup_read_check(sue_dup_pileup)

sun_genes = Counter()
 
sue_sun_filename = "Output/IsSUN_AaChr5.txt" # Debugging file
sue_notsun_filename = "Output/NotSUN_AaChr5.txt" # Debugging file

# SUN-finding algorithm that runs when finding SUNs in As. Looks for 50/50 Con/SNP in As
with open(sue_sun_filename, 'w') as sue_sun_file: # Debugging file
	with open(sue_notsun_filename, 'w') as sue_notsun_file: # Debugging file
		fisher_line = "Pos\tRef\tRefValue\tMaxSNP\tMaxSNPValue\tP-value\n"
		sue_sun_file.write(fisher_line)
		sue_notsun_file.write(fisher_line)
		for entry in sue_dup_pileup:
			(chr, pos), [ref, nuclist] = entry
			gene_name = dup_genes.get_gene_name(chr, pos)
			if gene_name in dup_genes:
				totalreads, consensus_num, maxSNPs = Pileup.Fisher_SNP_Info(ref,nuclist)
				if totalreads < 12 or sum([int(v) for v in maxSNPs.values()]) == 0: continue # Either the total usable reads (consensus + maxSNP)  is < 12, or no SNPs are present
				is_sun = True
				#				Actual		Ideal (for 50/50)
				#	Consensus	  a			  b
				#		SNPs	  c			  d
				maxSNP, maxSNPvalue = maxSNPs.popitem(False)
				oddsratio, pvalue = fisher_exact([[consensus_num, totalreads*0.5], [maxSNPvalue, totalreads*0.5]])
				if pvalue < 0: pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
				fisher_line = '\t'.join([str(chr) + ":" + str(pos),str(ref),str(consensus_num),str(maxSNP),str(maxSNPvalue),str(pvalue)]) + "\n"
				if pvalue < 0.10: is_sun = False # To a 90% confidence level, can reject the hypothesis that the observed SNP frequency is the same as 50/50 consensus/SNP
				
				if is_sun == True:
					if gene_name not in sun_genes:
						sun_genes[gene_name] = 1
					else:
						sun_genes[gene_name] += 1
					sue_sun_file.write(fisher_line) #Debugging file
				else:
					sue_notsun_file.write(fisher_line) # Debugging file

gene_dict = {}
for gene_info in dup_genes:
	gene_name, [chr, spos, epos] = gene_info
	if gene_name in sun_genes:
		gene_dict[gene_name] = [chr, spos, epos, sun_genes[gene_name]]
	else:
		gene_dict[gene_name] = [chr, spos, epos, 0]

# 3) Outputting a descending order list of top per-gene SUN counts for genes in the supposedly duplicated region of Aa-Chr.5 in A. suecica
outfilename = "Output/Sue_SUN_Genes_(Within_Aa-Chr5).txt"
line = "Gene\tChr\tStart_pos\tEnd_pos\tSUN_Count\n"
for entry in sorted(gene_dict.items(), key = lambda x: x[1][3], reverse=True):
	gene_name, [chr_num, spos, epos, SUN_count] = entry
	line += "\t".join([gene_name,chr_num,str(spos),str(epos),str(SUN_count)]) + "\n"
with open(outfilename, 'w') as outfile:
	outfile.write(line)