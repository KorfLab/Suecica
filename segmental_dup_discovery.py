#!/usr/bin/env python3
import argparse
from collections import Counter
import sys; sys.path.append('../Packages/')
from GenomeTools import GFF

def command_line():
	parser = argparse.ArgumentParser(description="SegmentalDupDiscovery.py, given a list of name-sorted duplicated genes, determines \
									 how many of them are present in series. Whether or not a gene is present in series is \
									 determined by a user-specified maximum distance within which two genes can be from one \
									 another and be considered part of the same segment. This distance value defaults to 7,500 bps.")
	parser.add_argument('-i', default="Data/DuplicatedGenes/2500bp_window/Sue1Set5.txt", type=str, help='File containing list of duplicated genes.', metavar='DupGenes')
	parser.add_argument('-d', default=7500, type=int, help='Maximum distance two genes can be apart to be called as part of the same segment.', metavar='SegmentDist')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file for the reference genome. Used to identify start and end positions for genes.', metavar='GFF_File')
	parser.add_argument('-o', default="Output/SegmentalDups/Sue1Set5SegDupHist.txt", type=str, help='Output filepath for expression comparison results', metavar='OutputPath')
	
	args = parser.parse_args()
	dups_file = args.i
	seg_dist = args.d
	gff_file = args.gff
	outpath = args.o
	
	return(dups_file, seg_dist, gff_file, outpath)

dups_file, seg_dist, gff_file, outpath = command_line()
dup_list = []
with open(dups_file) as infile:
	for line in infile:
		line = line.strip()
		dup_list.append(line)
gff_genes = GFF.Parse_genes(gff_file,create_nuc_dict=False)

segments = []
j = 0 # Creates j outside the scope of the following for loop
for i in range(0, len(dup_list)-1):
	# Compare current gene's location to its closest upstream neighbor in list of duplicated genes
	# If closest neighbor is present in segment, continue checking sequential upstream genes until
	# the full segment of genes is determined
	tmp_seg = []
	if i < j: continue # Skip genes which were found to be within a segmental duplication
	j = i
	end_loop = False
	while not end_loop:
		cur_gene = dup_list[j]
		if cur_gene not in tmp_seg: tmp_seg.append(cur_gene)
		cur_gene_loc = gff_genes.loc(cur_gene)
		cur_gene_chr, cur_gene_spos, cur_gene_epos = cur_gene_loc
		
		j += 1
		upstr_gene = dup_list[j]
		upstr_gene_loc = gff_genes.loc(dup_list[j])
		upstr_gene_chr, upstr_gene_spos, upstr_gene_epos = upstr_gene_loc
		if upstr_gene_chr == cur_gene_chr:
			# If upstream duplicated gene's first nucleotide is less than user-specific maximum distance from
			# the last nucleotide of our current duplicated gene, consider the upstream gene as part of a
			# segmental duplication and continue searching upstream for more segmentally duplicated genes
			if upstr_gene_spos - cur_gene_epos <= seg_dist:
				tmp_seg.append(upstr_gene)
			else:
				end_loop = True
		else:
			end_loop = True
		if j == len(dup_list) - 1:
			end_loop = True
	segments.append(tmp_seg)

# Output histogram of how many x-gene segments are present, as well as a list of duplicated genes
# newline-separated based on segments
seg_hist = Counter()
for entry in segments:
	seg_hist[len(entry)] += 1
with open(outpath, 'w') as outfile:
	outline = "Genes_in_segment\tSegment_count\n"
	outfile.write(outline)
	for entry in sorted(seg_hist.items()):
		outline = str(entry[0]) + "\t" + str(entry[1]) + "\n"
		outfile.write(outline)

outpath2 = outpath[:-4] + "_list" + outpath[-4:]
with open(outpath2, 'w') as outfile:
	for entry in segments:
		outline = '\n'.join(entry) + "\n" + "\n"
		outfile.write(outline)