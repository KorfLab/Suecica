#!/usr/bin/env python
import sys; sys.path.append('../Packages/')
from GenomeTools import Pileup, SAM
import argparse, gzip, re
import pdb # Delete me

def command_line():
	parser = argparse.ArgumentParser(description="Aa_Mismapping.py (1) helps identify positions in A. suecica which might contain \
									 mismapped A. arenosa reads, and (2) parses an A. suecica data file to remove reads which contribute \
									 mismapping-based SNPs. Option (1) is performed by outputting a list of SNPs found in a pileup file \
									 generated from an A. arenosa library mapped to TAIR9. Option (2) is performed by taking in the SNP \
									 list from Option 1 and removing reads from an A. suecica SAM file which contain such SNPs.")
	parser.add_argument('-p', default="Data/Care/2008_40bpPE/080801_SOLEXA1_FC3094DAAXX_combined_aln.pileup.gz", type=str, help='(For Option 1): Pileup file which was generated from A. arenosa reads mapping to TAIR9.', metavar='Pileup')
	parser.add_argument('-o1', default="Output/Aa_TAIR9_SNPs.txt.gz", type=str, help='(For Option 1): Output file that will contain a list of TAIR9 positions with A. arenosa mismapping SNPs.', metavar='Option1Output')
	parser.add_argument('-sam', default="Data/Sue/sue1/single_bp_preprocess/set5/libSUE1_set5_aln.sam.gz", type=str, help='(For Option 2): SAM file containing the A. suecica reads which will be filtered for A. arenosa mismapping SNPs.', metavar='SAM_File')
	parser.add_argument('-o2', default="Data/Sue/sue1/single_bp_preprocess/set5/libSUE1_set5_aln_Aa_Filtered.sam.gz", type=str, help='(For Option 2): Output file for A. arenosa-filtered SAM file', metavar='Option2Output')
	
	args = parser.parse_args()
	pileup_file = args.p
	o_1_file = args.o1
	sam_file = args.sam
	o_2_file = args.o2
	
	return(pileup_file, o_1_file, sam_file, o_2_file)

def arenosa_snps(pileup_file, o_1_file):
	Aa_SNPs = Pileup.parse(pileup_file, region_list="Chr1;Chr2;Chr3;Chr4;Chr5") # Only saves SNP positions by default
	with gzip.open(o_1_file, 'wb') as outfile:
		for entry in sorted(Aa_SNPs.pileup_dict.items(), key=lambda entry: (entry[0][0], entry[0][1])):
			(chr, pos), [ref, nuclist] = entry
			#nuclist_str = ""
			#for entry in nuclist:
			#	nuclist_str += str(entry) + "\t"
			#nuclist_str = nuclist_str[:-1]
			nuclist_str = '\t'.join([entry for entry in nuclist])
			outline = '\t'.join(["At" + str(chr),str(pos),str(ref),str(nuclist_str)]) + "\n"
			outfile.write(bytes(outline,"UTF-8"))
	return()
	
def filter_sam(sam_file, SNP_file, o_2_file):
	pattern = r"(\d+)([MSID])"
	recomp = re.compile(pattern)
	Aa_SNPs = Pileup.parse(SNP_file, simplified=True, only_save_SNPs=False, region_list="AtChr1;AtChr2;AtChr3;AtChr4;AtChr5")
	if sam_file[-3:] == ".gz":
		infile = gzip.open(sam_file, 'rb')
	else:
		infile = open(sam_file)
	with gzip.open(o_2_file, 'wb') as outfile:
		for line in infile:
			if sam_file[-3:] == ".gz": line = str(line, encoding='utf8')
			line_split = line.split("\t")
			if line_split[0][:1] != "@":
				chr = line_split[2]
				if chr == "*":
					outfile.write(bytes(line,"UTF-8"))
					continue
				skip_read = False
				for i in range(11,len(line_split)):
					if line_split[i] == "NM:i:0": # Read has a perfect sequence match
						outfile.write(bytes(line,"UTF-8"))
						skip_read = True
						break
				if skip_read == True: continue
				read_name = line_split[0] # Delete me
				s_pos = int(line_split[3])
				cigar = line_split[5]
				cigar_pieces = recomp.findall(cigar)
				nucs = line_split[9]
				current_pos = s_pos
				for piece in cigar_pieces:
					val, operation = piece
					val = int(val)
					if operation == "M": # A, C, T, or G
						for i in range(current_pos, current_pos + val): # Check all nucleotides for match w/ reference
							ref, nuclist = Aa_SNPs.PosInfo(chr, i)
							if ref != -1: # Ref == -1 when there isn't an A. arenosa SNP at position (chr, i)
								cur_nuc = nucs[i - s_pos:i - s_pos + 1]
								if cur_nuc != ref: # Read doesn't contain reference nucleotide
									total_nucs = Pileup.sum_nucs(nuclist, True)
									SNP = {}
									SNP["A"] = nuclist[0]
									SNP["C"] = nuclist[1]
									SNP["G"] = nuclist[2]
									SNP["T"] = nuclist[3]
									if cur_nuc in SNP and (SNP[cur_nuc] / total_nucs) >= 0.05:
										# SNP present in >=5% of A. arenosa reads mapped to TAIR9, so consider it a valid mismap; skip read
										print(str(cur_nuc), " (", str(ref), "): ", str(round(SNP[cur_nuc] / total_nucs,2)), "\t", str(chr), ":", str(i), "\t", str(read_name), sep='') # Delete me
										skip_read = True
										break
					else: # Insertion or deletion
						ref, nuclist = Aa_SNPs.PosInfo(chr, i)
						if ref != -1: # Ref == -1 when there isn't an A. arenosa SNP at position (chr, i)
							if operation == "I":
								indel_str = "+" + str(val)
							else:
								indel_str = "-" + str(val)
							for i in range(4,len(nuclist)):
								if nuclist[i][:2] == indel_str:
									total_nucs = Pileup.sum_nucs(nuclist, True)
									if ( int(nuclist[i][3+int(val):]) / total_nucs ) >= 0.05:
										# Indel present in >=5% of A. arenosa reads mapped to TAIR9, so consider it a valid mismap; skip read
										print(str(nuclist[i]), " (", str(ref), "): ", str(round(int(nuclist[i][3+int(val):]) / total_nucs,2)), str(chr), ":", str(i), "\t", str(read_name), sep='') # Delete me
										skip_read = True
										break
					if skip_read == True:
						break
					current_pos += val
				if skip_read == True:
						continue
			outfile.write(bytes(line,"UTF-8"))

pileup_file, o_1_file, sam_file, o_2_file = command_line()
if pileup_file != "" and o_1_file != "":
	arenosa_snps(pileup_file, o_1_file)
if o_1_file != "" and sam_file != "" and o_2_file != "":
	filter_sam(sam_file, o_1_file, o_2_file)