#!/usr/bin/env python3.2
import re, sys, argparse, time; from operator import itemgetter

def get_genes():
	filename = "Data/t9_sorted.gff"
	#Chr1	TAIR9	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
	pattern = r'^Chr([1-5])\t.*\t.*\t(\d*)\t(\d*)\t.*\t[+-.]\t.*\tID=(.*);Note=.*$'
	gene_nuc_dict = {}
	reads = {}
	with open(filename) as gene_file:
		linelist = gene_file.readlines()
		for line in linelist:
			match = re.search(pattern,line)
			if match:
				chr_num = str(match.group(1)); start_pos = int(match.group(2)); end_pos = int(match.group(3)); gene_name = str(match.group(4));
				reads_gene_info = [0, chr_num, start_pos]
				reads[gene_name] = reads_gene_info
				for i in range(start_pos, end_pos+1):
					gene_info = chr_num, i
					gene_nuc_dict[gene_info] = gene_name
	return (gene_nuc_dict, reads)

def options():
	parser = argparse.ArgumentParser(description='The -p and -f arguments specify the path and filename for the SAM file.')
	parser.add_argument('-p', metavar='Path', type=str, default="Data\SueAll\Sue-OurRESCAN\libsSue\ThalyrataAlignment", help='The path to the SAM file')
	parser.add_argument('-f', metavar='Filename', type=str, default="RESCANlibsSue_sxn.sam", help='The filename of the SAM file')
	args = parser.parse_args()
	return(args.p,args.f)

def read_counts(infilename,outfilename):
	#SOLEXA2_1222:7:1:7773:4588#0/2_CCAAG	16	AtChr4	3982487	0	30M	*	0	0	TTTAAGTTCTTATACTCAATCATACACATG	::BFFGGGGGG?FEGBEGGGGGGGGB?E:8	XT:A:R	NM:i:0	X0:i:235	XM:i:0	XO:i:0	XG:i:0	MD:Z:30
	pattern = r'^SOLEXA.*\tAtChr(\d{1})\t(\d*)\t.*$'
	recomp = re.compile(pattern)
	#SOLEXA3_0107:8:120:17811:15184#0/1_AGTGT	0	AlChr14	6984372	0	42M	*	0	0	CATGGTAGAAGCAGCCTAGAAGTTCCTTCTTCAGGGTGGTTG	B?:5BGGGGGGBGEEGBGEDGGDGGGGGGGGG@DGE@GD;<:	XT:A:R	NM:i:1	X0:i:3	X1:i:6	XM:i:1	XO:i:0	XG:i:0	MD:Z:21C20
	pattern2 = r'^SOLEXA.*\tAlChr[11-18]}\t\d*\t.*$'
	recomp2 = re.compile(pattern2)
	#SOLEXA2_1222:7:1:10178:3408#0/2_CCAAG	4	*	0	0	*	*	0	0	CATGTTTAGAATATAATGAAATTCAAGAGGTA	D?D?DGIII:GG:GGGGG@GEGGGEE@GGG=B
	pattern3 = r'^SOLEXA.*\t\*\t\d*\t.*$'
	recomp3 = re.compile(pattern3)
	gene_nuc_dict, reads = get_genes()
	thal_count = 0; thal_gene_count = 0; care_count = 0; nomap_count = 0
	with open(infilename) as samfile:
		for line in samfile:
			match = recomp.search(line)
			if match:
				thal_count += 1
				chr_num = str(match.group(1)); start_pos = int(match.group(2));
				if chr_num == '6': chr_num = 'C'
				elif chr_num == '7': chr_num = 'M'
				gene_info = chr_num, start_pos
				if gene_info in gene_nuc_dict:
					thal_gene_count += 1
					gene_name = gene_nuc_dict[gene_info]
					#readcount = int(reads[gene_name].pop(0)) + 1
					#reads[gene_name].insert(0,readcount)
					reads[gene_name][0] += 1
			else:
				match = recomp2.search(line)
				if match:
					care_count += 1
				else:
					match = recomp3.search(line)
					if match:
						nomap_count += 1
	
	# Incorporate mismapped reads that should truly map within At-Chr.3
	
	pattern4 = r'^>Query\d+\tChr3:(\d+)-(\d+)\t(.*)$'
	recomp4 = re.compile(pattern4)
	filename = "Data/AtReadsMismappedToAl.fasta"
	with open(filename) as infile:
		linelist = infile.readlines()
		for line in linelist:
			match = recomp4.match(line)
			if match:
				start_pos = int(match.group(1))
				end_pos = int(match.group(2))
				#reason = str(match.group(3))
				checked_genes = []
				for i in range(min(start_pos,end_pos),max(start_pos,end_pos)):
					if (3, i) in gene_nuc_dict:
						gene_name = gene_nuc_dict[(3, i)]
						if gene_name not in checked_genes:
							reads[gene_name][0] += 1
							checked_genes.append(gene_name)
	
	with open(outfilename, 'w') as outfile:
		totalcount = thal_count + care_count + nomap_count
		line = "Gene\tChr\tStart_Pos\tCounts\nTOTALREAD\t" + str(totalcount) + "\nTHALCOUNT\t" + str(thal_count) + "\nTHALGENEC\t" + str(thal_gene_count) + "\nLYRCOUNTS\t" + str(care_count) + "\nNOMAPPING\t" + str(nomap_count) + "\n"
		outfile.write(line)
		gene_list = sorted(reads)
		for gene_name in gene_list:
			[read_count, chr_num, start_pos] = reads[gene_name]
			line = str(gene_name) + "\t" + str(chr_num) + "\t" + str(start_pos) + "\t" + str(read_count) + "\n"
			outfile.write(line)

start = time.time()

path, filename = options()
infilename = str(path[1:]) + str(filename[1:])
outfilename = "Output/" + str(filename[1:-4]) + "_readcount.txt"
read_counts(infilename,outfilename)

total = (time.time() - start) / 60
print("\nDuplication search for " + str(filename[1:]) + " completed (" + str(total) + " minutes)")