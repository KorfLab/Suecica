#!/usr/bin/env python
import re, time

start = time.time()

infilename = "Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/RESCANlibsSue_sxn.sam"
outfilename = "Output/AsAlReadsMappedToAt/AsAlReads.fasta"
#SOLEXA3_0107:8:120:17811:15184#0/1_AGTGT	0	AlChr14	6984372	0	42M	*	0	0	CATGGTAGAAGCAGCCTAGAAGTTCCTTCTTCAGGGTGGTTG	B?:5BGGGGGGBGEEGBGEDGGDGGGGGGGGG@DGE@GD;<:	XT:A:R	NM:i:1	X0:i:3	X1:i:6	XM:i:1	XO:i:0	XG:i:0	MD:Z:21C20
pattern = r'^SOLEXA.*\tAlChr(\d{2})\t(\d*)\t.*\t.*\t.*\t.*\t.*\t([ACTGKMRYSWBVHDXN]*)\t.*$'
recomp = re.compile(pattern)
i = 0
outlinelist = []
with open(infilename) as samfile:
	linelist = samfile.readlines()
	for line in linelist:
		match = recomp.match(line)
		if match:
			i += 1
			chr = str(match.group(1))
			pos = str(match.group(2))
			nucs = str(match.group(3))
			line = ">SEQ_" + str(i) + "_Chr" + str(chr) + "_" + str(pos) + "\n" + str(nucs) + "\n"
			outlinelist.append(line)

with open(outfilename, 'w') as outfile:
	outlines = ''.join(outlinelist)
	outfile.write(outlines)

total = (time.time() - start) / 60
print("\nOutput of A. lyrata-mapping reads completed (" + str(total) + " minutes)")