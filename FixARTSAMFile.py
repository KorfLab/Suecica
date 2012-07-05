#!/usr/bin/env python3.2
import re, argparse, gzip, os

def command_line():
    parser = argparse.ArgumentParser(description="Fixes the SAM file output from the read simulation software ART")
    parser.add_argument('-i', default="Output\HMMCovWin_RandomAtChr3Reads\RandomAtChr3Reads.sam", type=str, help='Filename of SAM file to be modified (default="Output\HMMCovWin_RandomAtChr3Reads\RandomAtChr3Reads.sam")', metavar='Filename')
    parser.add_argument('-o', default="", type=str, help='Filename of output modified SAM file (default=InputFilename + "Fixed.sam.gz")', metavar='Filename')
    args = parser.parse_args()
    i_fn = args.i
    o_fn = args.o
    if o_fn == "": o_fn = SAM[:-4] + "Fixed" + ".sam"
    return(i_fn, o_fn)

i_fn, o_fn = command_line()
recomp = re.compile("=|X")

outlines = []
with open(i_fn) as infile:
    for line in infile:
        line = line.split()
        if line[0][:1] != "@":
            for m in recomp.finditer(line[5]):
                s_pos = m.start()
                line[5] = line[5][:s_pos] + "M" + line[5][s_pos+1:]
        outlines.append(line)

with open(o_fn, 'w') as outfile:
    with gzip.open(o_fn, 'wb') as outfile2:
        for line in outlines:
            outline = '\t'.join(line) + "\n"
            outfile.write(line)
            outfile2.write(bytes(outline,"UTF-8"))