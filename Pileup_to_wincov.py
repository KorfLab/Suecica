#!/usr/bin/env python3.2
import argparse
import sys; sys.path.append('/home/mattchat/Packages/')
from GenomeTools import Pileup
from collections import Counter
import os, gzip

def command_line():
    parser = argparse.ArgumentParser(description="Uses a pileup file to calculate and output total coverage along each chromosome in user-defined window intervals")
    parser.add_argument('-f', default='Data/comai_col/Col2x#1/libCol2x#1_pileup_simplified.txt.gz', type=str, help='Name of pileup file to be parsed (default="Data/comai_col/Col2x#1/libCol2x#1_pileup_simplified.txt.gz")', metavar='Pileup')
    parser.add_argument('-r', default='../at9_genome.fasta', type=str, help='Filename for reference genome (default="../at9_genome.fasta")', metavar='Reference')
    parser.add_argument('-w', default=10000, type=int, help='Window size to be used in bps (default=10000)', metavar='Window')
    args = parser.parse_args()
    pileup = args.f
    ref = args.r
    win = args.w
    return(pileup, ref, win)

def get_chr_len(ref):
    chr_lengths = Counter()
    chr = ""
    chr_len = 0
    with open(ref) as infile:
        for line in infile:
            linelist = line.split()
            if line[:1] == ">":
                if chr != "":
                    chr_lengths[chr] = chr_len
                chr = linelist[0][1:]
                chr_len = 0
            else:
                line = line.split()
                chr_len += len(line[0])
    chr_lengths[chr] = chr_len
    return(chr_lengths)

pileup, ref, win = command_line()
chr_lengths = get_chr_len(ref)

PileupFile = Pileup.parse(pileup,True)

#windows = [i for i in range(1,1501)]
windows = [i for i in range(1600,5001,100)]
for win in windows:
    outdir = "Output/HMMCovWin/" + str(win) + "bp-window/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for chr in sorted(chr_lengths):
        histogram = Counter()
        outlist = []
        for i in range(1,chr_lengths[chr]+1,win):
            loc = str(chr) + ":" + str(i) + "-" + str(i+win)
            cov = PileupFile.Coverage(loc,True)
            if cov in histogram: histogram[cov] += 1
            else: histogram[cov] = 1
            outlist.append(str(cov))
        
        outline = " ".join(outlist)
        outfilename1 = str(outdir) + str(chr) + "_obs.txt"
        with open(outfilename1, 'w') as outfile:
            outfile.write(outline)
        
        maxcov = max([cov for cov in histogram.keys()])
        for i in range(1,maxcov):
            if i not in histogram:
                histogram[i] = 0
        
        # Removing values from histogram starting where we have 4 consecutive 0's
        kmer_window = 4
        hist_len = len(histogram.keys())
        for i in range(0,maxcov):
            zerosum = 0
            for j in range(i,i+kmer_window):
                zerosum += histogram[j]
            #if zerosum == 0 and i >= 300: # For larger windows, about >900 bps, curve shifts right, and a lot of zeroes occur in the beginning. >= 300 coverage values ensures histogram is actually created instead of stopping at these early series of 0 coverage windows
            if zerosum == 0 and i >= 2000: # For larger windows, about >150 bps, curve shifts even further right, and a lot of zeroes occur in the beginning. >= 2000 coverage values ensures histogram is actually created instead of stopping at these early series of 0 coverage windows
                for j in range(i,maxcov+1):
                    del histogram[j]
                break
        
        outline = "\n".join([str(cov) + "\t" + str(freq) for cov, freq in sorted(histogram.items(), key = lambda cov: cov[0])])
        outfilename2 = str(outdir) + str(chr) + "_histogram.txt"
        with open(outfilename2, 'w') as outfile:
            outline2 = "Coverage_per_" + str(win) + "bp_window\tCounts\n"
            outfile.write(outline2)
            outfile.write(outline)