#!/usr/bin/env python3.2
from collections import Counter
import sys, os, argparse, re, subprocess

def command_line():
    parser = argparse.ArgumentParser(description="Uses histogram data files to calculate and output total coverage along each chromosome in user-defined window intervals")
    parser.add_argument('-d', default='Output/HMMCovWin/', type=str, help='Directory containing bp window folders (default="Output/HMMCovWin/")', metavar='Directory')
    parser.add_argument('-r', default='../at9_genome.fasta', type=str, help='Filename for reference genome (default="../at9_genome.fasta")', metavar='Reference')
    args = parser.parse_args()
    hist_folder = args.d
    ref = args.r
    return(hist_folder, ref)

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
    del chr_lengths["chloroplast"]  # Removing chromosomes unlikely to produce coverage useful for duplication study. Histograms do not work for windows larger than a few bps; too little sampling.
    del chr_lengths["mitochondria"] # Removing chromosomes unlikely to produce coverage useful for duplication study. Histograms do not work for windows larger than a few bps; too little sampling.
    return(chr_lengths)

hist_folder, ref = command_line()
chr_lengths = get_chr_len(ref)
folders = os.listdir(hist_folder)
pattern = r"(\d+)bp-window"
recomp = re.compile(pattern)
pattern2 = r"(\d+)bp-window/(Chr\d{1})_hist_distfits.txt$"
recomp2 = re.compile(pattern2)

bp_windows = []
for folder in folders:
    match = recomp.match(folder)
    if match:
        bp_windows.append(int(match.group(1)))
bp_windows = sorted(bp_windows)

filelist = []
for window in bp_windows:
    for chr in sorted(chr_lengths.keys()):
        filelist.append(''.join(["Output/HMMCovWin/",str(window),"bp-window/",chr,"_hist.txt"]))

outfilename = "Output/HMMInputfiles.txt"
line = '\n'.join(filelist)
with open(outfilename, 'w') as outfile:
    outfile.write(line)

Rpath = r"C:\Program Files\R\R-2.14.1\bin\x64\Rscript.exe"
params = ' '.join([str(Rpath) + " dist_fit.R", "--no-save"])
simulation = subprocess.Popen(params)
simulation.wait()

filelist = []
for window in bp_windows:
    for chr in sorted(chr_lengths.keys()):
        filelist.append(''.join(["Output/HMMCovWin/",str(window),"bp-window/",chr,"_hist_distfits.txt"]))

bestfit_dict = {}
for file in filelist:
    window_match = recomp2.search(file)
    keytuple = (window, chr) = window_match.groups()
    with open(file) as infile:
        inlines = infile.readlines()
        inlines = [line.split() for line in inlines]
        pois_lambda = float(inlines[1][1])
        pois_negtwologlik = float(inlines[2][1])
        expo_decay = float(inlines[3][1])
        expo_negtwologlik = float(inlines[4][1])
        geom_decay = float(inlines[5][1])
        geom_negtwologlik = float(inlines[6][1])
        norm_mean = float(inlines[7][1])
        norm_stdev = float(inlines[8][1])
        norm_negtwologlik = float(inlines[9][1])
        data = (pois_lambda, pois_negtwologlik, expo_decay, expo_negtwologlik, geom_decay, geom_negtwologlik, norm_mean, norm_stdev, norm_negtwologlik)
        bestfit_dict[keytuple] = data

for chr_key in sorted(chr_lengths.keys()):
    outline = ["bp_window\tchr\tPoisson_lambda\tPoisson_-2loglik\tExpo_decay\tExpo_-2loglik\tGeom_decay\tGeom_-2loglik\tNormal_mean\tNormal_stdev\tNormal-2loglik"]
    for entry in sorted(bestfit_dict.items(), key = lambda entry:(entry[1][1])):
        (window, chr), [pois_lambda, pois_negtwologlik, expo_decay, expo_negtwologlik, geom_decay, geom_negtwologlik, norm_mean, norm_stdev, norm_negtwologlik] = entry
        if chr == chr_key:
            line = '\t'.join([window, chr, str(pois_lambda), str(pois_negtwologlik), str(expo_decay), str(expo_negtwologlik), str(geom_decay), str(geom_negtwologlik), str(norm_mean), str(norm_stdev), str(norm_negtwologlik)])
            outline.append(line)
    outline = '\n'.join(outline)
    
    outfilename = ''.join(["Output/HMMCovWin/DistributionFits_", chr_key,".txt"])
    with open(outfilename,'w') as outfile:
        outfile.write(outline)

outline = ["bp_window\tchr\tPoisson_lambda\tPoisson_-2loglik\tExpo_decay\tExpo_-2loglik\tGeom_decay\tGeom_-2loglik\tNormal_mean\tNormal_stdev\tNormal-2loglik"]
for entry in sorted(bestfit_dict.items(), key = lambda entry:(entry[1][1])):
    (window, chr), [pois_lambda, pois_negtwologlik, expo_decay, expo_negtwologlik, geom_decay, geom_negtwologlik, norm_mean, norm_stdev, norm_negtwologlik] = entry
    line = '\t'.join([window, chr, str(pois_lambda), str(pois_negtwologlik), str(expo_decay), str(expo_negtwologlik), str(geom_decay), str(geom_negtwologlik), str(norm_mean), str(norm_stdev), str(norm_negtwologlik)])
    outline.append(line)
outline = '\n'.join(outline)

outfilename = ''.join(["Output/HMMCovWin/DistributionFits.txt"])
with open(outfilename,'w') as outfile:
    outfile.write(outline)