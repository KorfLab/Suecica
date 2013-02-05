#!/usr/bin/env python3
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import argparse, re, os
from collections import Counter

def Command_line():
    parser = argparse.ArgumentParser(description="Merge_Expression_Files.py combines the reads mapped per gene for all data files in the provided folder into a single table.")
    parser.add_argument('-i', default="Data/RNA-Seq/AsLeaf_vs_ColLeaf", type=str, help='Input folder for mRNA expression data files', metavar='Input')
    parser.add_argument('-o', default="Data/RNA-Seq/AsLeaf_vs_ColLeaf/mRNA_Expression.tsv", type=str, help='Output file for combined mRNA expression data files.', metavar='Output')
    
    args = parser.parse_args()
    indir = args.i
    outfile = args.o
    
    return(indir, outfile)

indir, outfile = Command_line()
pattern = r"(.*_gene_count.tsv$)"
recomp = re.compile(pattern)

gene_dict = {}
fname_list = []
for dirname, dirnames, filenames in os.walk(indir):
    for f_name in filenames:
        match = recomp.match(f_name)
        if match:
            fname = match.group(1)
            fname_list.append(fname)
            
for fname in fname_list:
    open_file_name = os.path.join(indir, str(fname))
    with open(open_file_name) as infile:
        for line in infile:
            line = line.strip()
            line = line.split("\t")
            gene = line[0]
            count = line[1]
            if line[0] != "Gene":
                try:
                    gene_dict[gene][fname] = str(line[1])
                except:
                    gene_dict[gene] = {x: str(0) for x in fname_list}
                    gene_dict[gene][fname] = str(count)

with open(outfile, 'w') as outfile:
    fname_outline = "Gene\t" + "\t".join(fname_list) + "\n"
    outfile.write(fname_outline)
    for gene in sorted(gene_dict):
        outline = str(gene) + "\t"
        for fname in fname_list:
            outline += gene_dict[gene][fname] + "\t"
        outline += "\n"
        outfile.write(outline)