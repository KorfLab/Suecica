#!/usr/bin/env python
import re, time

start = time.time()

infilename = "Data\TAIR9_seq_20090619"
outfilename = "Output\RESiteCounts.txt"
pattern = r'^>(.{9,12}) | .*$'
pattern2 = r'CATG|GTAC'
recomp = re.compile(pattern)
recomp2 = re.compile(pattern2)
gene_dict = {}
with open(infilename) as infile:
    gene_seq = ""
    for line in infile:
        match = recomp.search(line)
        if match:
            gene_name = str(match.group(1))
            for line in infile:
                match = recomp.search(line)
                if not match:
                    gene_seq += line
                else:
                    cut_count = 0
                    matches = recomp2.findall(gene_seq)
                    cut_count = len(matches)
                    gene_dict[gene_name] = cut_count
                    gene_name = str(match.group(1))
                    gene_seq = ""
                    break
        else:
            cut_count = 0
            gene_seq += line
            for line in infile:
                match = recomp.search(line)
                if not match:
                    gene_seq += line
                else:
                    matches = recomp2.findall(gene_seq)
                    cut_count = len(matches)
                    gene_dict[gene_name] = cut_count
                    gene_name = str(match.group(1))
                    gene_seq = ""
                    break

noisoforms = [key[0:9] for key in gene_dict]
isoformdict = {}
for key in gene_dict:
    if key[0:9] not in isoformdict: isoformdict[key[0:9]] = noisoforms.count(key[0:9])
        
global_count = 0
gene_count = 0
for key in sorted(isoformdict):
    isolist = []
    added = 0
    i = 0
    max = 0
    maxiso = ""
    while added < isoformdict[key]:
        i += 1
        newiso = str(key) + "." + str(i)
        if newiso in gene_dict:
            isolist.append(newiso)
            added += 1

    maxiso = isolist[0]
    for iso in isolist:
        if gene_dict[iso] > max:
            max = gene_dict[iso]
            maxiso = iso

    gene_dict[key] = gene_dict[maxiso]
    global_count += gene_dict[maxiso]
    gene_count += 1
    for iso in isolist:
        del gene_dict[iso]

with open(outfilename, 'w') as outfile:
    avg_count = global_count / gene_count
    line = "Average\t" + str(avg_count) + "\n"
    outfile.write(line)
    for gene in sorted(gene_dict):
        line = str(gene) + "\t" + str(gene_dict[gene]) + "\n"
        outfile.write(line)

total = (time.time() - start) / 60
print("\nPer-gene restriction enzyme sites counted (" + str(total) + " minutes)")