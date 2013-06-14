#!/usr/bin/env python3.2
import argparse, subprocess, os, re, types
from random import randrange, sample
from sys import executable as pypath
import SAM

# RandDupCov.py creates and outputs a randomly generated read profile for a given chromosome or chromosomes. Read lengths and overall read coverage are user-specified values,
#	  as are the number and length of duplicated regions to be randomly incorporated into the read profile. This function utilizes ART, a read simulation program,
#	  which is also used in the 1000 human genomes project.

def command_line():
    parser = argparse.ArgumentParser(description="This program can perform two functions:\n(1) Run a duplication search HMM on a SAM file in an attempt to detect duplications. \
                                     Through a pipeline, the SAM file is automatically broken up into user-specified window sizes, reads mapping per-window counted, and a Poisson \
                                     regression is performed on the resulting count data. This Poisson lambda value is then used to generate an HMM parameter file, and \
                                     the duplication search HMM is run using this parameter file. Alternatively, the user can provide a lambda value to be used and regression will be skipped.\n\
                                     (2) Create and output a randomly generated read profile for a given chromosome. Read lengths and overall read coverage are user-specified values, \
                                     as are the number and length of duplicated regions to be randomly incorporated into the read profile. This function utilizes ART, a read \
                                     simulation program, which is also used in the 1000 human genomes project.")
    parser.add_argument('-sam', default="Output/HMMCovWin_RandomAtChr3Reads/RandomAtChr3Reads.sam.gz", type=str, help='The output SAM filename. (Default="Output/HMMCovWin_RandomAtChr3Reads/RandomAtChr3Reads.sam")', metavar='SAMFilename')
    parser.add_argument('-ref', default="../at9_genome.fasta", help='Reference sequence from which reads shall be simulated. (Default="../at9_genome.fasta")', metavar='ReferenceSeq')
    parser.add_argument('-chr', default="All", help='Chromosome(s) from the reference sequence on which reads shall be simulated. (Default=Entire genome. To simulate >1 chromosome, comma-delimit chromosomes. E.g. -chr "Chr1,Chr3"', metavar='ReferenceChr')
    parser.add_argument('-rl', default=94, type=int, help='Length of reads after barcode removal (default = 94 bps).', metavar='ReadLength')
    parser.add_argument('-rc', default=5, type=float, help='Overall fold read coverage for the library (default = 5).', metavar='ReadCoverage')
    parser.add_argument('-dn', default=30, type=int, help='Number of duplicated regions to incorporate into randomly generated read profile (default=3).', metavar='NumberofDuplications')
    parser.add_argument('-dl', default=2500, type=int, help='The length in bps for which a duplicated region should run (default=10000).', metavar='DuplicationLength')
    
    args = parser.parse_args()
    sam = args.sam
    ref = args.ref
    chromosome = args.chr
    read_len = args.rl
    read_cov = args.rc
    dup_num = args.dn
    dup_len = args.dl
    
    return(sam, ref, chromosome, read_len, read_cov, dup_num, dup_len)

def get_chr(ref):
    chr_lengths = {}
    chromosome = ""
    chr_seq = []
    with open(ref) as infile:
        for line in infile:
            linelist = line.split()
            if line[:1] == ">":
                if chromosome != "":
                    chr_lengths[chromosome] = ''.join(chr_seq)
                chromosome = linelist[0][1:]
                chr_seq = []
            else:
                chr_seq.append(linelist[0])
    chr_lengths[chromosome] = ''.join(chr_seq)
    return(chr_lengths)

def simulate_reads(ref, chromosome, read_len, read_cov, dup_num, dup_len, sam, chr_len):
    final_SAM_fns = [] # List containing all SAM files to be combined together. First contains regular read simulation SAM file, with regional duplication SAM files added later
    chr_list = []
    pattern = "Chr\d{1,}"
    recomp = re.compile(pattern)
    if chromosome == "All":
        for chr_temp in chr_len:
            match = recomp.match(chr_temp)
            if match:
                chr_list.append(chr_temp)
    else:
        chr_list = [entry for entry in chromosome.split(',')]
    for chromosome in chr_list:
        for i in range(0,len(chr_len[chromosome]),10000):
            out_fn = str(sam[:-3]) + "_temp_" + str(chromosome) + "_" + str(i)
            new_fasta = ">" + str(chromosome) + "\n" + chr_len[chromosome][i:i+10000]
            temp_fasta = "regtemp_" + str(i) + ".fasta"
            with open(temp_fasta, 'w') as outfile:
                outfile.write(new_fasta)
            params = ' '.join(["art_illumina.exe", "-i", str(temp_fasta), "-l", str(read_len), "-f", str(read_cov), "-o", out_fn, "-sam", "-q"])
            #Example usage: art_illumina.exe -i ../at9_chr3.fasta -l 34 -f 2 -o Output/RandomAtReads -sam
            simulation = subprocess.Popen(params)
            simulation.wait()
            new_out_fn = out_fn + ".sam"
            if i > 0:
                reg_SAM = SAM.parse(new_out_fn, True)
                for j in range(0,len(reg_SAM.sam_list)):
                    reg_SAM.sam_list[j][3] = str(int(reg_SAM.sam_list[j][3]) + i) # Changing SAM file read numbers
                reg_SAM.output(new_out_fn)
            os.remove(temp_fasta)
            os.remove(out_fn + ".fq")
            os.remove(out_fn + ".aln")
            final_SAM_fns.append(new_out_fn) # Regular read profile files assigned here

    prev_sel = []
    for i in range(1,dup_num+1):
        spos = 0
        undupped_region = False
        ok_regions = 0
        while undupped_region == False:
            chromosome = ''.join(sample(chr_list,1))
            spos = randrange(0,len(chr_len[chromosome]))
            epos = spos + dup_len
            if epos > len(chr_len[chromosome]): continue
            for entry in prev_sel:
                (c, start, end) = entry
                if chromosome == c:
                    if not (int(spos) >= int(start) and int(spos) <= int(end)) or (int(epos) >= int(start) and int(epos) <= int(end)):
                        # Testing to make sure that newly selected duplicated region is not within a region already selected to be duplicated
                        ok_regions += 1
                    else:
                        break
                else:
                    ok_regions += 1
            if ok_regions == len(prev_sel):
                undupped_region = True
                prev_sel.append((chromosome,spos,epos))
        new_fasta = ">" + str(chromosome) + "\n" + chr_len[chromosome][spos:epos]
        temp_fasta = "duptemp_" + str(i) + ".fasta"
        with open(temp_fasta, 'w') as outfile:
            outfile.write(new_fasta)
        params = ' '.join(["art_illumina.exe", "-i", str(temp_fasta), "-l", str(read_len), "-f", str(read_cov), "-o", str(temp_fasta[:-6]), "-sam", "-q"])
        simulation = subprocess.Popen(params)
        simulation.wait()
        dup_fn = temp_fasta[:-6] + ".sam"
        dup_SAM = SAM.parse(dup_fn, True)
        for i in range(0,len(dup_SAM.sam_list)):
            dup_SAM.sam_list[i][3] = str(int(dup_SAM.sam_list[i][3]) + spos)
        dup_SAM.output(dup_fn[:-4] + "_temp.sam")
        os.remove(temp_fasta)
        os.remove(temp_fasta[:-6] + ".fq")
        os.remove(temp_fasta[:-6] + ".aln")
        os.remove(dup_fn)
        final_SAM_fns.append(dup_fn[:-4] + "_temp.sam")
    
    final_SAM = SAM.parse(final_SAM_fns, True)
    pattern2 = r"Chr\d{1}-(\d+)"
    recomp2 = re.compile(pattern2)
    
    max_read = final_SAM.sam_list[0][0]
    match = recomp2.match(max_read)
    max_num = int(match.group(1))
    min_met = False
    for i in range(0, len(final_SAM.sam_list)): # Ensure that read names from one SAM file do not coincide with read names from another SAM file
        r_name = final_SAM.sam_list[i][0]
        match = recomp2.match(r_name)
        if match:
            r_num = int(match.group(1))
            if r_num == 1 and min_met == False:
                min_met = True
                continue
            if min_met == True:
                max_num += 1
                line = final_SAM.sam_list[i]
                line[0] = r_name[:5] + str(max_num)
                final_SAM.sam_list[i] = line
    final_SAM.header = sorted(final_SAM.header, key = lambda read: read[1][6:])
    final_SAM.sam_list = sorted(final_SAM.sam_list, key = lambda read: int(read[0][5:]))
    for i in range(0,len(final_SAM.header)): # Set chromosome length in header portion of SAM file to correct length (will be 10,000 in temporary simulation SAM files)
        chromosome = final_SAM.header[i][1][3:]
        final_SAM.header[i][2] = "LN:" + str(len(chr_len[chromosome]))
    final_SAM.output(sam[:-3])
    for file in final_SAM_fns:
        os.remove(file)
    params = ' '.join([str(pypath) + " FixARTSAMFile.py", "-i", sam[:-3], "-o", sam])
    fixSAM = subprocess.Popen(params) # Replaces new CIGAR string format for matches, which uses = and X, to the old format M for both
    fixSAM.wait()
    print("Duplicated regions are located at:\n")
    for entry in sorted(prev_sel, key = lambda entry:(entry[0], entry[1])):
        (chromosome, spos, epos) = entry
        print(str(chromosome), ":", str(spos),"-", str(epos), sep="")

sam, ref, chromosome, read_len, read_cov, dup_num, dup_len = command_line()
chr_len = get_chr(ref)
simulate_reads(ref, chromosome, read_len, read_cov, dup_num, dup_len, sam, chr_len)