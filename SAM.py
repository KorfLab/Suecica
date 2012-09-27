#!/usr/bin/env python3.2
import re, gzip
from collections import Counter

# SAM.py package is used in parsing SAM files, capable of outputting coverage, coverage-per-window,
# trimming reads based on read length, and outputting read length histograms. This package can also
# combine multiple SAM files together into a single file.


def Chr_Lengths(filename):
    ''' Opens a SAM file and returns chromosome names and the lengths of those chromosomes.'''
    chr_len = {}
    SQ_header_began = False
    if filename[-3:] == ".gz":
        infile = gzip.open(filename, 'rb')
    else:
        infile = open(filename)
    for line in infile:
        if filename[-3:] == ".gz": line = str(line, encoding='utf8')
        line = line.split("\t")
        if line[0] == "@SQ":
            SQ_header_began = True
            chr = line[1][3:]
            len = int(line[2][3:])
            chr_len[chr] = len
        elif SQ_header_began == True:
            break
    infile.close()
    return(chr_len)

def Combine(filelist, outfilename):
    ''' Combines SAM files into a single file, each separated by a semicolon.
        Among header lines, only @SQ SN:ChrName  LN:NumberHere are kept from the first SAM file given. '''
    filelist = filelist.split(";")
    print("Combining ",len(filelist), " SAM files into ", str(outfilename),sep='')
    gzipped_in = False
    gzipped_out = False
    if filelist[0][-3:] == ".gz":
        outfile = gzip.open(outfilename, 'wb')
        gzipped_out = True
    else:
        outfile = open(outfilename, 'w')
    # Open first file separately from other SAM files; only grabbing header lines from first file
    if filelist[0][-3:] == ".gz":
        infile = gzip.open(filelist[0], 'rb')
        gzipped_in = True
    else:
        infile = open(filelist[0])
    in_header = True
    print("Outputting ", filelist[0],sep='')
    for line in infile:
        if gzipped_in == True: line = str(line, encoding='utf8')
        if in_header == False: # Non-header line
            if gzipped_out == True:
                outfile.write(bytes(line,"UTF-8"))
            else:
                outfile.write(line)
        else:
            if line[0:1] != "@":
                if line[0:3] == "@SQ":
                    if gzipped_out == True:
                        outfile.write(bytes(line,"UTF-8"))
                    else:
                        outfile.write(line)
            else:
                in_header = False
                if gzipped_out == True:
                    outfile.write(bytes(line,"UTF-8"))
                else:
                    outfile.write(line)
    infile.close()
    # Open remaining SAM files
    for i in range(1,len(filelist)):
        gzipped_in = False
        if filelist[i][-3:] == ".gz":
            infile = gzip.open(filelist[i], 'rb')
            gzipped_in = True
        else:
            infile = open(filelist[i])
        in_header = True
        print("Outputting ", filelist[i],sep='')
        for line in infile:
            if gzipped_in == True: line = str(line, encoding='utf8')
            if in_header == False: # Non-header line
                if gzipped_out == True:
                    outfile.write(bytes(line,"UTF-8"))
                else:
                    outfile.write(line)
            else:
                if line[0:1] != "@":
                    in_header = False
                    if gzipped_out == True:
                        outfile.write(bytes(line,"UTF-8"))
                    else:
                        outfile.write(line)
        infile.close()
    outfile.close()
    print("Combining of SAM files complete!")

def Coverage(sam_files):
    ''' Calculates coverage for entire SAM file.'''
    chr_len = {}
    sam_files = sam_files.split(";")
    for sam_file in sam_files:
        total_count = {}
        print("\nCalculating coverage for ",str(sam_file)," ...\n",sep='')
        if sam_file[-3:] == ".gz":
            infile = gzip.open(sam_file, 'rb')
        else:
            infile = open(sam_file)
        
        for line in infile:
            if sam_file[-3:] == ".gz": line = str(line, encoding='utf8')
            line = line.split("\t")
            if line[0][0] != "@":
                chr = line[2]
                if chr == "*": continue
                s_pos = line[3]
                CIGAR = line[5]
                matches = re.findall(r"(\d+)M",CIGAR)
                match_count = sum([int(match) for match in matches])
                try:
                    total_count[chr] += match_count
                except:
                    total_count[chr] = match_count
            else:
                if line[0] == "@SQ":
                    chr = line[1][3:]
                    len = int(line[2][3:])
                    chr_len[chr] = len
        infile.close()
        
        for (chr, total_pos) in chr_len.items():
            print(str(chr), " coverage: ", str(round(total_count[chr] /  total_pos,2)), "X", sep='')
        print("Overall coverage: ", str(round(sum(total_count.values()) / sum(chr_len.values()),2)), "X\n", sep='')

def Coverage_Window(sam_file, win_list, chr_list, control_use):
    ''' Returns histogram of read coverage on a per-bp-window basis for a given chromosome '''
    
    chr_list = chr_list.split(";")
    
    # Initialize a multi-dimensional per-window, per-chromosome series of read count HMM paths
    win_db = {int(win): {} for win in win_list}
    for win in win_list:
        win_db[win] = {chr: Counter() for chr in chr_list}
    
    # Read in SAM file and create reads-per-window histograms
    total_reads = 0
    if sam_file[-3:] == ".gz":
        infile = gzip.open(sam_file, 'rb')
    else:
        infile = open(sam_file)
    for line in infile:
        if sam_file[-3:] == ".gz": line = str(line, encoding='utf8')
        split_line = line.split("\t")
        chr = str(split_line[2])
        if chr not in chr_list: continue # Ignore if read does not map, or if mapped to a chromosome not analyzed (e.g. mito, chloro)
        quality = int(split_line[4])
        if quality == 0: continue # Ignore reads with a quality of 0, aka those which do not map uniquely
        total_reads += 1
        spos = int(split_line[3])
        read_len = len(split_line[9])
        epos = spos + read_len
        for win in win_list:
            spos_for_dict = int(spos // win)
            epos_for_dict = int(epos // win)
            try:
                win_db[win][chr][spos_for_dict] += 1
            except KeyError:
                win_db[win][chr][spos_for_dict] = 1
            if spos_for_dict != epos_for_dict:
                total_reads += 1
                try:
                    win_db[win][chr][epos_for_dict] += 1
                except KeyError:
                    win_db[win][chr][epos_for_dict] = 1
    infile.close()
        
    # For all read-per-window counts that did not exist, assign them a value of 0
    for chr in chr_list:
        for win in win_list:
            for i in range(0,max(win_db[win][chr].keys())):
                if i not in win_db[win][chr].keys(): win_db[win][chr][i] = 0
    
    # If file is to be used to find exp:control ratio, increment all read-per-window counts by 1 (to ensure no division by 0).
    # Then, calculate % of total reads for each window
    if control_use == "True":
        for chr in chr_list:
            for win in win_list:
                win_db[win][chr] = {i: (j+1)/total_reads for i, j in win_db[win][chr].items()}
    
    # Create multi-dimensional per-window, per-chromosome dictionaries for paths and histograms
    path_dict = {win: {} for win in win_list}
    hist_dict = {win: {} for win in win_list}
    
    # Fill in path and histogram dictionaries with data, then return them to the calling function
    for win in win_list:
        for chr in chr_list:
            path = ' '.join([str(cov) for (win_num, cov) in sorted(win_db[win][chr].items(), key = lambda win_num: win_num)])
            path_dict[win][chr] = path
            if control_use == "False":
                hist = Counter({cov: 0 for cov in range(0,max(win_db[win][chr].values() ) + 1) } )
                for count in win_db[win][chr].values():
                    hist[count] += 1
                hist_dict[win][chr] = hist
            else:
                hist_dict = {} # Only path_dict is used when control_use==True
    return(hist_dict, path_dict)

def Reads_Histogram(filelist):
    '''Outputs a histogram of read length in the given SAM files'''
    filelist = filelist.split(";")
    for filename in filelist:
        print("Outputting read length histogram for ", str(filename), " ...", sep='')
        hist = Counter()
        gzipped = False
        if filename[-3:] == ".gz":
            infile = gzip.open(filename, 'rb')
            gzipped = True
        else:
            infile = open(filename)
        
        for line in infile:
            if gzipped == True: line = str(line, encoding='utf8')
            splitlist = line.split()
            if len(splitlist) > 9:
                length = len(splitlist[9]) # Nucleotides in read
                try:
                    hist[length] += 1
                except KeyError:
                    hist[length] = 1
        infile.close()
        
        outline = "Read length\tCount\n"
        outline2 = ''.join([str(entry) + "\t" + str(hist[entry]) + "\n" for entry in sorted(hist.keys())])
        outfilename = filename[:-7] + "_" + filename[-6:] + "_histogram.txt"
        with open(outfilename, 'w') as outfile:
            outfile.write(outline)
            outfile.write(outline2)
    
def Trim_Reads(filelist,threshold):
    '''Trims a list of SAM files based on a minimum length threshold'''
    filelist = filelist.split(";")
    for filename in filelist:
        print("Trimming ",str(filename)," based on reads with a minimum length of ", str(threshold), "...", sep='')
        gzipped = False
        if filename[-3:] == ".gz":
            outfilename = filename[:-7] + "_trimmed" + filename[-7:]
            infile = gzip.open(filename, 'rb')
            outfile = gzip.open(outfilename, 'wb')
            gzipped = True
        else:
            outfilename = filename[:-4] + "_trimmed" + filename[-4:]
            infile = open(filename, 'r')
            outfile = open(outfilename, 'w')
        for line in infile:
            if gzipped == True: line = str(line, encoding='utf8')
            splitlist = line.split()
            if len(splitlist) > 9:
                if len(splitlist[9]) >= threshold: # Nucleotides in read
                    if gzipped == True:
                        outfile.write(bytes(line,"UTF-8"))
                    else:
                        outfile.write(line)
            else:
                if gzipped == True:
                    outfile.write(bytes(line,"UTF-8"))
                else:
                    outfile.write(line)
        infile.close()
        outfile.close()