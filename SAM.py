#!/usr/bin/env python3
import re, gzip, sys, os
from collections import Counter
import pdb

# SAM.py package is used in parsing SAM files, capable of calculating coverage, coverage-per-window,
# obtaining chromosome lengths, outputting histograms for read length, read mapping quality, and base
# quality, and trimming reads based on read length.

class Parse:
    ''' Parses the provided SAM file and can provide a number of various bits of information '''
    
    def __init__(self, f_in):
        ''' Prepares the provided SAM file for parsing through later function calls '''
        # Store filename
        self.f_in = f_in
        
        # Initialize data types
        self.reset()
        
    def base_quality(self, map_threshold=20, base_threshold=20, o_fn="default"):
        ''' Outputs a histogram of base qualities, the percentage of reads containing
        low-quality bases among reads above a mapping threshold, and the average number
        of low-quality bases observed per read containing at least one low-quality base '''
        self.base_qual_run = True
        self.base_qual_map_thr = map_threshold
        self.base_qual_base_thr = base_threshold
        self.base_qual_hist = Counter()
        self.total_low_reads_low_base_count = 0
        self.total_low_reads = 0
        self.base_qual_total_reads = 0
        
        if o_fn == "default":
            # Set default output filename
            dir_path, f_name_ext = os.path.split(self.f_in)
            f_name, f_ext = os.path.splitext(f_name_ext)
            f_outname = f_name + "_base_qual.txt"
            o_fn = os.path.join(dir_path, f_outname)
        self.base_qual_o_fn = o_fn
    
    def chr_lengths(self):
        ''' Opens a SAM file and returns chromosome names and the lengths of those chromosomes.'''
        self.chr_lengths_run = True
    
    def coverage(self, threshold=20):
        ''' Calculates coverage for each chromosome, as well as for the entire SAM file.'''
        self.coverage_run = True
        self.cov_qual_thr = threshold
        
    def get_chr_lengths(self):
        ''' Returns chromosome lengths following a SAM.start() run '''
        if self.chr_len == {}:
            print("Call chr_lengths() and start() first.")
            sys.exit(0)
        return(self.chr_len)
    
    def get_reads_per_gene(self):
        if len(self.gene_counter.keys()) == 0:
            print("Call reads_per_gene() and start() first.")
            sys.exit(0)
        return(self.gene_counter)
    
    def get_total_reads(self):
        ''' Returns total reads in SAM file '''
        if self.total_reads_count == 0:
            print("Call total_reads() and start() first.")
            sys.exit(0)
        return(self.total_reads_count)
    
    def map_quality(self, threshold=20, o_fn="default"):
        ''' Outputs a histogram of read mapping qualities '''
        self.map_qual_run = True
        self.map_qual_hist = Counter()
        self.map_qual_thr = threshold
        self.map_qual_low_reads = 0
        if o_fn == "default":
            # Set default output filename
            dir_path, f_name_ext = os.path.split(self.f_in)
            f_name, f_ext = os.path.splitext(f_name_ext)
            f_outname = f_name + "_map_qual.txt"
            o_fn = os.path.join(dir_path, f_outname)
        self.map_qual_o_fn = o_fn
    
    def _output(self):
        ''' Provide various forms of command-line and file output '''
        # Print out total reads mapping
        if self.total_reads_run:
            print("There are ", str(self.total_reads_count), " reads (q>=", str(self.map_qual_thr), ") mapping in this library.", sep='')
        
        # Print out library coverage
        if self.coverage_run:
            for (chromosome, total_pos) in sorted(self.chr_len.items()):
                print(str(chromosome), " coverage: ", str(round(self.total_count[chromosome] /  total_pos,2)), "X", sep='')
            print("Overall coverage: ", str(round(sum(self.total_count.values()) / sum(self.chr_len.values()),2)), "X\n", sep='')
        
        # Write out read length histogram
        if self.read_length_hist_run:
            outline = "Read length\tCount\n"
            outline2 = ''.join([str(entry) + "\t" + str(self.read_len_hist[entry]) + "\n" for entry in sorted(self.read_len_hist.keys())])
            with open(self.read_len_o_fn, 'w') as outfile:
                outfile.write(outline)
                outfile.write(outline2)
        
        # Write out read mapping quality histogram
        if self.map_qual_run:
            print("\nTotal reads:\t\t\t\t", str(self.map_qual_total_reads))
            print("Total low mapping quality reads (q<", str(self.map_qual_thr) , "):\t", str(self.map_qual_low_reads), " (",
                  str(round(self.map_qual_low_reads / self.map_qual_total_reads * 100,2)), "%)", sep='')
            
            with open(self.map_qual_o_fn, 'w') as outfile:
                outline = "MapQuality\tCount\n"
                outfile.write(outline)
                for qual, count in sorted(self.map_qual_hist.items()):
                    outline = str(qual) + "\t" + str(count) + "\n"
                    outfile.write(outline)
        
        # Write out base quality histogram
        if self.base_qual_run:
            print("\n\nTotal reads (map_q>=", str(self.base_qual_map_thr), "):\t\t\t\t", str(self.base_qual_total_reads))
            print("Total low quality base-containing reads (base_q<, ", str(self.base_qual_base_thr), "):\t",
                  str(self.total_low_reads), " (", str(round(self.total_low_reads / self.base_qual_total_reads * 100,2)), "%)", sep='')
            print("Average number of low quality bases in each read containing low quality bases:\t",
                  str(round(self.total_low_reads_low_base_count / self.total_low_reads,2)), sep='') 
            
            with open(self.base_qual_o_fn, 'w') as outfile:
                outline = "Quality\tCount\n"
                outfile.write(outline)
                for qual, count in sorted(self.base_qual_hist.items()):
                    outline = str(qual) + "\t" + str(count) + "\n"
                    outfile.write(outline)

    def _process_reads(self, line):
        ''' Grabs all the requested information from SAM file reads '''
        read_name, flag, chromosome, s_pos, map_q, cigar, chr_next, pos_next, temp_len, read_seq, base_q = line[:11]
        s_pos = int(s_pos)
        map_q = int(map_q)
        if chromosome == "*": return()
        # Total up number of reads in library
        if self.total_reads_run:
            if map_q >= self.map_qual_thr:
                self.total_reads_count += 1
        
        # Calculate coverage for each chromosome & library
        if (self.coverage_run or self.total_reads_run):
            add_read = False
            if self.coverage_run and map_q >= self.cov_qual_thr: add_read = True
            elif self.map_qual_run and map_q >= self.map_qual_thr: add_read = True
            if add_read:
                matches = re.findall(r"(\d+)M",cigar)
                match_count = sum([int(match) for match in matches])
                self.total_count[chromosome] += match_count
        
        # Add to read length histogram
        if self.read_length_hist_run:
            read_len = len(read_seq)
            self.read_len_hist[read_len] += 1
        
        # Add to read mapping quality histogram
        if self.map_qual_run:
            if map_q < self.map_qual_thr:
                self.map_qual_low_reads += 1
            self.map_qual_total_reads += 1
            self.map_qual_hist[map_q] += 1
        
        # Add to base quality histogram
        if self.base_qual_run:
            if map_q >= self.base_qual_map_thr:
                is_low_read = False
                for c in base_q:
                    base_qual = ord(c) - 33
                    if base_qual < self.base_qual_base_thr:
                        self.total_low_reads_low_base_count += 1
                        if is_low_read == False:
                            is_low_read = True
                            self.total_low_reads += 1
                    self.base_qual_hist[base_qual] += 1
                self.base_qual_total_reads += 1
                
        # Count up reads mapping per gene
        if self.reads_per_gene_run:
            if map_q >= self.reads_per_gene_map_thr:
                read_len = len(read_seq)
                e_pos = s_pos + read_len
                genes_mapped = []
                for i in range(s_pos, e_pos+1):
                    gene_name = self.gff_genes_dict.get_gene_name(chromosome, i)
                    if gene_name is not None:
                        for gene in gene_name:
                            if gene not in genes_mapped:
                                self.gene_counter[gene] += 1
                                genes_mapped.append(gene)
    
    def read_length_hist(self, o_fn="default"):
        ''' Outputs a histogram of read length in the given SAM files '''
        self.read_length_hist_run = True
        self.read_len_hist = Counter()
        if o_fn == "default":
            # Set default output filename
            dir_path, f_name_ext = os.path.split(self.f_in)
            f_name, f_ext = os.path.splitext(f_name_ext)
            f_outname = f_name + "_read_len_hist.txt"
            o_fn = os.path.join(dir_path, f_outname)
        self.read_len_o_fn = o_fn
    
    def reads_per_gene(self, gff_genes_dict, threshold=20):
        ''' Returns a dictionary histogram of the number of reads mapped to each gene provided in gff_genes_dict.
        Must be provided with a GFF class with create_nuc_dict set to True. '''
        self.reads_per_gene_run = True
        self.gene_counter = Counter(gff_genes_dict.gene_list())
        if len(gff_genes_dict.gene_nuc_dict.keys()) == 0:
            print("Parse GFF file with create_nuc_dict set to True")
            sys.exit(0)
        self.gff_genes_dict = gff_genes_dict
        self.reads_per_gene_map_thr = threshold
    
    def reset(self):
        ''' Reset all to-run variables, as well as gene count and read count variables '''
        self.base_qual_run = False
        self.chr_lengths_run = False
        self.coverage_run = False
        self.map_qual_run = False
        self.read_length_hist_run = False
        self.reads_per_gene_run = False
        self.total_reads_run = False
        self.trim_reads_run = False
        
        self.base_qual_hist = Counter()
        self.base_qual_total_reads = 0
        self.chr_len = {}
        self.gene_counter = Counter()
        self.map_qual_hist = Counter()
        self.map_qual_low_reads = 0
        self.map_qual_total_reads = 0
        self.total_count = Counter()
        self.total_low_reads = 0
        self.total_low_reads_low_base_count = 0
        self.total_reads_count = 0
    
    def total_reads(self, threshold=20):
        ''' Calculates the total number of reads in the SAM file '''
        self.total_reads_run = True
        self.total_reads_count = 0
        self.map_qual_thr = threshold
    
    def trim_reads(self, threshold, o_fn="default"):
        ''' Trims a list of SAM files based on a minimum length threshold '''
        self.trim_reads_run = True
        if o_fn == "default":
            # Set default output filename
            dir_path, f_name_ext = os.path.split(self.f_in)
            f_name, f_ext = os.path.splitext(f_name_ext)
            f_outname = f_name + "_trimmed"
            o_fn = os.path.join(dir_path, f_outname, f_ext)
        self.trim_reads_o_fn = o_fn
    
    def start(self):
        ''' Gather all the requested information on the provided SAM file '''
        if self.chr_lengths_run == False and self.coverage_run == False and \
           self.read_length_hist_run == False and self.map_qual_run == False and \
           self.base_qual_run == False and self.reads_per_gene_run == False and \
           self.total_reads_run == False and self.trim_reads_run == False:
            print("Provide at least one function to be run on the SAM file.")
            sys.exit(0)
        
        # Open SAM file and start reading in lines
        gzipped = False
        past_header = False
        if self.f_in.endswith(".gz"):
            infile = gzip.open(self.f_in, 'rb')
            gzipped = True
        else:
            infile = open(self.f_in)
        for line in infile:
            if gzipped: line = line.decode('utf-8')
            line = line.split("\t")
            if past_header:
                # Looking at reads
                self._process_reads(line)
            else:
                # Looking at header information
                if line[0] == "@SQ" and (self.chr_lengths_run or self.coverage_run):
                    # Gather chromosome names and the lengths of those chromosomes
                    chromosome = line[1][3:]
                    length = int(line[2][3:])
                    self.chr_len[chromosome] = length
                elif line[0][0] != "@" and not past_header:
                    # We are no longer in the header
                    past_header = True
                    # Quit reading SAM file if we only want chromosome lengths
                    if self.coverage_run == False and \
                    self.read_length_hist_run == False and self.map_qual_run == False and \
                    self.base_qual_run == False and self.reads_per_gene_run == False and \
                    self.total_reads_run == False and self.trim_reads_run == False:
                        break
                    else:
                        # Make sure to process this read too
                        self._process_reads(line)
        infile.close()
        
        # Output relevant information to the screen or into files
        self._output()

def coverage_window(sam_file, win_list, chr_list):
    ''' Returns histogram of read coverage on a per-chromosome, per-bp-window basis '''
    
    chr_list = chr_list.split(",")
    
    # Initialize a multi-dimensional per-window, per-chromosome series of read count HMM paths
    win_db = {int(win): {} for win in win_list}
    for win in win_list:
        win_db[win] = {chromosome: Counter() for chromosome in chr_list}
    
    # Read in SAM file and create per-chromosome, reads-per-window histograms
    total_reads = 0
    if sam_file[-3:] == ".gz":
        infile = gzip.open(sam_file, 'rb')
    else:
        infile = open(sam_file)
    for line in infile:
        if sam_file[-3:] == ".gz": line = line.decode('utf-8')
        split_line = line.split("\t")
        chromosome = str(split_line[2])
        # Ignore read if it does not map (*) or if mapped to a chromosome not being considered
        if chromosome not in chr_list: continue
        quality = int(split_line[4])
        if quality < 20: continue # Ignore reads with a quality of 0, aka those which do not map uniquely, as well as low-quality mapping reads
        total_reads += 1
        spos = int(split_line[3])
        read_len = len(split_line[9])
        epos = spos + read_len
        for win in win_list:
            spos_for_dict = int(spos // win)
            epos_for_dict = int(epos // win)
            win_db[win][chromosome][spos_for_dict] += 1
            if spos_for_dict != epos_for_dict:
                total_reads += 1
                win_db[win][chromosome][epos_for_dict] += 1
    infile.close()
        
    # For all read-per-window counts that did not exist, assign them a value of 0
    for chromosome in chr_list:
        for win in win_list:
            for i in range(0,max(win_db[win][chromosome].keys())):
                if i not in win_db[win][chromosome].keys(): win_db[win][chromosome][i] = 0
    
    # Create multi-dimensional per-window, per-chromosome dictionaries for paths and histograms
    path_dict = {win: {} for win in win_list}
    hist_dict = {win: {} for win in win_list}
    
    # Fill in path and histogram dictionaries with data, then return them
    for win in win_list:
        for chromosome in chr_list:
            path = ' '.join([str(cov) for (win_num, cov) in sorted(win_db[win][chromosome].items(), key = lambda win_num: win_num)])
            path_dict[win][chromosome] = path
            hist = Counter({cov: 0 for cov in range(0,max(win_db[win][chromosome].values() ) + 1) } )
            for count in win_db[win][chromosome].values():
                hist[count] += 1
            hist_dict[win][chromosome] = hist
    return(hist_dict, path_dict)

def trim_reads(filelist, threshold=20):
    '''Trims a list of SAM files based on a minimum length threshold'''
    filelist = filelist.split(";")
    for filename in filelist:
        print("Trimming ",str(filename)," based on reads with a minimum length of ", str(threshold), "...", sep='')
        gzipped = False
        if filename.endswith(".gz"):
            outfilename = filename[:-7] + "_trimmed" + filename[-7:]
            infile = gzip.open(filename, 'rb')
            outfile = gzip.open(outfilename, 'wb')
            gzipped = True
        else:
            outfilename = filename[:-4] + "_trimmed" + filename[-4:]
            infile = open(filename, 'r')
            outfile = open(outfilename, 'w')
        for line in infile:
            if gzipped: line = line.decode('utf-8')
            splitlist = line.split()
            if len(splitlist) > 9:
                # Reads
                if len(splitlist[9]) >= threshold: # Nucleotides in read
                    if gzipped:
                        outfile.write(bytes(line,"UTF-8"))
                    else:
                        outfile.write(line)
            else:
                # Header
                if gzipped:
                    outfile.write(bytes(line,"UTF-8"))
                else:
                    outfile.write(line)
        infile.close()
        outfile.close()