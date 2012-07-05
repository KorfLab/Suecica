#!/usr/bin/env python
import re, time, threading, os, sys

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML, Record

from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
	
from scipy.stats import fisher_exact
import math
import sys

# Determines whether or not a read should map to A. thaliana or A. lyrata. This determination is made by:
#
# 1a) Looking at all SNPs in an Al read which is mapped to Al, and determining the number of SNPs which are present
# at a position that is indicative of a mismapping (a position with 0:100 Con/SNP distribution across all the original Al-mapping reads from BWA alignment).
# 1b) Looking at all SNPs in an Al read which is mapped to At, and determining the number of SNPs which are present
# at a position that is indicative of a SUN (a position with 50:50 Con/SNP distribution across all the original At-mapping reads from BWA alignment).
#
# If all of a read's SNPs are indicative of a SUN when mapped to At, or indicative of mismapping when mapped to Al,
# then the read shall be considered to be correctly mapped to A. thaliana.
#
# 1c) If there is no indication of correct mapping in At or mismapping in Al according to the above criteria, then a read
# will be determined to be correctly mapped to A. thaliana if the e-value for a read's alignment is lower with A. thaliana than with A. lyrata

class BLAST(threading.Thread):
	def __init__(self, reads, outfile, locks, genome):
		self.reads = reads
		self.outfile = outfile
		self.locks = locks
		self.genome = genome
		threading.Thread.__init__(self)
	
	def run(self):
		global finishcount
		global totalthreads
		finished = 'F'
		while True:
			if finished == 'T': break
			self.counter = -1
			for lock in self.locks:
				self.counter += 1
				if finished == 'T': break
				if not lock.locked():
					with lock:
						tempfile = ""
						if self.genome == "t":
							tempfile = "Output/AsAlReads/MappedToAt/AtReads_" + str(self.counter) + ".fasta"
						elif self.genome == "l":
							tempfile = "Output/AsAlReads/MappedToAl/AlReads_" + str(self.counter) + ".fasta"
						with open(tempfile, 'w') as outfile:
							outline = ''.join(self.reads)
							outfile.write(outline)
						
						if self.genome == "t":
							#self.blastn_cline = NcbiblastnCommandline(cmd='blastn', query=str(tempfile), db="at9db", evalue=0.0001, word_size=9, outfmt=5, out=str(self.outfile)) # For my computer (Windows)
							self.blastn_cline = NcbiblastnCommandline(cmd='/share/apps/blast+/bin/blastn', query=str(tempfile), db="../at9db", evalue=0.0001, word_size=9, outfmt=5, out=str(self.outfile)) # For Linux server
						elif self.genome == "l":
							#self.blastn_cline = NcbiblastnCommandline(cmd='blastn', query=str(tempfile), db="lyratadb", evalue=0.0001, word_size=9, outfmt=5, out=str(self.outfile)) # For my computer (Windows)
							self.blastn_cline = NcbiblastnCommandline(cmd='/share/apps/blast+/bin/blastn', query=str(tempfile), db="../lyratadb", evalue=0.0001, word_size=9, outfmt=5, out=str(self.outfile)) # For Linux server
						print(self.blastn_cline)
						stdout, stderr = self.blastn_cline()
						if len(stdout) > 0 or len(stderr) > 0: print(stdout,stderr,sep="\n")
						while True:
							try:
								os.remove(tempfile)
								break
							except:
								pass
						finished = 'T'
			time.sleep(1)
		finishcount += 1
		percentdone = round(finishcount / totalthreads * 100,1)
		print(str(percentdone) + "% done\n")

def BLAST_lyrata_reads(genome):
	# Align all A. lyrata-mapping reads to either A. thaliana (TAIR9) or A. lyrata using blastn, wordsize=9
	start = time.time()
	global finishcount
	global totalthreads
	pattern13 = r'^([ACTG]+)$'
	recomp13 = re.compile(pattern13)
	readsfilename = "Output/AsAlReads/AsAlReads.fasta"
	linelist = []
	with open(readsfilename) as infile:
		linelist = infile.readlines()
	
	num_threads = 8
	threadlist = []
	locks = []
	for i in range(0,num_threads):
		lock = threading.Lock()
		locks.append(lock)
	
	totalthreads = len(range(0,len(linelist),40000))
	
	for i in range(0,len(linelist),40000):
		nonambig_subset = []
		outfilename = ""
		if i + 40001 < len(linelist):
			if genome == "t":
				outfilename = "Output/AsAlReads/MappedToAt/AsAlReadsMappedToAt_" + str(int((i+40001)/2)-19999) + "-" + str(int((i+40001)/2)) + ".xml"
			elif genome == "l":
				outfilename = "Output/AsAlReads/MappedToAl/AsAlReadsMappedToAl_" + str(int((i+40001)/2)-19999) + "-" + str(int((i+40001)/2)) + ".xml"
			for entry in linelist[i:i+40000]:
				match = recomp13.match(entry)
				if entry[0] == ">" or match:
					nonambig_subset.append(entry)
		else:
			if genome == "t":
				outfilename = "Output/AsAlReads/MappedToAt/AsAlReadsMappedToAt_" + str(int((i+40001)/2)-19999) + "-" + str(int(len(linelist)/2)) + ".xml"
			elif genome == "l":
				outfilename = "Output/AsAlReads/MappedToAl/AsAlReadsMappedToAl_" + str(int((i+40001)/2)-19999) + "-" + str(int(len(linelist)/2)) + ".xml"
			for entry in linelist[i:]:
				match = recomp13.match(entry)
				if entry[0] == ">" or match:
					nonambig_subset.append(entry)
		p_thread = BLAST(nonambig_subset, outfilename, locks, genome)
		p_thread.start()
		threadlist.append(p_thread)

	for p_thread in threadlist:
		p_thread.join()
		
	total = round((time.time() - start) / 60,2)
	if genome == "t":
		print("All A. lyrata-mapping A. suecica reads have been mapped to A. thaliana (TAIR9)  (" + str(total) + " minutes)")
	elif genome == "l":
		print("All A. lyrata-mapping A. suecica reads have been mapped to A. lyrata (" + str(total) + " minutes)")

def compare_alignments():
	
	start = time.time()
	
	thal_genome, thal_cor_mapped = genome_info("t") # thal_genome = Dictionary containing thaliana genome. Key = Chr #. thal_cormapped = Dictionary containing positions with 50/50 Con/SNP ratio, indicative of correct mapping.
	lyr_genome, lyr_mismapped = genome_info("l") # lyr_genome = Dictionary containing lyrata genome. Key = Chr #. lyr_mismapped = Dictionary containing positions with 0/100 Con/SNP ratio, indicative of mismapping.
	
	pattern = r"Chr3"
	recomp = re.compile(pattern)
	pattern2 = r" "
	recomp2 = re.compile(pattern2)
	pattern3 = r"^SEQ_(\d+)_Chr\d{2}_\d+$"
	recomp3 = re.compile(pattern3)
	#pattern11 = r"-"
	#recomp11 = re.compile(pattern11)
	pattern17 = r"scaffold_(\d{1})"
	recomp17 = re.compile(pattern17)
	
	results = open("Output/AsAlReads/MappedToAt/AsAlReadsMappedToAt.xml")
	blast_records = NCBIXML.parse(results)
	
	chr3_reads = chr3_dup_reads = at_reads = at_dup_reads = 0
	next_record = False
	final_output = []
	at_read_db = {}
	
	for blast_record in blast_records: # Parse BLAST XML file for Al reads aligned to At to obtain mismapping information for reads mapped to At-Chr3 in the duplicated region
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				match = recomp.search(alignment.title)
				if match: # Read maps to At-Chr.3
					start_pos = hsp.sbjct_start
					end_pos = hsp.sbjct_end
					if min(start_pos, end_pos) > 3300050 and max(start_pos, end_pos) < 3400015:
						chr3_dup_reads += 1 # Keep count on how many Al reads map to within the duplicated region on At-Chr3
					else:
						chr3_reads += 1 # Keep count on how many Al reads map to At-Chr3, but not within duplicated region
					match = recomp3.match(blast_record.query)
					
					seq_num = int(match.group(1)) # Unique identifier for each A. lyrata-mapping read from the original BWA alignment of A. suecica to Thalyrata
					#query = Seq(recomp11.sub('',hsp.query), unambiguous_dna) # The read with gaps removed
					query = Seq(hsp.query, unambiguous_dna) # The read
					matching = hsp.match # String of matches (|) or mismatches ( ) to A. thaliana
					dbseq = Seq(hsp.sbjct, unambiguous_dna)
					evalue = hsp.expect
					mismatches = hsp.align_length - hsp.identities # Total number of mismatches when read is aligned to A. thaliana
					at_read_db[seq_num] = query, matching, dbseq, start_pos, end_pos, mismatches, evalue # Store relevant information for use when parsing the Al BLAST XML output
					
					next_record = True
					break
			if next_record == True:
				next_record = False
				break
	
	total = round((time.time() - start) / 60,2)
	print("\nFinished parsing the BLAST XML file containing A. lyrata-mapping reads mapped to A. thaliana, Chr. 3. Starting to parse the BLAST XML file containing A. lyrata-mapping reads mapped to A. lyrata and determining the best mapping for reads which have HSPs in both species (" + str(total) + " minutes)")
	results = open("Output/AsAlReads/MappedToAl/AsAlReadsMappedToAl.xml")
	blast_records = NCBIXML.parse(results)
	
	for blast_record in blast_records: # Parse BLAST XML file for Al reads aligned to Al to obtain mismapping information for reads mapped to Al. Will compare these mappings to At mappings to find the most likely correct mapping
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				match = recomp17.search(alignment.title)
				if match:
					al_chr = int(match.group(1))
					match = recomp3.match(blast_record.query)
					seq_num = int(match.group(1)) # Unique identifier for each A. lyrata-mapping read from the original BWA alignment of A. suecica to Thalyrata
					if seq_num in at_read_db: # If this sequence is found in this At dictionary, then it was found to be mapped to Chr. 3 in A. thaliana. Will now determine whether At or Al alignment is better
						al_mismatches = hsp.align_length - hsp.identities # Number of mismatches when read is aligned to A. lyrata genome
						if al_mismatches == 0: # Alignment to A. lyrata is perfect
							next_record = True
							break
						al_start_pos = hsp.sbjct_start
						al_end_pos = hsp.sbjct_end
						#al_query = Seq(recomp11.sub('',hsp.query), unambiguous_dna) # The read with gaps removed
						al_query = Seq(hsp.query, unambiguous_dna) # The read
						al_match = hsp.match
						al_dbseq = Seq(hsp.sbjct, unambiguous_dna)
						al_evalue = hsp.expect
						al_refseq = Seq(hsp.sbjct, unambiguous_dna)
						
						at_query, at_match, at_dbseq, at_start_pos, at_end_pos, at_mismatches, at_evalue = at_read_db[seq_num]
						at_good_SNPs = 0
						if at_start_pos > at_end_pos: # Query mapped to genome's reverse complement. To obtain correct genome nucleotide positions, must apply reverse complement again to get original genome sequence
							at_query = at_query.reverse_complement()
							at_match = at_match[::-1]
							at_dbseq = at_dbseq.reverse_complement()
							
						if al_start_pos > al_end_pos: # Query mapped to genome's reverse complement. To obtain correct genome nucleotide positions, must apply reverse complement again to get original genome sequence
							al_query = al_query.reverse_complement()
							al_match = al_match[::-1]
							al_dbseq = al_dbseq.reverse_complement()
						
						for m in recomp2.finditer(at_match):
							q_nuc = at_query[m.start()]
							db_pos = at_start_pos+m.start()
							if (3, db_pos) in thal_cor_mapped:
								if q_nuc == thal_cor_mapped[3, db_pos]:
									if at_start_pos > at_end_pos:#
										print("Reverse complement query (At):")#
										print(at_query,at_match,at_dbseq,at_start_pos,m.start(),q_nuc,thal_cor_mapped[3, db_pos],sep="\n")#
									at_good_SNPs += 1
						
						al_mismap_SNPs = 0
						for m in recomp2.finditer(al_match):
							q_nuc = al_query[m.start()]
							db_pos = al_start_pos+m.start()
							if (al_chr, db_pos) in lyr_mismapped:
								if q_nuc == lyr_mismapped[al_chr, db_pos]:
									if al_start_pos > al_end_pos:#
										print("Reverse complement query (Al):")#
										print(al_query,al_match,al_dbseq,al_start_pos,m.start(),q_nuc,lyr_mismapped[al_chr, db_pos],sep="\n")#
									al_mismap_SNPs += 1
						
						at_success = 0
						if at_mismatches > 0 and at_mismatches - at_good_SNPs == 0:
							at_success = 1
						elif al_mismatches > 0 and al_mismatches - al_mismap_SNPs == 0:
							at_success = 2
						elif at_evalue < al_evalue:
							at_success = 3
						if at_success != 0:
							if min(at_start_pos, at_end_pos) > 3300050 and max(at_start_pos, at_end_pos) < 3400015:
								at_dup_reads += 1
							else:
								at_reads += 1
							read_line = ">Query" + str(at_dup_reads + at_reads) + " (Read #" + str(seq_num) + ")\tChr3:" + str(at_start_pos) + "-" + str(at_end_pos)
							if at_success == 1:
								read_line = str(read_line) + "\t(" + str(at_good_SNPs) + " Good At SNPs"
							elif at_success == 2:
								read_line = str(read_line) + "\t(" + str(al_mismap_SNPs) + " Mismapped Al SNPs"
							else:
								read_line = str(read_line) + "\t(" + str(at_evalue) + "<" + str(al_evalue) + ""
							#at_query = recomp11.sub('',str(at_query))#
							read_line = str(read_line) + ")\n" + str(at_query) + "\n"
							final_output.append(read_line)
						
						next_record = True
						break
			if next_record == True:
				next_record = False
				break
	
	outline = ''.join(final_output)
	#with open("Output/AtReadsMismappedToAl-No_evalue_reads.fasta", 'w') as outfile:
	with open("Output/AtReadsMismappedToAl.fasta", 'w') as outfile:
		outfile.write(outline)
	print(chr3_dup_reads, "reads mapped to the duplicated region")
	print(at_dup_reads, "were more likely to map to the duplicated region")
	print(chr3_dup_reads - at_dup_reads, "were more likely to map to A. lyrata\n\n")
	print(chr3_reads, "reads mapped to Chr. 3, outside of the duplicated region")
	print(at_reads, "were more likely to map to Chr. 3, outside of the duplicated region")
	print(chr3_reads - at_reads, "were more likely to map to A. lyrata")
	
	total = round((time.time() - start) / 60,2)
	print("\nFinished comparing A. lyrata-mapping reads which map to within the duplicated region of At-Chr3 (" + str(total) + " minutes)")

def compile_BLAST_XMLs(genome):
	# Compiles the hundreds of individual BLAST XML results files into a single "AsAlReadsMappedToA(t|l).xml" file
	path = pattern5 = outfilename = ""
	if genome == "t":
		path = "Output/AsAlReads/MappedToAt/"
		pattern5 = r'AsAlReadsMappedToAt_\d{1,}-\d{1,}.xml'
		outfilename = "Output/AsAlReads/MappedToAt/AsAlReadsMappedToAt.xml"
	elif genome == "l":
		path = "Output/AsAlReads/MappedToAl/"
		pattern5 = r'AsAlReadsMappedToAl_\d{1,}-\d{1,}.xml'
		outfilename = "Output/AsAlReads/MappedToAl/AsAlReadsMappedToAl.xml"
	dirlist = os.listdir(path)
	filelist = []
	recomp5 = re.compile(pattern5)
	for filename in dirlist:
		match = recomp5.search(filename)
		if match:
			filelist.append(str(path) + str(filename))
	
	pattern6 = r'<BlastOutput_query-ID>Query_\d{1,}</BlastOutput_query-ID>'
	recomp6 = re.compile(pattern6)
	pattern7 = r'(\d+)'
	recomp7 = re.compile(pattern7)
	firstfile = True
	lastline = ""
	i = 1
	with open(outfilename, 'w') as outfile:
		for file in filelist:
			with open(file) as infile:
				linelist = infile.readlines()
				if firstfile == False:
					i += 1
					match = recomp6.search(linelist[7])
					if match:
						linelist[7] = recomp7.sub(str(i), linelist[7], 1) # Ex: <BlastOutput_query-ID>Query_1</BlastOutput_query-ID> changes to <BlastOutput_query-ID>Query_2</BlastOutput_query-ID>
					outfile.writelines(linelist[7:-1])
				else:
					outfile.writelines(linelist[:-1])
					lastline = linelist[-1] #Last line looks like this: </BlastOutput>
					firstfile = False
			print(str(round(i / len(filelist)*100,1)) + "% done with compilation of BLAST XMLs")
		outfile.write(lastline)
	
	for file in filelist:
		while True:
			try:
				os.remove(file)
				break
			except:
				pass
	
def genome_info(species):
	start = time.time()
	
	#Put all chromosome sequences into memory, as well as information on individual nucleotide positions with relevant SNP frequencies
	filename = pattern14 = recomp14 = pattern4 = recomp4 = ""
	if species == "t":
		filename ="../at9_genome.fasta"
		#AtChr3	22906118	T	12	............	EGEGEI>HHGHC	88FFFFFFFFFF
		pattern14 = r"^AtChr3\t(\d+)\t\w{1}\t\d+\t(.*)\t.*\t.*$"
		recomp14 = re.compile(pattern14)
		pattern4 = r"^>Chr(\d{1})"
		recomp4 = re.compile(pattern4)
	elif species == "l":
		filename = "../lyrata_8chr.fasta"
		#AlChr18	22906118	T	12	............	EGEGEI>HHGHC	88FFFFFFFFFF
		pattern14 = r"^AlChr1(\d{1})\t(\d+)\t\w{1}\t\d+\t(.*)\t.*\t.*$"
		recomp14 = re.compile(pattern14)
		pattern4 = r"^>scaffold_(\d{1})$"
		recomp4 = re.compile(pattern4)
	#-1A, +2C, or just a nucleotide without a preceding [+-]\d{1}
	pattern15 = r"([+-]\d{1})?([ATCG]{1})"
	recomp15 = re.compile(pattern15)
	pattern16 = r"[,.]"
	recomp16 = re.compile(pattern16)
	
	genome_list = []
	genome_dict = {} # Contains genome, with each key = chr. #, value = entire chr. sequence
	chr = 0
	with open(filename) as infile:
		for line in infile:
			match = recomp4.match(line)
			if match:
				chr = int(match.group(1))
				if (species == "l" and chr > 1) or (species == "t" and chr == 4):
					genome_dict[chr-1] = ''.join(genome_list)
				genome_list = []
			elif species == "l" or (species == "t" and chr == 3):
				genome_list.append(line.strip())
	if species == "l": genome_dict[chr] = ''.join(genome_list) # Necessary to add an entry for the last A. lyrata chromosome
	
	genome_pileup = {} # Will contain positions in genome which, across all reads, have approximately 0:100 Con/SNP ratio in A. lyrata or 50:50 in A. thaliana, indicative of a mismapping and correct mapping, respectively
	filename = "Data/SueAll/Sue-OurRESCAN/libsSue/ThalyrataAlignment/cleaned-RESCANlibsSue_sxn_pileup.txt"
	with open(filename) as infile:
		for line in infile:
			match = recomp14.match(line)
			if match:
				chr = pos = 0
				nucs = ""
				if species == "l":
					chr = int(match.group(1))
					pos = int(match.group(2))
					nucs = str(match.group(3)).upper()
					if len(nucs) >= 12:
						mismap = SNP_count(nucs, species, recomp15, recomp16)
						if mismap != "-":
							genome_pileup[(chr,pos)] = mismap
				elif species == "t":
					pos = int(match.group(1))
					nucs = str(match.group(2)).upper()
					if len(nucs) >= 12 and pos >= 3300050 and pos <= 3400015:
						cormap = SNP_count(nucs, species, recomp15, recomp16)
						if cormap != "-":
							genome_pileup[(3,pos)] = cormap
	
	total = round((time.time() - start) / 60,2)
	if species == "l":
		print("\nFinished compiling dictionary containing nucleotide positions which are indicative of mismapping (0:100 Con/SNP ratios) in A. lyrata (" + str(total) + " minutes)")
	elif species == "t":
		print("\nFinished compiling dictionary containing nucleotide positions which are indicative of correct mapping in a duplicated region (50:50 Con/SNP ratios) in A. thaliana (" + str(total) + " minutes)")
	return(genome_dict, genome_pileup)

def SNP_count(nucs, species, recomp15, recomp16):
	SNPs = {x: 0 for x in ('A','T','C','G')}
	maxSNP = ""
	maxSNP_value = 0
	countahead = 0
	indel = ""
	consensus_num = len(recomp16.findall(nucs))
	for m in recomp15.finditer(nucs):
		prev = str(m.group(1))
		SNP = str(m.group(2))
		if prev != "None": # nucleotide belongs to an indel
			countahead = int(prev[1]) - 1
			indel = str(prev) + str(SNP)
		elif countahead > 0: # continue adding nucleotides 1 at a time to indel
			indel += str(SNP)
			countahead -= 1
			if countahead == 0:
				if indel in SNPs:
					SNPs[indel] += 1
				else:
					SNPs[indel] = 1
		else: # SNP not part of an indel
			SNPs[SNP] += 1
	SNPcount = sum([v for v in SNPs.values()])
	maxSNP = max([k for k, v in SNPs.items()],key=lambda k: SNPs[k])
	maxSNP_value = SNPs[maxSNP]
	totalusedreads = consensus_num + maxSNP_value
	totalreads = consensus_num + SNPcount
	oddsratio = pvalue = nullvalue = 0
	mapping = "-"
	if species == "t":
		# Fisher's Exact Test for correctly mapped reads in a duplicated region on A. thaliana. Testing for 50/50 Con:SNP distribution.
		#
		#				Actual		Ideal (for 50/50)
		#	Consensus	  a			  b
		#		SNPs	  c			  d
		if consensus_num == totalreads and maxSNP_value == 0:
			return(mapping) # Time-saver
		else:
			oddsratio, pvalue = fisher_exact([[consensus_num, math.ceil(totalreads*0.5)], [maxSNP_value, math.ceil(totalreads*0.5)]])
	elif species == "l":
		# Fisher's Exact Test for ped reads on A. lyrata. Testing for 0/100 Con:SNP distribution.
		#
		#				Actual		Ideal (for 0/100)
		#	Consensus	  a			  0
		#		SNPs	  c			  d
		if consensus_num == 0: # Scipy's implementation of Fisher's Exact Test will error if this is 0; a variable in scipy becomes infinity, and it isn't handled properly. Therefore, increment all variables in the test by 1
			consensus_num += 1
			maxSNP_value += 1
			totalreads += 1
			nullvalue += 1
		if consensus_num == totalreads and maxSNP_value == 0:
			return(mapping) # Time-saver
		else:
			oddsratio, pvalue = fisher_exact([[consensus_num, nullvalue], [maxSNP_value, totalreads]])
	if pvalue < 0: pvalue = 0 # An error sometimes occurs for ratios enormously different from the ideal where p-value is slightly negative. Probably a floating point error
	if pvalue >= 0.10:
		# To a 90% confidence level, cannot reject the hypothesis that the observed SNP frequency is the same as 0/100 (A. lyrata) or 50/50 (A. thaliana) consensus/SNP
		mapping = maxSNP
	return(mapping)

def thaliana_BLAST_speedtest():
	# Used just for testing purposes to see how quickly BLAST will take given x number of query sequences in a batch.
	# Time/memory required for BLAST to finish increases quadratically with x. 20,000 sequences per batch appears to be
	# the highest reasonable amount achievable while the time taken to finish is still approximately linear.
	# Hence, parallel BLASTing of small 20,000 sequence batches is used in this script.
	filename = "Output/AlReadsAFew.fasta" # Contains 100,000 or so sequences
	filename2 = "Output/AlReadsSpeedTest.fasta" # Output containing time taken to BLAST x number of sequences
	linelist = []
	times = []
	with open(filename) as infile:
		linelist = infile.readlines()
	maxtime = 0
	maxlines = 0
	while maxtime < 3600: # Continue BLASTing increasingly larger batches until the time a BLAST batch takes exceeds an hour
		#maxlines += 1000
		maxlines += 20000
		with open(filename2, 'w') as outfile:
			outline = ''.join(linelist[0:maxlines-1])
			outfile.write(outline)
		start = time.time()
		blastn_cline = NcbiblastnCommandline(query=filename2, db="at9db", evalue=0.0001, outfmt=5, out="Output/AsAlReadsMappedToAt-Test.xml")
		stdout, stderr = blastn_cline()
		print(stdout,stderr,sep="\n")
		total = time.time()
		maxtime = total - start
		outstr = str(maxlines) + ": " + str(total) + "\n"
		times.append(outstr)
		print(maxlines,"done")
	print(times)
	sys.exit(0)

#thaliana_BLAST_speedtest() # Used just for testing purposes. See function for details.

#finishcount = 0
#totalthreads = 0
#BLAST_lyrata_reads("t") # Blast A. lyrata-mapping sequences against A. thaliana (t) or A. lyrata (l; done for more accurate alignment w/ blastn wordsize=9). Outputs into a bunch of XML files.
#compile_BLAST_XMLs("t") # Compiles the hundreds of individual BLAST XML results files into a single "AsAlReadsMappedToA(t|l).xml" file
#BLAST_lyrata_reads("l") # Blast A. lyrata-mapping sequences against A. thaliana (t) or A. lyrata (l; done for more accurate alignment w/ blastn wordsize=9). Outputs into a bunch of XML files.
#compile_BLAST_XMLs("l") # Compiles the hundreds of individual BLAST XML results files into a single "AsAlReadsMappedToA(t|l).xml" file

compare_alignments() # Opens AsAlReadsMappedToAt.xml and AsAlReadsMappedToAl.xml, parses the BLAST alignments, and compares each read which maps to At-Chr.3 to its best-mapped counterpart in Al. See function for details.