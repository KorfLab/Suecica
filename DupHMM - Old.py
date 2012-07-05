#!/usr/bin/env python3.2
#from scipy.stats import poisson
import argparse, re, os, subprocess
from collections import OrderedDict
from scipy.stats import poisson
import sys; sys.path.append('../Packages/')
from GenomeTools import SAM, Viterbi

# DupHMM.py runs a duplication search HMM on a SAM file in an attempt to detect duplications. The SAM file is first broken up into user-specified window sizes,
#	  reads mapping per-window counted, and a Poisson regression is performed on the resulting count data. This Poisson lambda value is then used to generate
#	  a HMM parameter file, and the duplication search HMM is run using this parameter file. Alternatively, the user can provide a lambda value to be used
#	  and regression will be skipped.

def command_line():
	parser = argparse.ArgumentParser(description="DupHMM.py runs a duplication search HMM on a SAM file in an attempt to detect duplications. The SAM file is \
									 first broken up into user-specified window sizes, reads mapping per-window counted, and a Poisson regression is performed \
									 on the resulting count data. This Poisson lambda value is then used to generate a HMM parameter file, and the duplication \
									 search HMM is run using this parameter file. Alternatively, the user can provide a lambda value to be used and regression will be skipped.")
	#parser.add_argument('-sam', default="Data/weigel_col-0/SRR013327,SRR013328_7s_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/weigel_col-0/SRR013327,SRR013328_7s_aln.sam").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/weigel_bur-0/SRR013331,SRR013333_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/weigel_bur-0/SRR013331,SRR013333_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/weigel_tsu-1/SRR013335,SRR013337_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/weigel_tsu-1/SRR013335,SRR013337_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/Sue/Sue1/single_bp_preprocess/set1/090225_combined_sue_set1_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Sue/Sue1/single_bp_preprocess/set1/090225_combined_sue_set1_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/Sue/Sue1/single_bp_preprocess/set2/090601_sue_set2_PE_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Sue/Sue1/single_bp_preprocess/set2/090601_sue_set2_PE_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/Sue/Sue1/single_bp_preprocess/set2.5/sue_set2_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Sue/Sue1/single_bp_preprocess/set2.5/sue_set2_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/Sue/Sue1/single_bp_preprocess/set4/sue_set4_all_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Sue/Sue1/single_bp_preprocess/set4/sue_set4_all_aln.sam.gz").', metavar='SAMFilename')
	#parser.add_argument('-sam', default="Data/Sue/Sue1/single_bp_preprocess/set5/libSUE1_set5_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Sue/Sue1/single_bp_preprocess/set5/libSUE1_set5_aln.sam.gz").', metavar='SAMFilename')
	parser.add_argument('-sam', default="Data/C24/C24_all_reads_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/C24/C24_all_reads_aln.sam.gz").', metavar='SAMFilename')
	parser.add_argument('-chr', default="All", type=str, help='Chromosome to run HMM on (Default=HMM run on all chromosomes (separately) ). To look at >1 chromosome but not all chromosomes, provide them in a comma-delimited format: e.g. "Chr1,Chr2"', metavar='Chromosome')
	parser.add_argument('-w', default=5, type=str, help='Window size the HMM will use, in bps (Default=5). To look at multiple windows, provide them in a comma-delimited format: e.g. "1,5,10".', metavar='Window')
	parser.add_argument('-scn', default="1,2,5", type=str, help='Copy numbers for HMM states, separated by commas and listed in increasing numerical order (Default="1,2,5").', metavar='StateCopyNumbers')
	parser.add_argument('-sn', default="S,D,H", type=str, help='Names corresponding to copy numbers HMM states, separated by commas and listed in increasing numerical order (Default="S,D,H").', metavar='StateNames')
	parser.add_argument('-l', default=0, type=float, help='Lambda value for the poisson distribution representing a single copy region. (Default=Value determined by poisson regression on per-window histogram data)', metavar='LambdaMean')
	parser.add_argument('-trans', default=-100, type=float, help='Log probability for moving from same state->different state (Default=-100). Same state->Same state probability is unchangeable from 1. (Suggested:-10,-100,-200)', metavar='TransitionProb')
	parser.add_argument('-thr', default=0, type=str, help='Threshold value for read counts. If a window contains more reads than the threshold, the previous window\'s state and probability are retained (Default=No threshold). Multiple comma-separated values can be provided here, e.g. "10,20".', metavar='ReadThreshold')
	parser.add_argument('-thrp', default=0, type=float, help='Threshold value for read count per-window probability. Looking at the highest duplication state, the minimum per-window read count that has a Poisson probability which is below this threshold probability is determined. All per-window read counts with equal or higher counts for all states will retain the previous window\'s state and probability  (Default=No probability threshold)', metavar='ProbThreshold')
	#parser.add_argument('-o', default="Output/HMMCovWin_SRR013327,SRR013328/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_SRR013327,SRR013328/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_SRR013331,SRR013333/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_SRR013331,SRR013333")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_SRR013335,SRR013337/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_SRR013335,SRR013337/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_suecica_set1/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_suecica_set1/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_suecica_set2/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_suecica_set2/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_suecica_set2.5/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_suecica_set2.5/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_suecica_set4/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_suecica_set4/")', metavar='OutputDir')
	#parser.add_argument('-o', default="Output/HMMCovWin_suecica_set5/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_suecica_set5")', metavar='OutputDir')
	parser.add_argument('-o', default="Output/HMMCovWin_C24/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_C24")', metavar='OutputDir')
	
	args = parser.parse_args()
	sam = args.sam
	chr = args.chr
	win_list = args.w
	win_list = [entry for entry in win_list.split(',')]
	state_cns = args.scn
	state_names = args.sn
	pois_lambda = args.l
	trans_prob = args.trans
	t_count = (args.thr).split(",")
	t_prob = args.thrp
	outdir = args.o
	
	return(sam, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, t_prob, outdir)

def Dist_to_params(states, pois_lambda, trans_prob, t_count, t_prob, hist, win, chr, sam, outdir):
	# Calculating emission probabilities for X-fold copy states based on Poisson probabilities
	prob_dict = {}
	max_value = max(hist)
	if t_count != 0: # If there is a read count threshold t, discard read counts >= t
		max_value = int(t_count)
	
	j = 0
	for copy_num in states.keys():
		prob_list = []
		for i in range(0,max_value+1):
			prob = poisson.pmf(i, pois_lambda*copy_num)
			if prob == 0 and j == len(states.keys())-1: # Maximum probability for highest copy-number state
				left_tailed = True
				for entry in prob_list: # Check to make sure this isn't a left-tailed 0
					if entry != 0:
						left_tailed = False
						break
				if left_tailed == False: # Not a left-tailed 0, so stop calculating Poisson values
					max_value = i-1 # Don't consider count values with P=0 for lowest copy number state
					t_count = max_value
					break
			prob_list.append(prob)
		prob_dict[copy_num] = prob_list
		j += 1
		
	# If there is a probability threshold Pt, remove read counts that have a P <= Pt
	if t_prob != 0:
		for read_count in range(0, len(prob_dict[max(prob_dict.keys())])): # Looking at highest order state for per-window read count which is the closest to being <= t_prob
			if prob_dict[len(prob_dict)-1][read_count] <= t_prob:
				for i in states.keys():
					for j in range(len(prob_dict[list(states.keys())[i]])-1,read_count-1,-1):
						print("prob_dict[",i,"][",j,"]")
						del prob_dict[i][j]
				break
	
	outline = Dist_to_params_outline(states, trans_prob, prob_dict)
	t_outdir = os.path.join(outdir,str(win) + "bp-window/" + str(t_count) + "_read_threshold/")
	if not os.path.exists(t_outdir): os.makedirs(t_outdir)
	outfilename = str(t_outdir) + str(chr) + "_params.txt"
	with open(outfilename, 'w') as outfile:
		outfile.write(outline)
	return(max_value)

def Dist_to_params_outline(states, trans_prob, prob_dict):
	# States and start probabilities
	outline = "States:\n"
	outline += '\n'.join(states.values())
	outline += "\n\nStart Probabilities:\n"
	for i in states.values():
		outline += str(i) + "\t" + str(1/len(states)) + "\n"
	
	# Transition probabilities
	outline += "\nTransition Probabilities:"
	for i in range(0,len(states)):
		temp_trans = 1 if i == 0 else 10**trans_prob
		outline += "\n" + str(states[list(states.keys())[i]]) + "\t" + str(states[list(states.keys())[0]]) + "\t" + str(temp_trans)
		for j in range(1,len(states)):
			temp_trans = 1 if i == j else 10**trans_prob
			outline += "\n\t" + str(states[list(states.keys())[j]]) + "\t" + str(temp_trans)
	
	# Part of output containing emission probabilities
	outline += "\n\nEmission Probabilities:\n"
	for copy_num in states.keys():
		prob_list = prob_dict[copy_num]
		outline += str(states[copy_num])
		outline2 = ''.join(["\t" + str(i) + "\t" + str(prob_list[i]) + "\n" for i in range(0,len(prob_dict[copy_num]))])
		outline += str(outline2) + "\n"
	return outline

def HMM_Dup_Search(folder, chr, win, t_count, pois_lambda, states):
	outdir = str(folder) + str(win) + "bp-window/"
	param_file = str(outdir) + str(t_count) + "_read_threshold/" + str(chr) + "_params.txt"
	obs_file = str(outdir) + str(chr) + "_obs.txt"
	
	run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),"Viterbi.exe")
	params = ' '.join([str(run_script), "-p", str(param_file), "-o", str(obs_file)])
	HMMRun = subprocess.Popen(params, shell=False, stdout=subprocess.PIPE)
	path = HMMRun.communicate()[0].decode("utf-8")
	path = path[::-1]
	obs_list = []
	with open(obs_file) as infile:
		inline = infile.readlines()
		inline = inline[0].split(" ")
		obs_list = [int(inline[i]) for i in range(0,len(inline))]
	
	all_outlines = []
	for state in states:
		pattern = str(states[state]) + r"{1,}"
		recomp = re.compile(pattern)
		
		#outfilename = outdir + str(t_count) + "_read_threshold/" + str(chr) + "_" + str(state) + "XCopy.txt"
		#outlines = []
		for m in recomp.finditer(path):
			s_pos = m.start()*win
			e_pos = m.end()*win
			avg_cov = round(sum([obs_list[i] for i in range(m.start(),m.end()+1) if obs_list[i] < t_count]) / (m.end() - m.start()) / pois_lambda,2)
			#outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, avg_cov) )
			all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), avg_cov) )
		#if len(outlines) == 0: continue
		#outlines = sorted(outlines, key = lambda entry: entry[3], reverse = True)
		
		#with open(outfilename, 'w') as outfile:
		#	line = "Chr\tStart_pos\tEnd_pos\tDistance\n"
		#	outfile.write(line)
		#	for line in outlines:
		#		line = "\t".join([str(x) for x in line]) + "\n"
		#		#line = str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3])  + "\n"
		#		outfile.write(line)
	
	outfilename = outdir + str(t_count) + "_read_threshold/" + str(chr) + "_Results.txt"
	all_outlines = sorted(all_outlines, key = lambda entry: (entry[0], entry[1]))
	with open(outfilename, 'w') as outfile:
		line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tAvg_cov\n"
		outfile.write(line)
		for line in all_outlines:
			line = "\t".join([str(x) for x in line]) + "\n"
			outfile.write(line)

def run_HMM(sam, chr_len, chr, win, state_cns, state_names, pois_lambda, trans_prob, t_count, t_prob, outdir):
	# Create list of chromosomes to run HMM on - ignore A. arenosa chromosomes called "scaffoldX"
	pattern = "Chr\d{1,}"
	recomp = re.compile(pattern)
	chr_list = []
	if chr == "All":
		for chr_temp in chr_len:
			match = recomp.match(chr_temp)
			if match:
				chr_list.append(chr_temp)
	else:
		chr_list = [entry for entry in chr.split(',')]
	
	for chr in chr_list:
		# Create an HMM observation file (contains the number of reads per window in sequential order) and reads-per-window histogram table
		dirpath = os.path.dirname(os.path.abspath(__file__))
		i_folder = str(outdir) + str(win) + "bp-window/"
		i_fn_1 = str(outdir) + str(win) + "bp-window/" + str(chr) + "_obs.txt"
		hist = {}
		if not os.path.exists(i_fn_1): # Don't need to make the HMM observation file or histogram if they've already been created
			SAMFile = SAM.Parse(sam)
			hist = SAM_to_wincov(SAMFile, chr, win, outdir)
			if hist == 1:
				print(str(chr) + " does not exist in the provided SAM file")
				sys.exit()
		
		# Determines the poisson lambda value for a distribution of reads given the above reads-per-window histogram table
		# if the lambda value was not provided by the user
		singlecopy_lambda = pois_lambda # Variable for the 1X copy lambda. Here, it is kept separate from user-supplied lambda. 1X lambda needed for threshold calculations later
		if pois_lambda == 0: # If the user does not supply a lambda value, have R calculate the lambda that best fits the distribution
			o_fn = "Output/HMMInputfiles.txt"
			i_fn_2 = str(outdir) + str(win) + "bp-window/" + str(chr) + "_hist_distfits.txt"
			if not os.path.exists(i_fn_2):
				i_fn_3 = os.path.join(outdir, str(win) + "bp-window/" + str(chr) + "_hist.txt")
				line = str(i_fn_3) + "\n"
				with open(o_fn, 'w') as outfile:
					outfile.write(line)
				if "linux" in sys.platform:
					Rpath = r"/share/apps/R-2.14.2/bin/Rscript"
					params = ' '.join([str(Rpath) + " dist_fit.R", "--no-save"])
					simulation = subprocess.Popen(params, shell=True)
					simulation.wait()
				else:
					Rpath = r"C:/Program Files/R/R-2.15.0/bin/x64/Rscript.exe"
					params = ' '.join([str(Rpath) + " dist_fit.R", "--no-save --slave"])
					simulation = subprocess.Popen(params)
					simulation.wait()
			with open(i_fn_2) as infile:
				inlines = infile.readlines()
				inlines = [line.split() for line in inlines]
				singlecopy_lambda = float(inlines[1][1])
		
		# Output histogram
		if hist == {}: # Histogram file already exists, so histogram wasn't made from SAM file. Thus, read in histogram file.
			hist_i_fn = str(outdir) + str(win) + "bp-window/" + str(chr) + "_hist.txt"
			with open(hist_i_fn) as infile:
				inlines = infile.readlines()
				for line in inlines[1:]:
					cov, freq = line.split("\t")
					hist[int(cov)] = int(freq)
		
		temp_t_count = t_count[:]
		if t_count == [0]: # If no threshold was given, assign this to maximum coverage
			temp_t_count = [max(hist)]
		else: # Calculate read count threshold based on X-fold value of single-copy region lambda value
			for i in range(0,len(t_count)):
				if t_count[i][-1:] == "X":
					temp_t_count[i] = int(float(t_count[i][:-1]) * singlecopy_lambda)
		# States and their corresponding character codes
		state_cns_list = state_cns.split(',')
		state_names_list = state_names.split(',')
		states = OrderedDict({int(state_cns_list[i]): state_names_list[i] for i in range(0,len(state_cns_list))})
		
		# Cycle over each threshold value provided
		for threshold_count in temp_t_count:
			## Create an HMM parameter file (poisson probabilities of given # of reads in a window based on duplication status of that window)
			threshold_count = Dist_to_params(states, singlecopy_lambda, trans_prob, threshold_count, t_prob, hist, win, chr, sam, outdir)
			
			## Run HMM using generated parameter file and observation file
			HMM_Dup_Search(outdir, chr, win, threshold_count, singlecopy_lambda, states)

def SAM_to_wincov(SAMFile, chr, win, outdirpart):
	# Check to see if an observation file has already been made for the provided window
	outdir = os.path.join(outdirpart,str(win) + "bp-window/")
	if not os.path.exists(outdir): os.makedirs(outdir)

	hist, path = SAMFile.Coverage_Window(int(win),chr)
	if hist == 0 and path == 0: return(1) # Chromosome didn't exist in SAM file
	outfilename1 = os.path.join(outdir, str(chr) + "_obs.txt")
	with open(outfilename1, 'w') as outfile:
		outfile.write(path)
	
	outline = "\n".join([str(cov) + "\t" + str(freq) for cov, freq in sorted(hist.items(), key = lambda cov: cov[0])])
	outfilename2 = os.path.join(outdir, str(chr) + "_hist.txt")
	with open(outfilename2, 'w') as outfile:
		outline2 = "Coverage_per_" + str(win) + "bp_window\tCounts\n"
		outfile.write(outline2)
		outfile.write(outline)
	
	return(hist)

sam, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, t_prob, outdir = command_line()
chr_len = SAM.chr_lengths(sam)
for win in win_list:
	run_HMM(sam, chr_len, chr, int(win), state_cns, state_names, pois_lambda, trans_prob, t_count, t_prob, outdir)