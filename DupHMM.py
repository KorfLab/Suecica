#!/usr/bin/env python3.2
import argparse, re, os, subprocess
from collections import OrderedDict, Counter
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
	parser.add_argument('-sam', default="Data/Ler-0/ERR031544,SRR279136_TAIR9_aln.sam.gz", type=str, help='SAM filename containing reads for the HMM to run on (Default="Data/Ler-0/ERR031544,SRR279136_TAIR9_aln.sam.gz").', metavar='SAMFilename')
	parser.add_argument('-csam', default="", type=str, help='SAM filename for an optional control data set. Duplications for a window will be called in the experimental SAM file on the basis of having two times (% of total reads in window) relative to the control. Default=No control.', metavar='ControlSAMFilename')
	parser.add_argument('-gff', default="../TAIR9_GFF3_genes_transposons.gff", type=str, help='GFF file for the reference genome. Used to automatically annotate regions containing transposons to make spotting false positives easier (Default="../TAIR9_GFF3_genes_transposons.gff").', metavar='GFFFile')
	parser.add_argument('-chr', default="All", type=str, help='Chromosome to run HMM on (Default=Run HMM on all chromosomes). To look at >1 chromosome but not all chromosomes, provide them in a comma-delimited format: e.g. "Chr1,Chr2"', metavar='Chromosome')
	parser.add_argument('-w', default=2500, type=str, help='Window size the HMM will use, in bps (Default=2500). To look at multiple windows, provide them in a comma-delimited format: e.g. "1000,2000,2500,3000".', metavar='Window')
	parser.add_argument('-scn', default="1,2,5", type=str, help='Copy numbers for HMM states, separated by commas and listed in increasing numerical order (Default="1,2,5").', metavar='StateCopyNumbers')
	parser.add_argument('-sn', default="S,D,H", type=str, help='Single character names corresponding to copy numbers HMM states, separated by commas and listed in increasing numerical order (Default="S,D,H", short for Single,Double,High).', metavar='StateNames')
	parser.add_argument('-l', default=0, type=float, help='Manually provide a lambda value for the poisson distribution representing a single copy region. (Default=Value determined by poisson regression on per-window histogram data)', metavar='LambdaMean')
	parser.add_argument('-trans', default=-50, type=float, help='Log probability for moving from same state->different state (Default=-50, Suggested:-50 or -100). Same State->Same State probability is unchangeable from 1.', metavar='TransitionProb')
	parser.add_argument('-thr', default=0, type=str, help='Threshold value for read counts. If a window contains more reads than the threshold, the previous window\'s state and probability are retained (Default=No threshold). Multiple comma-separated values can be provided here, e.g. "100,200". Threshold can also be set relative to the lambda value for a 1X region, e.g. "2X,10X,20X"', metavar='ReadThreshold')
	parser.add_argument('-o', default="Output/HMMCovWin_Ler-0/", type=str, help='Output directory for HMM results (Default="Output/HMMCovWin_Ler-0")', metavar='OutputDir')
	
	args = parser.parse_args()
	sam = args.sam
	control = args.csam
	gff = args.gff
	chr = args.chr
	win_list = args.w
	win_list = [int(entry) for entry in win_list.split(',')]
	state_cns = args.scn
	state_names = args.sn
	pois_lambda = args.l
	trans_prob = args.trans
	t_count = (args.thr).split(",")
	outdir = args.o
	
	return(sam, control, gff, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir)

def Dist_to_params(states, pois_lambda, trans_prob, t_count, hist, win, chr, sam, outdir):
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
					max_value = t_count = i-1 # Don't consider count values with P=0 for lowest copy number state
					break
			prob_list.append(prob)
		prob_dict[copy_num] = prob_list
		j += 1
	
	# Trim histograms such that the copy number states do not have an excessive number of right-tailed 0's
	# HMM treats missing reads-per-window counts as 0, so right-sided trimming histograms to the maximum read-per-window read count
	# with a non-zero probability for the highest copy number state helps save parameter file disk space
	i = 0
	for copy_num in states.keys():
		prob_dict[copy_num] = {i: prob_dict[copy_num][i] for i in range(0,t_count+1)}
		i += 1
		if i == len(states.keys()) - 1: break # Do not need to re-trim highest copy-number state
		
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

def Grab_Transposons(gff, chr_len):
	pattern = r"(transposable_element_gene|transposable_element)"
	recomp = re.compile(pattern)
	pattern2 = r"^ID=(\S+?);"
	#pattern2 = r"^ID=(AT\S{1}G\d{5});"
	recomp2 = re.compile(pattern2)
	tp_positions = {}
	for chr in chr_len:
		tp_positions[chr] = []
	with open(gff) as infile:
		for line in infile:
			line = line.split("\t")
			match = recomp.match(line[2])
			match2 = recomp2.match(line[8])
			if match and match2:
				chr = str(line[0])
				name = match2.group(1)
				spos = int(line[3])
				epos = int(line[4])
				type = 2 # Neither transposable_element nor transposable_element_gene
				if match.group(1) == "transposable_element_gene": type = 1
				else: type = 0
				tp_positions[chr].append((name, spos,epos, type))
	return tp_positions

def HMM_Dup_Search(folder, tp_positions, chr, win, t_count, pois_lambda, states, path_dict_exp_control_ratio={}):
	# Run HMM using C++ Viterbi algorithm. Wait for completion and grab HMM path from output, then parse and output results.
	outdir = str(folder) + str(win) + "bp-window/"
	param_file = str(outdir) + str(t_count) + "_read_threshold/" + str(chr) + "_params.txt"
	obs_file = str(outdir) + str(chr) + "_obs.txt"
	
	run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),"Viterbi.exe")
	params = ' '.join([str(run_script), "-p", str(param_file), "-o", str(obs_file)])
	HMMRun = subprocess.Popen(params, shell=False, stdout=subprocess.PIPE)
	path = HMMRun.communicate()[0].decode("utf-8") # Grab path from Viterbi.exe output
	path = path[::-1]
	
	# Open up observation file and grab observations for average fold coverage calculations across multiple windows
	obs_list = []
	with open(obs_file) as infile:
		inline = infile.readlines()
		inline = inline[0].split(" ")
		obs_list = [int(inline[i]) for i in range(0,len(inline))]
	
	# Figure out which windows contain transposons
	tp_windows = {}
	for (name, spos, epos, type) in tp_positions[chr]:
		spos_for_dict = int(spos // win)
		epos_for_dict = int(epos // win)
		try:
			tp_windows[spos_for_dict].append((name, type))
		except KeyError:
			tp_windows[spos_for_dict] = [(name, type)]
		if spos_for_dict != epos_for_dict:
			try:
				tp_windows[epos_for_dict].append((name,type))
			except KeyError:
				tp_windows[epos_for_dict] = [(name,type)]
	
	# For all windows that did not contain a transposon, assign them an empty list
	for i in range(0,len(obs_list)):
		if i not in tp_windows: tp_windows[i] = []
	
	# Write out HMM path along with average fold coverage and transposon count for each window
	# If using a control, include exp:control ratio in output as well
	all_outlines = []
	for state in states:
		pattern = str(states[state]) + r"{1,}"
		recomp = re.compile(pattern)
		
		for m in recomp.finditer(path):
			s_pos = m.start()*win
			e_pos = m.end()*win
			fold_cov = round(sum([obs_list[i] for i in range(m.start(),m.end()+1) if obs_list[i] < t_count]) / (m.end() - m.start()) / pois_lambda, 2)
			tp_union_list = set()
			for i in range(m.start(),m.end()+1):
				for (transposon_name, type) in tp_windows[i]:
					tp_union_list.add((transposon_name, type))
			if len(path_dict_exp_control_ratio.keys()) == 0: # Not using a control
				transposable_elements = len([name for name, type in tp_union_list if type == 0])
				transposable_element_genes = len([name for name, type in tp_union_list if type == 1])
				all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, transposable_element_genes, transposable_elements) )
			else: # Using a control
				exp_control_ratio_avg = round(sum([float(path_dict_exp_control_ratio[win][chr][i]) for i in range(m.start(),m.end()+1)]) / (m.end() - m.start()), 2)
				all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, transposable_element_genes, transposable_elements, exp_control_ratio_avg) )
	
	outfilename = outdir + str(t_count) + "_read_threshold/" + str(chr) + "_Results.txt"
	all_outlines = sorted(all_outlines, key = lambda entry: (entry[0], entry[1]))
	with open(outfilename, 'w') as outfile:
		if len(path_dict_exp_control_ratio.keys()) == 0: # Not using a control
			line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tTE Genes\tTEs\n"
		else:
			line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tTE Genes\tTEs\tExp:Control\n"
		outfile.write(line)
		for line in all_outlines:
			line = "\t".join([str(x) for x in line]) + "\n"
			outfile.write(line)
	print(outfilename, "complete")

def prepare_HMM(sam, control, tp_positions, chr_len, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir):
	# Create list of chromosomes to run HMM on
	pattern = "\S*Chr\d{1,}"
	recomp = re.compile(pattern)
	chr_list = []
	if chr == "All":
		for chr_temp in chr_len:
			match = recomp.match(chr_temp)
			if match:
				chr_list.append(chr_temp)
	else:
		chr_list = [entry for entry in chr.split(',')]
	
	# Create HMM observation files (contains # of reads per window in sequential order) and reads-per-window histogram tables
	hist_dict, path_dict = SAM_to_wincov(sam, win_list, chr_list, outdir, "False")
	
	# Determine the poisson lambda values for a distribution of reads given the aforementioned reads-per-window histogram tables.
	# If the user supplied a 1X copy region lambda value, use that instead.
	singlecopy_lambda = {win: {} for win in win_list}
	for win in win_list:
		for chr in chr_list:
			singlecopy_lambda[win][chr] = pois_lambda # Dictionary holding 1X copy lambdas, kept separate from user-supplied value for threshold calculations later
	if pois_lambda == 0: # If user didn't supply 1X lambda, calculate it given histogram
		for win in win_list:
			for chr in chr_list:
				# First create text file for R script containing the filename for the read histogram
				o_fn = "Output/HMMInputfiles.txt"
				i_fn_1 = str(outdir) + str(win) + "bp-window/" + str(chr) + "_hist_distfits.txt"
				if not os.path.exists(i_fn_1):
					i_fn_2 = os.path.join(outdir, str(win) + "bp-window/" + str(chr) + "_hist.txt")
					line = str(i_fn_2) + "\n"
					with open(o_fn, 'w') as outfile:
						outfile.write(line)
					# Now call R script "dist_fit.R" to perform Poisson regression
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
				# Read in "dist_fit.R" results file and grab Poisson regression lambda value
				with open(i_fn_1) as infile:
					inlines = infile.readlines()
					inlines = [line.split() for line in inlines]
					singlecopy_lambda[win][chr] = float(inlines[1][1])
	
	# Read in histograms if hist_dict is empty for a particular window (which happens when histograms were previously generated for all chromosomes for that window)
	# This is just a time-saving feature when performing re-runs on the same data set.
	for win in win_list:
		if not hist_dict[win]: # If this dictionary entry is empty, grab histogram from file
			for chr in chr_list:
				hist_i_fn = str(outdir) + str(win) + "bp-window/" + str(chr) + "_hist.txt"
				with open(hist_i_fn) as infile:
					inlines = infile.readlines()
					hist = Counter()
					for line in inlines[1:]:
						cov, freq = line.split("\t")
						hist[int(cov)] = int(freq)
				hist_dict[win][chr] = hist
	
	# Assign read count thresholds for window/chromosome combinations
	temp_t_count = {win: {} for win in win_list}
	for win in win_list:
		for chr in chr_list:
			temp_t_count[win][chr] = t_count[:]
			if t_count == [0]: # If no threshold was given, assign this to maximum coverage
				temp_t_count = [max(hist_dict[win][chr])]
			else: # Calculate read count threshold based on X-fold value of single-copy region lambda value
				for i in range(0,len(t_count)):
					if t_count[i][-1:] == "X":
						temp_t_count[win][chr][i] = int(float(t_count[i][:-1]) * singlecopy_lambda[win][chr])
	
	# Assign states and their corresponding character codes
	state_cns_list = state_cns.split(',')
	state_names_list = state_names.split(',')
	states = OrderedDict({int(state_cns_list[i]): state_names_list[i] for i in range(0,len(state_cns_list))})

	# If using a control, divide experimental % total read counts by control, store in path_dict_exp
	path_dict_exp = {}
	path_dict_control = {}
	if control != "":
		# Find % total read counts for experimental without re-reading the SAM file
		path_dict_exp = {win: {} for win in win_list}
		
		# Read in paths if path_dict is empty for a particular window (which happens when paths were previously generated for all chromosomes for that window)
		# This is just a time-saving used to calculate total number of reads mapping without re-reading SAM file
		for win in win_list:
			if path_dict[win] == {}: # If this dictionary entry is empty
				for chr in chr_list:
					path_i_fn = str(outdir) + str(win) + "bp-window/" + str(chr) + "_obs.txt"
					with open(path_i_fn) as infile:
						path = infile.readline()
						path = path.split(' ')
						path_dict[win][chr] = path
		total_reads = 0
		for chr in chr_list:
			for win in win_list:
				total_reads += sum([int(i) for i in path_dict[win][chr]])
		
		for chr in chr_list:
			for win in win_list:
				pct_list = [(int(i)+1)/total_reads for i in path_dict[win][chr]]
				path_dict_exp[win][chr] = pct_list
		
		hist_dict_control, path_dict_control = SAM_to_wincov(control, win_list, chr_list, outdir, "True")
		for win in win_list:
			for chr in chr_list:
				path_dict_control[win][chr] = path_dict_control[win][chr].split(' ')
				path_dict_exp[win][chr] = {i: float(path_dict_exp[win][chr][i]) / float(path_dict_control[win][chr][i]) for i in range(0,len(path_dict_exp[win][chr]))}
		path_dict_control = {}
	
	# Cycle over each window, chromosome, and threshold value combination provided and run HMM for each combination
	for win in win_list:
		for chr in chr_list:
			for threshold_count in temp_t_count[win][chr]:
				## Create an HMM parameter file (poisson probabilities of given # of reads in a window based on duplication status of that window)
				threshold_count = Dist_to_params(states, singlecopy_lambda[win][chr], trans_prob, threshold_count, hist_dict[win][chr], win, chr, sam, outdir)
				
				## Run HMM using generated parameter file and observation file
				HMM_Dup_Search(outdir, tp_positions, chr, win, threshold_count, singlecopy_lambda[win][chr], states, path_dict_exp)

def SAM_to_wincov(sam_file, win_list, chr_list, outdirpart, control_use):
	# Check to see if, for a given window size, all chromosomes have already been examined. If so, do not re-create data files for that window size
	win_list_needed = win_list
	if control_use == "False":
		win_list_needed = []
		for win in win_list:
			not_needed = 0
			for chr in chr_list:
				i_fn = str(outdirpart) + str(win) + "bp-window/" + str(chr) + "_obs.txt"
				if os.path.exists(i_fn): not_needed += 1
			if not_needed != len(chr_list): win_list_needed.append(win)
	
	hist_dict = {win: {} for win in win_list}
	path_dict = {win: {} for win in win_list}
	chr_string = ';'.join(chr_list)
	if win_list_needed: # If list is not empty, create paths for each chromosome/window pair
		# Get read count histograms and paths for the given window sizes
		hist_dict, path_dict = SAM.Coverage_Window(sam_file, win_list_needed, chr_string, control_use)
		if control_use == "True":
			pass # Only need path_dict for a file being used only to find the exp:control ratio of % total reads
		elif control_use == "False":
			# Check to see if folders have already been made for the provided windows
			for chr in chr_list:
				for win in win_list:
					outdir = os.path.join(outdirpart,str(win) + "bp-window/")
					if not os.path.exists(outdir): os.makedirs(outdir)
			
			for win in win_list:
				outdir = os.path.join(outdirpart,str(win) + "bp-window/")
				for chr in chr_list:
					outfilename1 = os.path.join(outdir, str(chr) + "_obs.txt")
					outfilename2 = os.path.join(outdir, str(chr) + "_hist.txt")
					with open(outfilename1, 'w') as outfile:
						outfile.write(path_dict[win][chr])
					
					outline = "\n".join([str(cov) + "\t" + str(freq) for cov, freq in sorted(hist_dict[win][chr].items(), key = lambda cov: cov[0])])
					with open(outfilename2, 'w') as outfile:
						outline2 = "Coverage_per_" + str(win) + "bp_window\tCounts\n"
						outfile.write(outline2)
						outfile.write(outline)
	
	return(hist_dict, path_dict)

sam, control, gff, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir = command_line()
chr_len = SAM.Chr_Lengths(sam)
tp_positions = Grab_Transposons(gff, chr_len)
prepare_HMM(sam, control, tp_positions, chr_len, chr, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir)