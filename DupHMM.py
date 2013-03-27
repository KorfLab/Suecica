#!/usr/bin/env python3
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import argparse, re, random, os, subprocess, sys
from collections import Counter, OrderedDict
from copy import deepcopy
from math import exp, factorial, log, lgamma, fabs
import SAM
verbose = False

# DupHMM.py runs a duplication-finding HMM on a SAM file. The SAM file is first broken up into user-specified window sizes,
#	  reads mapping per-window counted, and a Poisson regression is performed on the resulting histogram of count data. This Poisson lambda
#	  value is then used to generate a HMM parameter file, which the duplication searching HMM then uses to identify duplicated windows.

def cluster_group_assign_data_points(cluster_group, hist_dict, max_reads):
	# This function groups data points to the appropriate cluster point and returns the data point-cluster point distance sum
	
	data_count_dict = Counter({cluster: 0 for cluster in cluster_group})
	optimum_cluster_dists = {cluster: [] for cluster in cluster_group} # Value = List containing (data point-cluster distance**2 * freq)
	for data_point, freq in hist_dict.items():
		# Cycle through each data point, finding which cluster point each is closest to
		# and storing (distance * freq) to that cluster's list of data points
		min_cluster = -1												# Will hold closest cluster to data point. -1 is a default 'not assigned yet' value
		min_cluster_dist = (int(fabs(data_point - cluster_group[0]))+2)**2+1	# Will hold closest cluster's distance from data point. Default is larger than 1st entry
		for cluster_point in cluster_group:
			data_cluster_dist = (int(fabs(data_point - cluster_point))+2)**2
			if data_cluster_dist < min_cluster_dist:
				min_cluster_dist = data_cluster_dist
				min_cluster = cluster_point
		optimum_cluster_dists[min_cluster].append(min_cluster_dist * freq)
		data_count_dict[min_cluster] += freq
	
	# Sum up all the data point-cluster point distances for each cluster.
	# This sum will be used as a 'fitness' value for this group of clusters.
	dist_sum = 0
	for cluster_point in cluster_group:
		dist_sum += sum( [min_cluster_dist for min_cluster_dist in optimum_cluster_dists[cluster_point]] )
	return(dist_sum, data_count_dict)

def cluster_mutate(cluster_group, small_mut_freq, large_mut_freq, max_reads):
	# With low frequency , shift cluster points a bit
	# Maximum small mutation shift distance is 0.5% of maximum number of reads per window possible
	# Maximum large mutation shift distance is 1% of maximum number of reads per window possible
	small_mut_dist = 0.005
	large_mut_dist = 0.01
	
	new_cluster_group = cluster_group
	for i in range(len(cluster_group)):
		random_value = random.random()
		shift_dist = 0
		if random_value < large_mut_freq:
			shift_dist = int(random.randrange(int(2 * max_reads*large_mut_dist)) - max_reads*large_mut_dist)
		elif random_value < small_mut_freq:
			shift_dist = int(random.randrange(int(2 * max_reads*small_mut_dist)) - max_reads*small_mut_dist)
		if shift_dist > 0:
			if new_cluster_group[i] + shift_dist <= max_reads:
				new_cluster_group[i] += shift_dist
			else:
				new_cluster_group[i] = max_reads
		elif shift_dist < 0:
			if new_cluster_group[i] + shift_dist > 0:
				new_cluster_group[i] += shift_dist
			else:
				new_cluster_group[i] = 0

	new_cluster_group = sorted(new_cluster_group)
	return new_cluster_group

def cluster_next_generation(cluster_group_list, hist_dict, breeding_individuals, total_individuals, small_mut_freq, large_mut_freq, k, max_reads):
	# Create a new generation of cluster points given a breeding population of clusters
	
	# From the previous generation's total_individuals groups of k # of clusters, pick the top cluster groups which
	# have the lowest data point-cluster point distance sum to populate the next generation
	next_gen_cluster_group_list = sorted(cluster_group_list, key = lambda cluster_group_list: cluster_group_list[1])[:breeding_individuals]
	
	# Change mutation frequency to be scaled to the number of clusters in a group. This is done since as k increases
	# it becomes increasingly more difficult to reach convergence with rapidly shifting cluster groups, typically becoming
	# quite difficult around k=9
	small_mut_freq /= k
	large_mut_freq /= k
	
	# Create (total_individuals-#_of_breeding_clusters) individuals. Their cluster points will result from 'crossovers' among the
	# breeding individuals. The next generation individuals are also subject to 'mutational' position drifting
	cluster_list = []
	for cluster_group, dist_sum, data_count_dict in next_gen_cluster_group_list:
		for cluster in cluster_group:
			cluster_list.append(cluster)
	cluster_set = set(cluster_list)
	
	for i in range(breeding_individuals, total_individuals):
		unique_indiv_created = False
		while unique_indiv_created == False:
			new_cluster_group = random.sample(cluster_set,k)
			new_cluster_mutated_group = cluster_mutate(new_cluster_group, small_mut_freq, large_mut_freq, max_reads)
			new_cluster_mutated_set = set(new_cluster_mutated_group) # Sets do not allow for duplicate clusters
			if len(new_cluster_mutated_set) == len(new_cluster_mutated_group):
				new_dist_sum, new_data_count_dict = cluster_group_assign_data_points(new_cluster_mutated_group, hist_dict, max_reads)
				next_gen_cluster_group_list.append([new_cluster_mutated_group, new_dist_sum, new_data_count_dict]) # Note: mutated group may or may not be mutated
				unique_indiv_created = True
	
	next_gen_cluster_group_list = sorted(next_gen_cluster_group_list, key = lambda next_gen_cluster_group_list: next_gen_cluster_group_list[1])
	
	return next_gen_cluster_group_list

def command_line():
	parser = argparse.ArgumentParser(description="DupHMM.py runs a duplication-finding HMM on a SAM file. The SAM file is first broken up into user-specified window sizes, \
									 reads mapping per-window counted, and a Poisson regression is performed on the resulting histogram of count data. This Poisson lambda \
									 value is used to generate a HMM parameter file, which the duplication searching HMM then uses.")
	parser.add_argument('-m', default=1, type=int, choices=[1,2], help="Mode in which DupHMM will run. DupHMM can attempt to find CNVs via two means: (1) Fitting a Poisson distribution to a chromosome & assuming that Poisson fit represents a single copy region from which other CNV states' Poisson fits are determined, and (2) Using a k-means clustering algorithm to identify all major CNVs. Method (1) is the default used, and should be used for normal chromosomes, while Method (2) is more appropriate for abnormal chromosomes, e.g. cancer data sets, shattered chromosomes", metavar="ProgramMode")
	parser.add_argument('-sam', default="", type=str, help='SAM filename containing reads for the HMM to run on.', metavar='SAMFilename')
	parser.add_argument('-o', default="", type=str, help='Output directory for HMM results', metavar='OutputDir')
	parser.add_argument('-gff', default="", type=str, help='GFF file for the reference genome. Used to annotate regions in the results which contain transposons. This makes spotting false positives easier.', metavar='GFFFile')
	parser.add_argument('-chr', default="All", type=str, help='Chromosome to run HMM on (Default=Run HMM on all chromosomes). To look at >1 chromosome but not all chromosomes, provide them in a comma-delimited format: e.g. "Chr1,Chr2"', metavar='Chromosome')
	parser.add_argument('-w', default="", type=str, help='Window size the HMM will use, in bps (Default=2500). To look at multiple windows, provide them in a comma-delimited format: e.g. "1000,2000,2500,3000".', metavar='Window')
	parser.add_argument('-scn', default="1,2,5", type=str, help='Mode 1 Feature: Copy numbers for HMM states, separated by commas and listed in increasing numerical order (Default="1,2,5").', metavar='StateCopyNum')
	parser.add_argument('-sn', default="", type=str, help='Mode 1 Feature (Optional): Single character names corresponding to copy number HMM states, separated by commas and listed in increasing numerical order. Example: If using 1,2,5 as CNV states, try "S,D,H", short for Single,Double,High. (Default="A,B,C...")', metavar='StateNames')
	parser.add_argument('-l', default=0, type=float, help='Mode 1 Feature: Manually provide a lambda value for the poisson distribution representing a single copy region. Entering this skips Poisson regression over the per-window read count histogram.)', metavar='LambdaMean')
	parser.add_argument('-trans', default=-50, type=float, help='Mode 1 Feature: Log probability for moving from same state->different state (Default=-50, Suggested:-50 or -100). Same State->Same State probability is unchangeable from 1.', metavar='TransProb')
	parser.add_argument('-thr', default="10X", type=str, help='Mode 1 Feature: Threshold value for read counts. If a window contains more reads than the threshold, the previous window\'s state and probability are copied to the read spike window, effectively ignoring read spikes. Multiple comma-separated values can be provided here for multiple re-runs, e.g. "100,200". Threshold can also be set relative to the lambda value for a 1X region, e.g. "10X,20X". Default="10X,20X"', metavar='ReadThreshold')
	parser.add_argument('-v', action='store_true', help='Turns on verbose output for k-means clustering')
	
	args = parser.parse_args()
	mode = args.m
	sam = args.sam
	gff = args.gff
	chromosome = args.chr
	win_list = args.w
	win_list = [int(entry) for entry in win_list.split(',')]
	state_cns = args.scn
	state_names = args.sn
	if state_names == "":
		state_names = ",".join([chr(x) for x in range(65,65+len(state_cns))])
	pois_lambda = args.l
	trans_prob = args.trans
	t_count = (args.thr).split(",")
	outdir = args.o
	
	return(mode, sam, gff, chromosome, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir)

def dist_to_params_mode_one(states, pois_lambda, trans_prob, t_count, hist, win, chromosome, sam, outdir):
	# Distribution_to_parameters_file
	# Converts per-window read count histogram with a Poisson fit to it into
	# parameters that an HMM can utilize
	
	# Calculating emission probabilities for X-fold copy states based on Poisson probabilities
	prob_dict = {}
	max_value = max(hist)
	if t_count != 0: # If there is a read count threshold t, discard read counts >= t
		max_value = int(t_count)
	
	j = 0
	for copy_num in states.keys():
		prob_list = []
		for k in range(0,max_value+1):
			mu = pois_lambda*copy_num
			#nonlogprob = exp(-mu) * mu**k / factorial(k)	# Fails when k is too large. Must be handled in log-space
			prob = exp(k*log(mu)-lgamma(k+1) - mu)			# Above equation handled in log-space
			if prob == 0 and j == len(states.keys())-1:
				# Find the maximum non-zero probability for highest copy-number state
				left_tailed = True
				for entry in prob_list: # Check to make sure this isn't a left-tailed 0
					if entry != 0:
						left_tailed = False
						break
				if left_tailed == False:
					# Not a left-tailed 0, so stop calculating Poisson values beyond this maximum non-zero probability
					# Want to next remove count values for lower copy number states where it is guaranteed that P=0
					max_value = t_count = k-1
					break
			prob_list.append(prob)
		prob_dict[copy_num] = prob_list
		j += 1
	
	# Trim histograms such that the copy number states do not have an excessive number of right-tailed 0's
	# HMM treats missing reads-per-window counts as 0, so right-sided trimming histograms to the maximum read-per-window read count
	# with a non-zero probability for the highest copy number state helps save parameter file disk space, time, and memory
	i = 0
	for copy_num in states.keys():
		prob_dict[copy_num] = {i: prob_dict[copy_num][i] for i in range(0,t_count+1)}
		# Do not need to re-trim highest copy-number state
		i += 1
		if i == len(states.keys()) - 1: break
		
	# Create parameter file
	# Call dist_to_params_outline_mode_one to create the line of text to be written out
	outline = dist_to_params_outline_mode_one(states, trans_prob, prob_dict)
	t_outdir = os.path.join(outdir,str(win) + "bp_window/" + str(t_count) + "_read_threshold/")
	if not os.path.exists(t_outdir): os.makedirs(t_outdir)
	outfilename = os.path.join(t_outdir, str(chromosome) + "_params.txt")
	with open(outfilename, 'w') as outfile:
		outfile.write(outline)
	return(max_value)

def dist_to_params_mode_two(hist, k_cluster_count, k_loc, k_data_count, win, chromosome, sam, outdir):
	# Distribution_to_parameters_file
	# Converts per-window read count histogram with k-clusters into
	# parameters that an HMM can utilize
	
	# Calculating emission probabilities for CNV states based on Poisson probabilities
	prob_dict = {}
	max_value = max(hist)

	j = 0
	for cluster_pos in k_loc:
		prob_list = []
		for k in range(0,max_value+1):
			mu = cluster_pos
			#nonlogprob = exp(-mu) * mu**k / factorial(k)	# Fails when k is too large. Must be handled in log-space
			prob = exp(k*log(mu)-lgamma(k+1) - mu)			# Above equation handled in log-space
			if prob == 0 and j == len(k_loc)-1:
				# Find the maximum non-zero probability for highest copy-number state
				left_tailed = True
				for entry in prob_list: # Check to make sure this isn't a left-tailed 0
					if entry != 0:
						left_tailed = False
						break
				if left_tailed == False:
					# Not a left-tailed 0, so stop calculating Poisson values beyond this maximum non-zero probability
					# Want to next remove count values for lower copy number states where it is guaranteed that P=0
					max_value = k-1
					break
			prob_list.append(prob)
		prob_dict[cluster_pos] = prob_list
		j += 1
	
	# Trim histograms such that the copy number states do not have an excessive number of right-tailed 0's
	# HMM treats missing reads-per-window counts as 0, so right-sided trimming histograms to the maximum read-per-window read count
	# with a non-zero probability for the highest copy number state helps save parameter file disk space, time, and memory
	i = 0
	for cluster_pos in k_loc:
		prob_dict[cluster_pos] = {i: prob_dict[cluster_pos][i] for i in range(0,max_value+1)}
		# Do not need to re-trim highest copy-number state
		i += 1
		if i == len(k_loc) - 1: break
	
	# Create transition probabilities between states based on relative percentages of each state
	total_data_points = sum([point_sum for point_sum in k_data_count.values()])
	trans_prob_dict = {cluster: {} for cluster in k_loc}
	for cluster in k_loc:
		trans_prob_dict[cluster] = {x: 10**(-75 - (-75 * k_data_count[x] / total_data_points)) for x in k_loc}
	
	# Create parameter file
	# Call dist_to_params_outline_mode_two to create the line of text to be written out
	outline = dist_to_params_outline_mode_two(k_loc, prob_dict, trans_prob_dict)
	t_outdir = os.path.join(outdir,str(win) + "bp_window/")
	if not os.path.exists(t_outdir): os.makedirs(t_outdir)
	outfilename = os.path.join(t_outdir, str(chromosome) + "_params.txt")
	with open(outfilename, 'w') as outfile:
		outfile.write(outline)

def dist_to_params_outline_mode_one(states, trans_prob, prob_dict):
	# Distribution_to_parameters_file_outline
	# Formats HMM parameters for output into a format that the HMM can read in
	
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

def dist_to_params_outline_mode_two(k_loc, prob_dict, trans_prob_dict):
	# Distribution_to_parameters_file_outline
	# Formats HMM parameters for output into a format that the HMM can read in
	
	# States and start probabilities
	# ASCII 65 = 'A'
	outline = "States:\n"
	states = [chr(x) for x in range(65,65+len(k_loc))]
	outline += '\n'.join(states)
	outline += "\n\nStart Probabilities:\n"
	for x in states:
		outline += str(x) + "\t" + str(1/len(states)) + "\n"
	
	# Transition probabilities
	outline += "\nTransition Probabilities:"
	for i in range(0,len(states)):
		outline += "\n" + str(states[i]) + "\t" + str(states[0]) + "\t" + str(trans_prob_dict[k_loc[i]] [k_loc[0]] )
		for j in range(1,len(states)):
			outline += "\n\t" + str(states[j]) + "\t" + str(trans_prob_dict[k_loc[i]] [k_loc[j]] )
	
	# Part of output containing emission probabilities
	outline += "\n\nEmission Probabilities:\n"
	for i in range(65, 65+len(k_loc)):
		cluster = k_loc[i-65]
		prob_list = prob_dict[cluster]
		outline += str(states[i-65])
		outline2 = ''.join(["\t" + str(i) + "\t" + str(prob_list[i]) + "\n" for i in range(0,len(prob_dict[cluster]))])
		outline += str(outline2) + "\n"
	return outline

def grab_transposons(gff, chr_len):
	pattern = r"(transposable_element_gene|transposable_element|transposon)"
	recomp = re.compile(pattern)
	pattern2 = r"^ID=(\S+?);"
	recomp2 = re.compile(pattern2)
	tp_positions = {}
	if gff != "":
		# If a GFF file was provided, fill tp_positions dictionary
		# with transposon positions
		for chromosome in chr_len:
			tp_positions[chromosome] = []
		with open(gff) as infile:
			for line in infile:
				line = line.split("\t")
				match = recomp.match(line[2])
				match2 = recomp2.match(line[8])
				if match and match2:
					chromosome = str(line[0])
					name = match2.group(1)
					spos = int(line[3])
					epos = int(line[4])
					type = 2 # Neither transposable_element nor transposable_element_gene
					if match.group(1) == "transposable_element_gene" or match.group(1) == "transposon": type = 1
					elif match.group(1) == "transposable_element": type = 0
					tp_positions[chromosome].append((name, spos,epos, type))
	return tp_positions

def hmm_dup_search(mode, folder, tp_positions, chromosome, win, pois_lambda, states, t_count):
	# HMM_Dup_search runs the HMM using a Viterbi algorithm written in C++ .
	# It waits for the HMM run to complete, grabs the HMM results from Viterbi's output, then parses and outputs results.
	
	outdir = os.path.join(folder, str(win) + "bp_window/")
	param_file = os.path.join(outdir, str(t_count) + "_read_threshold/" + str(chromosome) + "_params.txt")
	obs_file = os.path.join(outdir, str(chromosome) + "_obs.txt")
	path_file = os.path.join(outdir, str(chromosome) + "_path.txt")
	
	# Run Viterbi algorithm written in C++
	run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),"viterbi")
	params = ' '.join([str(run_script), "-p", str(param_file), "-o", str(obs_file)])
	HMMRun = subprocess.Popen(params, shell=True, stdout=subprocess.PIPE)
	path = HMMRun.communicate()[0].decode("utf-8") # Grab path from Viterbi output
	path = path[::-1]
	
	# Write path file
	with open(path_file, 'w') as outfile:
		outline = "\n".join(x for x in path)
		outfile.write(outline)
	
	# Open up observation file and grab observations for average fold coverage calculations across multiple windows
	obs_list = []
	with open(obs_file) as infile:
		inline = infile.readlines()
		inline = inline[0].split(" ")
		obs_list = [int(inline[i]) for i in range(0,len(inline))]
	
	# If a GFF file was provided, figure out which windows contain transposons
	tp_windows = {}
	if tp_positions != {}:
		for (name, spos, epos, type) in tp_positions[chromosome]:
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
	
	if mode == 2: # K-means clustering mode
		# Write cluster file
		cluster_file = os.path.join(outdir, str(chromosome) + "_clusters.txt")
		k_cluster_group_dict = deepcopy(states)
		# ASCII 65 = 'A'
		states = [chr(x) for x in range(65,65+len(k_cluster_group_dict))]
		with open(cluster_file, 'w') as outfile:
			outline = "HMM_State\tPosition\n"
			outfile.write(outline)
			for i in range(0, len(states)):
				outline = str(states[i]) + "\t" + str(k_cluster_group_dict[i]) + "\n"
				outfile.write(outline)
		
	# Write out HMM path along with average fold coverage and transposon count for each window
	all_outlines = []
	for state in states:
		if mode == 1:   # Poisson fit mode
			pattern = str(states[state]) + r"{1,}"
		elif mode == 2: # K-means clustering mode
			pattern = str(state) + "+"
		recomp = re.compile(pattern)
		
		for m in recomp.finditer(path):
			s_pos = m.start()*win
			e_pos = m.end()*win
			fold_cov = round(sum([obs_list[i] for i in range(m.start(),m.end()+1) if obs_list[i] < t_count]) / (m.end() - m.start()) / pois_lambda, 2)
			tp_union_list = set()
			# If a GFF file was provided, include per-window transposon counts in results
			if tp_windows != {}:
				for i in range(m.start(),m.end()+1):
					for (transposon_name, type) in tp_windows[i]:
						# Only counting the number of unique transposon types present in a window
						tp_union_list.add((transposon_name, type))
				transposable_elements = len([name for name, type in tp_union_list if type == 0])
				transposable_element_genes = len([name for name, type in tp_union_list if type == 1])
				all_outlines.append( (chromosome,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, transposable_element_genes, transposable_elements) )
			else:
				# If a GFF file was NOT provided, DO NOT include per-window transposon counts in results
				all_outlines.append( (chromosome,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov) )
	
	if mode == 1:   # Poisson fit mode
		outfilename = os.path.join(outdir, str(t_count) + "_read_threshold/" + str(chromosome) + "_Results.txt")
	elif mode == 2: # K-means clustering mode
		outfilename = os.path.join(outdir, str(chromosome) + "_Results.txt")
	all_outlines = sorted(all_outlines, key = lambda entry: (entry[0], entry[1]))
	with open(outfilename, 'w') as outfile:
		if tp_windows != {}:
			# If a GFF file was provided, include per-window transposon counts in results
			line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tTE Genes\tTEs\n"
		else:
			# If a GFF file was NOT provided, DO NOT include per-window transposon counts in results
			line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\n"
		outfile.write(line)
		for line in all_outlines:
			line = "\t".join([str(x) for x in line]) + "\n"
			outfile.write(line)
	print(outfilename, "complete")

def k_means_clustering(hist_dict):
	# k_means_clustering, when given a reads-per-window histogram, attempts to find the best possible
	# clustering of this data to some number of groups 'k'. This is accomplished by using a genetic
	# algorithm. First, a randomized group of k clustering points is created total_individuals times, representing the
	# 'initial population'. These groups of clusters are scored based on the combined squared sum of the
	# data point-cluster point distances. Groups of clustering points are considered to be the 'most fit'
	# based on how low this combined sum is. From the initial population of total_individuals cluster groups, the top few are
	# chosen to pass on to the next generation. A population of total_individuals groups is created for the 2nd generation,
	# with the points in each group being randomly chosen cluster points from the top cluster groups in the 1st generation.
	# Some 'mutation' may occur: with low frequency, the positioning of a cluster can change up to +/- 0.1% of the
	# maximum number of reads-per-window possible. With some lower frequency, the positioning of a cluster can
	# change up to a larger value of +/- 1%. This process of picking the most fit individuals to continue on to a new
	# generation proceeds until each cluster point is virtually the same among the top few individuals, or until
	# 30 generations have passed, at which point the program will try again with a new starting population.
	# If each new starting population fails to converge after a total of three tries, the program aborts, informs
	# the user of the inability to reach an optimum grouping of data points, and prints the best available clusters.
	
	# The correct value of k to be used is determined by starting with k=1, finding the optimum grouping with
	# just one cluster point, and then incrementing k until the most optimum k value is found. This is determined
	# by (a) looking to see if the per-cluster point sum of data point-cluster point distances is no longer
	# decreasing when compared to the previous generation, and by (b) comparing the locations of clustering
	# points after the genetic algorithm has reached a consensus, and noting that if two clustering points
	# are excessively close to one another (<5 apart), we have likely exceeded the optimum k value.
	
	global verbose
	
	# Tweakable variables
	small_mut_freq = 0.18
	large_mut_freq = 0.02
	breeding_individuals = 20
	total_individuals = 300
	total_generations = 500
	minimum_generations = 3
	
	# Will store the optimum number of k cluster points to be used for each win/chromosome combo
	k_cluster_count_dict = {win: Counter() for win in hist_dict} 
	for win in k_cluster_count_dict:
		k_cluster_count_dict[win] = {chromosome: 0 for chromosome in hist_dict[win]}
	
	# Will store a list of genetic algorithm-optimized cluster groups for each win/chromosome combo
	# as well as the dist_sum value and the number of data points that belong to each cluster
	k_cluster_group_dict = {win: {} for win in hist_dict}
	for win in k_cluster_group_dict:
		k_cluster_group_dict[win] = {chromosome: [] for chromosome in hist_dict[win]}
	
	# Find optimum data clustering for each window & chromosome combination
	for win in hist_dict:
		for chromosome in hist_dict[win]:
			# Remove histogram elements in which a coverage value appears with a frequency of 0, then
			# calculate the minimum and maximum reads-per-window values that can be found on this chromosome
			hist_dict[win][chromosome] = Counter({cov: freq for cov, freq in hist_dict[win][chromosome].items() if freq != 0})
			min_reads = min(hist_dict[win][chromosome].keys())
			max_reads = max(hist_dict[win][chromosome].keys())
			max_reads_95th = sorted(hist_dict[win][chromosome].keys())[int(len(hist_dict[win][chromosome])*.95)] # 95th percentile reads-per-window value
			
			# Find the optimum number of k clustering points. The following iteratively runs a genetic algorithm
			# to optimize the grouping of data points to k clusters for as many k values as is necessary. This
			# is done until the breeding population for the next generation is the same as the previous generation,
			# or until the distance between close cluster points is insignificant. If a value of k yields
			# insignificant cluster point distances (defined as < 1% max # of possible reads mapping per window)
			# between two cluster points, k-1 is the ideal number of clusters. Alternatively, if the average data
			# point-cluster point sum, multiplied by k, is higher than the sum*k given k-1 clusters, k-1 would
			# again be set as the ideal number of clusters. This multiplication by k of the sum is simply a penalty
			# score that prevents k from increasing far too high.
			
			# best_cluster_group_for_each_k will hold optimum cluster group and the number of assigned data points for each k value
			# until optimum k is determined. Optimum k's clusters will be stored in k_cluster_group_dict
			best_cluster_group_for_each_k = {}
			
			k = 0
			while k < 50:
				k += 1
				break_for_next_k = False
				initial_pop_attempt = 1
				past_gen_top_breeders = {} # Will be used to determine when genetic algorithm has converged on the same answer for 3 generations. Key = gen #, Value = Set of top breeders
				while True:
					# Generate an initial population of total_individuals groups of k # of clusters, with each cluster assigned a
					# random integer position between the minimum and maximum possible values for reads mapping per window
					cluster_group_list = []
					for i in range(total_individuals):
						cluster_group = []
						for j in range(k):
							cluster_group.append(random.randrange(min_reads,max_reads))
						cluster_group = sorted(cluster_group)
						dist_sum = 0 # Will hold total distance of each cluster point from their data points
						data_count_dict = Counter()
						if k == 1:
							# Code optimized for when only one cluster point is used
							for cluster_point in cluster_group:
								dist_sum += sum( [(int(fabs(data_point - cluster_point))+2)**2 * freq for data_point, freq in hist_dict[win][chromosome].items()] )
								data_count_dict[cluster_point] = sum( [freq for freq in hist_dict[win][chromosome].values()] )
						else:
							# Code for when >1 clusters are present. Data points must be placed into the closest cluster
							dist_sum, data_count_dict = cluster_group_assign_data_points(cluster_group, hist_dict[win][chromosome], max_reads)
						cluster_group_list.append([cluster_group, dist_sum, data_count_dict])
					cluster_group_list = sorted(cluster_group_list, key = lambda cluster_group_list: cluster_group_list[1])
					
					# Iteratively create new generations given the most successful individuals in the previous generation.
					# Do not let this iterative process continue beyond the total allowable generations. If this limit is
					# reached, a new initial population will be created and the iterative process will begin again.
					gen_number = 1
					past_gen_top_breeders[gen_number] = cluster_group_list[:breeding_individuals]
					while gen_number < total_generations:
						gen_number += 1
						cluster_group_list = cluster_next_generation(cluster_group_list, hist_dict[win][chromosome], breeding_individuals, total_individuals, small_mut_freq, large_mut_freq, k, max_reads)
						past_gen_top_breeders[gen_number] = cluster_group_list[:breeding_individuals]
						if gen_number > 3:
							del past_gen_top_breeders[gen_number-3] # Save memory and time
						
						# Check several criteria to see if the genetic algorithm has reached a consensus on cluster location
						# If so, save the single highest scoring cluster for this k value, then check to see if previous k values
						# scored better. If so, set k-1 as the optimum clustering value. If not, continue increasing k.
						
						# 1) For the beginning few generations, check to see if the top breeding individuals are all the same and
						#	 are identical to the previous two generations
						if gen_number < total_generations / 3 and gen_number > 2:
							if gen_number >= minimum_generations and past_gen_top_breeders[gen_number] == past_gen_top_breeders[gen_number-1] and past_gen_top_breeders[gen_number] == past_gen_top_breeders[gen_number-2]:
								#and all(x == past_gen_top_breeders[gen_number][0] for x in past_gen_top_breeders[gen_number])
								best_cluster_group = cluster_group_list[0]
								best_cluster_group_for_each_k[k] = best_cluster_group
								break_for_next_k = True
								if verbose == True:
									print("Consensus reached for k=", str(k), " in ", gen_number, " generations for ", str(win), " bp windows, ", str(chromosome), " due to identity among breeding individuals to the previous 2 generations.\n",
											str(best_cluster_group[0]), "\t", str(best_cluster_group[1]), "\t", str(sorted(best_cluster_group[2].items())), "\t", str(log(best_cluster_group_for_each_k[k][1] * k)), "\n",sep='')
								break
						
						# 2) If the program has been running for a while without finding a perfect consensus among the top breeding individuals,
						#	 start checking to see if the average distance between cluster points groups among the breeding individuals is <2%
						#	 of the 95th percentile number of reads possible. If so, this indicates very little change is occurring among the
						#	 top individuals, so we may as well consider this list converged.
						elif gen_number > total_generations / 3:
							avg_dist = 0
							for i in range(1,breeding_individuals):
								for j in range(k):
									cur_cluster = int(cluster_group_list[i][0][j])
									first_cluster = int(cluster_group_list[0][0][j])
									avg_dist += int(fabs(cur_cluster - first_cluster))
							avg_dist /= ( (breeding_individuals - 1) * k )
							if avg_dist < 0.02 * max_reads_95th and gen_number >= minimum_generations:
								best_cluster_group = cluster_group_list[0]
								best_cluster_group_for_each_k[k] = best_cluster_group
								break_for_next_k = True
								if verbose == True:
									print("Consensus reached for k=", str(k), " in ", gen_number, " generations for ", str(win), " bp windows, ", str(chromosome), " due to negligible cluster point changes.\n",
											str(best_cluster_group[0]), "\t", str(best_cluster_group[1]), "\t", str(sorted(best_cluster_group[2].items())), "\t",
											str(log(best_cluster_group_for_each_k[k][1] * k)), "\n",sep='')
								break
					
					if break_for_next_k == True:
						# If break_for_next_k == True:
						# gen_number < total_generations. The clusters converged.
						break
					
					# If break_for_next_k == False:
					# gen_number == total_generations. The clusters did not converge, so create a new initial population and try again.
					initial_pop_attempt += 1
					print("Failed to converge within ", str(total_generations), " generations. Starting Try #", str(initial_pop_attempt), " with a new initial population.")
				
				if k > 1 and log(best_cluster_group_for_each_k[k-1][1] * (k-1)) < log(best_cluster_group_for_each_k[k][1] * k):
					# Determine whether or not k-1 dist_sum (normalized for the number of clusters used) was the optimal value for k.
					# If so, halt genetic algorithm and return the k-1 cluster groups to be used for the creation of HMM probabilities.
					k_cluster_count_dict[win][chromosome] = k-1
					k_cluster_group_dict[win][chromosome] = best_cluster_group_for_each_k[k-1]
					print("k=", str(k-1), " is optimal for ", str(win), " bp windows, ", str(chromosome), ":\n", str(best_cluster_group_for_each_k[k-1][0]), "\t", str(best_cluster_group_for_each_k[k-1][1]), "\t", sorted(best_cluster_group_for_each_k[k-1][2].items()), "\n", sep='')
					break
			if k == 50:
				print("Failed to converge with 50 clusters or less. Program will now end.")
				sys.exit(0)
	
	return(k_cluster_count_dict, k_cluster_group_dict)

def prepare_hmm(mode, sam, tp_positions, chr_len, chromosome, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir):
	# prepare_hmm creates (1) per-window read count histograms and (2) observation files, as well as
	# (3) details for the HMM parameter file such as read count thresholds and copy number state names.
	#
	# At the end of this function, dist_to_params is called to output this data to files, and finally
	# hmm_dup_search is called to run the HMM on these files and output the duplication search results.
	
	# Create list of chromosomes to run HMM on
	chromosome_list = []
	if chromosome == "All":
		chromosome_list = [chromosome_entry for chromosome_entry in chr_len]
	else:
		chromosome_list = [chromosome_entry for chromosome_entry in chromosome.split(',')]
	
	# Create the HMM observation files (path_dict, a dictionary containings # of reads per window in sequential order)
	# and reads-per-window histograms (hist_dict, a dictionary whose values are chromosome-specific Counters representing reads-per-window)
	hist_dict, path_dict = sam_to_wincov(sam, win_list, chromosome_list, outdir)
	
	# Read in histograms if hist_dict is empty for a particular window. This happens when histograms
	# had already been generated for all chromosomes for that window, and is simply a time-saving feature
	# when performing re-runs on the same data set. Histogram and observation files are not re-created if
	# they already exist, as this is the most time-consuming part of the entire pipeline.
	for win in win_list:
		if not hist_dict[win]:
			# If this dictionary entry is empty for some window, grab those histograms from files
			for chromosome in chromosome_list:
				hist_i_fn = os.path.join(outdir, str(win) + "bp_window/" + str(chromosome) + "_hist.txt")
				with open(hist_i_fn) as infile:
					inlines = infile.readlines()
					hist = Counter()
					for line in inlines[1:]:
						cov, freq = line.split("\t")
						hist[int(cov)] = int(freq)
				hist_dict[win][chromosome] = hist

	if mode == 1: # Poisson fit mode	
		# Determine the poisson lambda values for a distribution of reads given the aforementioned reads-per-window histogram tables.
		# If the user supplied a 1X copy region lambda value, use that instead.
		singlecopy_lambda = {win: {} for win in win_list}
		for win in win_list:
			for chromosome in chromosome_list:
				# Dictionary holding 1X copy lambdas, used for threshold calculations later
				singlecopy_lambda[win][chromosome] = pois_lambda
		
		# If user didn't supply 1X lambda, calculate it from chromosomal per-window read histograms
		if pois_lambda == 0:
			for win in win_list:
				for chromosome in chromosome_list:
					# First create text file for R script containing the filename for the read histogram
					o_fn = "HMMInputfiles.txt"
					i_fn_1 = os.path.join(outdir,str(win) + "bp_window/" + str(chromosome) + "_hist_distfits.txt")
					if not os.path.exists(i_fn_1):
						i_fn_2 = os.path.join(outdir, str(win) + "bp_window/" + str(chromosome) + "_hist.txt")
						line = str(i_fn_2) + "\n"
						with open(o_fn, 'w') as outfile:
							outfile.write(line)
						# Now call R script "dist_fit.R" to perform Poisson regression
						Rpath = "Rscript"
						params = ' '.join([str(Rpath) + " dist_fit.R", "--no-save"])
						simulation = subprocess.Popen(params, shell=True)
						simulation.wait()
						# Try to remove temporary file "HMMInputfiles.txt"
						try:
							os.remove(o_fn)
						except:
							print("Due to an issue with permissions, the temporary file 'HMMInputfiles.txt' could not be removed.")
					# Read in "dist_fit.R" results file and grab Poisson regression lambda value
					with open(i_fn_1) as infile:
						inlines = infile.readlines()
						inlines = [line.split() for line in inlines]
						singlecopy_lambda[win][chromosome] = float(inlines[1][1])
	
		# Assign read count thresholds for window/chromosome combinations
		temp_t_count = {win: {} for win in win_list}
		for win in win_list:
			for chromosome in chromosome_list:
				temp_t_count[win][chromosome] = t_count[:]
				if t_count == [0]:
					# If no threshold was given, assign this to maximum coverage
					temp_t_count[win][chromosome] = [max(hist_dict[win][chromosome])]
				else:
					# Calculate read count threshold based on X-fold value of single-copy region lambda value
					for i in range(0,len(t_count)):
						if t_count[i][-1:] == "X":
							temp_t_count[win][chromosome][i] = int(float(int(t_count[i][:-1]) * singlecopy_lambda[win][chromosome]))
	
		# Assign states and their corresponding character codes
		state_cns_list = state_cns.split(',')
		state_names_list = state_names.split(',')
		states = OrderedDict()
		for i in range(0,len(state_cns_list)):
			states[float(state_cns_list[i])] = state_names_list[i]
		
		# Cycle over each window, chromosome, and threshold value combination provided and run HMM for each combination
		for win in win_list:
			for chromosome in chromosome_list:
				for threshold_count in temp_t_count[win][chromosome]:
					# Create an HMM parameter file (poisson probabilities of given # of reads in a window based on duplication status of that window)
					threshold_count = dist_to_params_mode_one(states, singlecopy_lambda[win][chromosome], trans_prob, threshold_count, hist_dict[win][chromosome], win, chromosome, sam, outdir)
					
					# Run HMM using generated parameter file and observation file
					hmm_dup_search(mode, outdir, tp_positions, chromosome, win, singlecopy_lambda[win][chromosome], states, threshold_count)
	
	elif mode == 2: # K-means clustering mode
		# Pass histograms to k-means clustering function w/ genetic algorithm-based positioning of groups
		k_cluster_count_dict, k_cluster_group_dict = k_means_clustering(hist_dict)
		
		# Cycle over each window, chromosome, and threshold value combination provided and run HMM for each combination
		for win in win_list:
			for chromosome in chromosome_list:
				# Create an HMM parameter file (poisson probabilities of given # of reads in a window based on duplication status of that window)
				dist_to_params_mode_two(hist_dict[win][chromosome], k_cluster_count_dict[win][chromosome], k_cluster_group_dict[win][chromosome][0], k_cluster_group_dict[win][chromosome][2], win, chromosome, sam, outdir)
				
				# Run HMM using generated parameter file and observation file
				k_cluster_group_dict_temp = k_cluster_group_dict[win][chromosome][0]
				k_data_count_dict = k_cluster_group_dict[win][chromosome][2]
				# Assign cluster group with the most data points as the '1X' region, around which all other cluster groups
				# will have their relative fold coverage calculated. May not actually represent 1X region, but this at
				# least gives all other CNVs a baseline to be compared to.
				single_copy_lambda = sorted(k_data_count_dict.items(), key = lambda entry: entry[1], reverse=True)[0][0]
				hmm_dup_search(mode, outdir, tp_positions, chromosome, win, single_copy_lambda, k_cluster_group_dict_temp, 0)

def sam_to_wincov(sam_file, win_list, chromosome_list, outdirpart):
	# Check to see if, for a given window size, all chromosomes have already had their
	# per-window read count histogram generated. If so, do not re-create the existing
	# data files necessary for the HMM to run.
	win_list_needed = []
	for win in win_list:
		not_needed = 0
		for chromosome in chromosome_list:
			# Check here to see if an observation file exists. If so, assume all files have
			# already been generated for this chromosome & window combination
			i_fn = os.path.join(outdirpart, str(win) + "bp_window/" + str(chromosome) + "_obs.txt")
			if os.path.exists(i_fn): not_needed += 1
		# If we did not see an observation file for all chromosomes to be run on, re-create all
		# the HMM observation and histogram files
		if not_needed != len(chromosome_list): win_list_needed.append(win)
	
	hist_dict = {win: {} for win in win_list}
	path_dict = {win: {} for win in win_list}
	chr_string = ','.join(chromosome_list)
	# If list is not empty, create paths for each chromosome/window pair
	# If list is empty, we don't need to re-create existing data files, and thus the following 'if' block is skipped
	if win_list_needed:
		# Get per-window read count histograms and paths for the given window sizes
		hist_dict, path_dict = SAM.coverage_window(sam_file, win_list_needed, chr_string)
		
		# Check to see if folders have already been made for the provided windows. If not, make them.
		for chromosome in chromosome_list:
			for win in win_list:
				outdir = os.path.join(outdirpart,str(win) + "bp_window/")
				print(outdirpart)
				print(outdir)
				if not os.path.exists(outdir): os.makedirs(outdir)
		
		# Write out per-window read count histograms and paths
		for win in win_list:
			outdir = os.path.join(outdirpart,str(win) + "bp_window/")
			for chromosome in chromosome_list:
				outfilename1 = os.path.join(outdir, str(chromosome) + "_obs.txt")
				outfilename2 = os.path.join(outdir, str(chromosome) + "_hist.txt")
				with open(outfilename1, 'w') as outfile:
					outfile.write(path_dict[win][chromosome])
				
				outline = "\n".join([str(cov) + "\t" + str(freq) for cov, freq in sorted(hist_dict[win][chromosome].items(), key = lambda cov: cov[0])])
				with open(outfilename2, 'w') as outfile:
					outline2 = "Coverage_per_" + str(win) + "bp_window\tCounts\n"
					outfile.write(outline2)
					outfile.write(outline)
	
	return(hist_dict, path_dict)

# Grab command line options
mode, sam, gff, chromosome, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir = command_line()

# Determine the names of chromosomes and their lengths
sam_parser = SAM.Parse(sam)
sam_parser.chr_lengths()
sam_parser.start()
chr_len = sam_parser.get_chr_lengths()

# If a GFF file is provided, determine the location of all transposons
# This is useful for quickly identifying false positive duplication calls
tp_positions = grab_transposons(gff, chr_len)

# Generate the necessary data to perform an HMM-based duplication search,
# then run the HMM on these parameters and observations.
# Detailed commentary is provided within sub-routines
prepare_hmm(mode, sam, tp_positions, chr_len, chromosome, win_list, state_cns, state_names, pois_lambda, trans_prob, t_count, outdir)