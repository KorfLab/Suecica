#!/usr/bin/env python3.2
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import argparse, re, random, os, subprocess, sys
from collections import Counter, OrderedDict
from math import exp, factorial, log, lgamma, fabs
import SAM
verbose = False

# DupHMM.py runs a duplication-finding HMM on a SAM file. The SAM file is first broken up into user-specified window sizes,
#	  reads mapping per-window counted, and a Poisson regression is performed on the resulting histogram of count data. This Poisson lambda
#	  value is then used to generate a HMM parameter file, which the duplication searching HMM then uses to identify duplicated windows.

def Cluster_Group_Assign_Data_Points(cluster_group, hist_dict, max_reads):
	# This function groups data points to the appropriate cluster point and returns the data point-cluster point distance sum
	
	data_count_dict = Counter({cluster: 0 for cluster in cluster_group})
	optimum_cluster_dists = {cluster: [] for cluster in cluster_group} # Value = List containing (data point-cluster distance**2 * freq)
	for data_point, freq in hist_dict.items():
		# Cycle through each data point, finding which cluster point each is closest to
		# and storing (distance * freq) to that cluster's list of data points
		min_cluster = -1												# Will hold closest cluster to data point. -1 is a default 'not assigned yet' value
		min_cluster_dist = (fabs(data_point - cluster_group[0])+1)**2+1	# Will hold closest cluster's distance from data point. Default is larger than 1st entry
		#min_cluster_dist = (data_point - cluster_group[0])**2+1	# Will hold closest cluster's distance from data point. Default is larger than 1st entry
		for cluster_point in cluster_group:
			data_cluster_dist = (fabs(data_point - cluster_point)+1)**2
			#data_cluster_dist = (data_point - cluster_point)**2
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

def Cluster_Mutate(cluster_group, small_mut_freq, large_mut_freq, max_reads):
	# With low frequency , shift cluster points a bit
	# Maximum small mutation shift distance is 1% of maximum number of reads per window possible
	# Maximum large mutation shift distance is 2.5% of maximum number of reads per window possible
	
	new_cluster_group = cluster_group
	for i in range(len(cluster_group)):
		random_value = random.random()
		shift_dist = 0
		if random_value < large_mut_freq:
			shift_dist = int(random.randrange(int(2 * max_reads*0.025)) - max_reads*0.025)
		elif random_value < small_mut_freq:
			shift_dist = int(random.randrange(int(2 * max_reads*0.01)) - max_reads*0.01)
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

def Cluster_Next_Generation(cluster_group_list, hist_dict, breeding_individuals, total_individuals, small_mut_freq, large_mut_freq, k, max_reads):
	# Create a new generation of cluster points given a breeding population of clusters
	
	# From the previous generation's total_individuals groups of k # of clusters, pick the top cluster groups which
	# have the lowest data point-cluster point distance sum to populate the next generation
	next_gen_cluster_group_list = sorted(cluster_group_list, key = lambda cluster_group_list: cluster_group_list[1])[:breeding_individuals]
	
	# With low frequency, shift breeding population cluster points a bit
	for i in range(len(next_gen_cluster_group_list)):
		cluster_group, dist_sum, data_count_dict = next_gen_cluster_group_list[i]
		new_cluster_group = Cluster_Mutate(cluster_group, small_mut_freq, large_mut_freq, max_reads)
		if cluster_group != new_cluster_group:
			new_cluster_set = set(new_cluster_group) # Sets do not allow for duplicate clusters
			if len(new_cluster_set) == len(new_cluster_group):
				new_dist_sum, new_data_count_dict = Cluster_Group_Assign_Data_Points(new_cluster_group, hist_dict, max_reads)
				next_gen_cluster_group_list[i] = [new_cluster_group, new_dist_sum, new_data_count_dict]
	
	# Create (total_individuals-#_of_breeding_clusters) individuals. Their cluster points will result from 'crossovers' among the
	# breeding individuals. The next generation individuals are also subject to 'mutational' position drifting
	cluster_list = []
	for cluster_group, dist_sum, data_count_dict in next_gen_cluster_group_list:
		for cluster in cluster_group:
			cluster_list.append(cluster)
	
	for i in range(breeding_individuals, total_individuals):
		unique_indiv_created = False
		while unique_indiv_created == False:
			new_cluster_group = random.sample(cluster_list,k)
			new_cluster_mutated_group = Cluster_Mutate(new_cluster_group, small_mut_freq, large_mut_freq, max_reads)
			new_cluster_mutated_set = set(new_cluster_mutated_group) # Sets do not allow for duplicate clusters
			if len(new_cluster_mutated_set) == len(new_cluster_mutated_group):
				new_dist_sum, new_data_count_dict = Cluster_Group_Assign_Data_Points(new_cluster_mutated_group, hist_dict, max_reads)
				next_gen_cluster_group_list.append([new_cluster_mutated_group, new_dist_sum, new_data_count_dict]) # Note: mutated group may or may not be mutated
				unique_indiv_created = True
	
	next_gen_cluster_group_list = sorted(next_gen_cluster_group_list, key = lambda next_gen_cluster_group_list: next_gen_cluster_group_list[1])
	
	return next_gen_cluster_group_list

def Command_line():
	parser = argparse.ArgumentParser(description="DupHMM.py runs a duplication-finding HMM on a SAM file. The SAM file is first broken up into user-specified window sizes, \
									 reads mapping per-window counted, and a Poisson regression is performed on the resulting histogram of count data. This Poisson lambda \
									 value is used to generate a HMM parameter file, which the duplication searching HMM then uses.")
	parser.add_argument('-sam', default="Data/Sue/sue1/single_bp_preprocess/set5/libSUE1_set5_aln_Aa_Filtered.sam.gz", type=str, help='SAM filename for DupHMM to run on.', metavar='SAMFilename')
	parser.add_argument('-gff', default="../Thalyrata.gff", type=str, help='GFF file for the reference genome. Used to annotate regions in the results which contain transposons. This makes spotting false positives easier.', metavar='GFFFile')
	parser.add_argument('-chr', default="AtChr5", type=str, help='Chromosome to run HMM on (Default=Run HMM on all chromosomes). To look at >1 chromosome but not all chromosomes, provide them in a comma-delimited format: e.g. "Chr1,Chr2"', metavar='Chromosome')
	parser.add_argument('-w', default="500", type=str, help='Window size the HMM will use, in bps (Default=500). To look at multiple windows, provide them in a comma-delimited format: e.g. "500,2000,2000,2500".', metavar='Window')
	parser.add_argument('-v', action='store_true', help='Turns on verbose output for k-means clustering')
	parser.add_argument('-o', default="Output/DupHMM_sue1_set5_AaFiltered_Kmeans_clustering", type=str, help='Output directory for DupHMM', metavar='OutputDir')
	parser.add_argument('-R', default="", type=str, help="Directory where R is located. This only needs to be entered once, as the path is then saved in 'Rpath.txt'. NOTE: Your R installation must have the 'MASS' package installed.", metavar='RPath')
	
	args = parser.parse_args()
	sam = args.sam
	gff = args.gff
	chr = args.chr
	win_list = args.w
	win_list = [int(entry) for entry in win_list.split(',')]
	verbose = args.v # PROPOGATE THIS THROUGHOUT THE PROGRAM!
	outdir = args.o
	if outdir[-1:] == "\"": outdir = outdir[:-1]
	Rpath = args.R
	if Rpath[-1:] == "\"": Rpath = Rpath[:-1]
	#Rpath = r"/share/apps/R-2.15.0/bin/Rscript"
	#Rpath = r"C:/Program Files/R/R-2.15.2/bin/x64/Rscript.exe"
	
	# Determine the path for R installation
	if Rpath == "":
		try:
			with open("Rpath.txt") as infile:
				Rpath = infile.readline()
		except:
			print("A path to R was not provided, and 'Rpath.txt' does not exist.")
			sys.exit(0)
	else:
		with open("Rpath.txt", 'w') as outfile:
			outfile.write(Rpath)
	
	return(sam, gff, chr, win_list, outdir, Rpath)

def Dist_To_Params(hist, k_cluster_count, k_loc, k_data_count, win, chr, sam, outdir):
	# Distribution_to_parameters_file
	# Converts per-window read count histogram with a Poisson fit to it into
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
					max_value = t_count = k-1
					break
			prob_list.append(prob)
		prob_dict[cluster_pos] = prob_list
		j += 1
	
	# Trim histograms such that the copy number states do not have an excessive number of right-tailed 0's
	# HMM treats missing reads-per-window counts as 0, so right-sided trimming histograms to the maximum read-per-window read count
	# with a non-zero probability for the highest copy number state helps save parameter file disk space, time, and memory
	i = 0
	for cluster_pos in k_loc:
		prob_dict[cluster_pos] = {i: prob_dict[cluster_pos][i] for i in range(0,t_count+1)}
		# Do not need to re-trim highest copy-number state
		i += 1
		if i == len(k_loc) - 1: break
	
	# Create transition probabilities between states based on relative percentages of each state
	total_data_points = sum([point_sum for point_sum in k_data_count.values()])
	trans_prob_dict = {cluster: {} for cluster in k_loc}
	for cluster in k_loc:
		trans_prob_dict[cluster] = {x: 10**(k_data_count[x] / total_data_points * 75) for x in k_loc}
	for cluster in k_loc:
		trans_prob_dict[cluster][cluster] = 1 - (len(k_loc) - 1) * 10**(k_data_count[x] / total_data_points * 75)
	
	# Create parameter file
	# Call Dist_To_Params_Outline to create the line of text to be written out
	outline = Dist_To_Params_Outline(k_loc, prob_dict, trans_prob_dict)
	t_outdir = os.path.join(outdir,str(win) + "bp-window/" + str(t_count) + "_read_threshold/")
	if not os.path.exists(t_outdir): os.makedirs(t_outdir)
	outfilename = os.path.join(t_outdir, str(chr) + "_params.txt")
	with open(outfilename, 'w') as outfile:
		outfile.write(outline)
	return(max_value)

def Dist_To_Params_Outline(k_loc, prob_dict, trans_prob_dict):
	# Distribution_to_parameters_file_outline
	# Formats HMM parameters for output into a format that the HMM can read in
	
	# States and start probabilities
	# ASCII 65 = 'A'
	outline = "States:\n"
	states = [chr(x) for x in range(65,len(k_loc))]
	outline += '\n'.join(states)
	outline += "\n\nStart Probabilities:\n"
	for i in states:
		outline += str(i) + "\t" + str(1/len(states)) + "\n"
	
	# Transition probabilities
	outline += "\nTransition Probabilities:"
	for i in range(0,len(states)):
		outline += "\n" + str(states[i]) + "\t" + str(states[0]) + "\t" + str(trans_prob_dict[k_loc[i]] [k_loc[0]] )
		for j in range(1,len(states)):
			outline += "\n\t" + str(states[j]) + "\t" + str(trans_prob_dict[k_loc[i]] [k_loc[j]] )
	
	# Part of output containing emission probabilities
	outline += "\n\nEmission Probabilities:\n"
	for cluster in k_loc:
		prob_list = prob_dict[cluster]
		outline += str(states[cluster])
		outline2 = ''.join(["\t" + str(i) + "\t" + str(prob_list[i]) + "\n" for i in range(0,len(prob_dict[cluster]))])
		outline += str(outline2) + "\n"
	return outline

def Grab_Transposons(gff, chr_len):
	pattern = r"(transposable_element_gene|transposable_element)"
	recomp = re.compile(pattern)
	pattern2 = r"^ID=(\S+?);"
	recomp2 = re.compile(pattern2)
	tp_positions = {}
	if gff != "":
		# If a GFF file was provided, fill tp_positions dictionary
		# with transposon positions
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
					elif match.group(1) == "transposable_element": type = 0
					tp_positions[chr].append((name, spos,epos, type))
	return tp_positions

def HMM_Dup_Search(folder, tp_positions, chr, win, t_count, k_loc_dict, k_data_count_dict):
	# HMM_Dup_search runs the HMM using a Viterbi algorithm written in C++ .
	# It waits for the HMM run to complete, grabs the HMM results from Viterbi's output, then parses and outputs results.
	
	outdir = os.path.join(folder, str(win) + "bp-window/")
	param_file = os.path.join(outdir, str(t_count) + "_read_threshold/" + str(chr) + "_params.txt")
	obs_file = os.path.join(outdir, str(chr) + "_obs.txt")
	
	# Run Viterbi algorithm written in C++
	if sys.platform in ["linux","linux2"] or sys.platform in ["darwin", "os2", "os2emx"]:
		# Linux and Mac
		run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),"Viterbi")
		params = ' '.join([str(run_script), "-p", str(param_file), "-o", str(obs_file)])
		HMMRun = subprocess.Popen(params, shell=True, stdout=subprocess.PIPE)
	elif sys.platform in ["win32","win64"]:
		# Windows
		run_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),"Viterbi.exe")
		params = ' '.join([str(run_script), "-p", str(param_file), "-o", str(obs_file)])
		HMMRun = subprocess.Popen(params, shell=False, stdout=subprocess.PIPE)
	else:
		# User is operating on an unknown OS, end program
		print("Working with an unknown OS, ", sys.platform, ", unable to run Viterbi algorithm. Please report this problem.", sep='')
		sys.exit(0)
	path = HMMRun.communicate()[0].decode("utf-8") # Grab path from Viterbi.exe output
	path = path[::-1]
	
	# Open up observation file and grab observations for average fold coverage calculations across multiple windows
	obs_list = []
	with open(obs_file) as infile:
		inline = infile.readlines()
		inline = inline[0].split(" ")
		obs_list = [int(inline[i]) for i in range(0,len(inline))]
	
	# If a GFF file was provided, figure out which windows contain transposons
	tp_windows = {}
	if tp_positions != {}:
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
	# ASCII 65 = 'A'
	states = [chr(x) for x in range(65,len(k_loc))]
	single_copy_lambda = sorted(k_data_point_dict.items(), key = lambda entry: entry[1], reverse=True)[0][0]
	for state in states:
		pattern = str(states[state]) + r"{1,}"
		recomp = re.compile(pattern)
		
		for m in recomp.finditer(path):
			s_pos = m.start()*win
			e_pos = m.end()*win
			fold_cov = round(sum([obs_list[i] for i in range(m.start(),m.end()+1) if obs_list[i] < t_count]) / (m.end() - m.start()) / single_copy_lambda, 2)
			tp_union_list = set()
			# If a GFF file was provided, include per-window transposon counts in results
			if tp_windows != {}:
				for i in range(m.start(),m.end()+1):
					for (transposon_name, type) in tp_windows[i]:
						# Only counting the number of unique transposon types present in a window
						tp_union_list.add((transposon_name, type))
				transposable_elements = len([name for name, type in tp_union_list if type == 0])
				transposable_element_genes = len([name for name, type in tp_union_list if type == 1])
				if len(path_dict_exp_control_ratio.keys()) == 0:
					# Not using a control
					all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, transposable_element_genes, transposable_elements) )
				else:
					# Using a control
					exp_control_ratio_avg = round(sum([float(path_dict_exp_control_ratio[win][chr][i]) for i in range(m.start(),m.end()+1)]) / (m.end() - m.start()), 2)
					all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, transposable_element_genes, transposable_elements, exp_control_ratio_avg) )
			else:
				# If a GFF file was NOT provided, DO NOT include per-window transposon counts in results
				if len(path_dict_exp_control_ratio.keys()) == 0:
					# Not using a control
					all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov) )
				else:
					# Using a control
					exp_control_ratio_avg = round(sum([float(path_dict_exp_control_ratio[win][chr][i]) for i in range(m.start(),m.end()+1)]) / (m.end() - m.start()), 2)
					all_outlines.append( (chr,s_pos,e_pos,e_pos - s_pos, str(state), fold_cov, exp_control_ratio_avg) )
	
	outfilename = os.path.join(outdir, str(t_count) + "_read_threshold/" + str(chr) + "_Results.txt")
	all_outlines = sorted(all_outlines, key = lambda entry: (entry[0], entry[1]))
	with open(outfilename, 'w') as outfile:
		if tp_windows != {}:
			# If a GFF file was provided, include per-window transposon counts in results
			if len(path_dict_exp_control_ratio.keys()) == 0:
				# Not using a control
				line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tTE Genes\tTEs\n"
			else:
				# Using a control
				line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tTE Genes\tTEs\tExp:Control\n"
		else:
			# If a GFF file was NOT provided, DO NOT include per-window transposon counts in results
			if len(path_dict_exp_control_ratio.keys()) == 0:
				# Not using a control
				line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\n"
			else:
				# Using a control
				line = "Chr\tStart_pos\tEnd_pos\tDistance\tState\tFold_Coverage\tExp:Control\n"
		outfile.write(line)
		for line in all_outlines:
			line = "\t".join([str(x) for x in line]) + "\n"
			outfile.write(line)
	print(outfilename, "complete")

def K_Means_Clustering(hist_dict):
	# K_Means_Clustering, when given a reads-per-window histogram, attempts to find the best possible
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
	
	# Tweakable variables
	small_mut_freq = 0.02
	large_mut_freq = 0.01
	breeding_individuals = 20
	total_individuals = 300
	total_generations = 300
	
	# Will store the optimum number of k cluster points to be used for each win/chr combo
	k_cluster_count_dict = {win: Counter() for win in hist_dict} 
	for win in k_cluster_count_dict:
		k_cluster_count_dict[win] = {chr: 0 for chr in hist_dict[win]}
	
	# Will store a list of genetic algorithm-optimized cluster locations for each win/chr combo
	# as well as the number of data points that belong to each cluster
	k_loc_dict = {win: {} for win in hist_dict}
	for win in k_loc_dict:
		k_loc_dict[win] = {chr: [] for chr in hist_dict[win]}
	
	# Find optimum data clustering for each window & chromosome combination
	for win in hist_dict:
		for chr in hist_dict[win]:
			# Remove histogram elements in which a coverage value appears with a frequency of 0, then
			# calculate the minimum and maximum reads-per-window values that can be found on this chromosome
			hist_dict[win][chr] = Counter({cov: freq for cov, freq in hist_dict[win][chr].items() if freq != 0})
			min_reads = min(hist_dict[win][chr].keys())
			max_reads = max(hist_dict[win][chr].keys())
			max_reads_95th = sorted(hist_dict[win][chr].keys())[int(len(hist_dict[win][chr])*.95)] # 95th percentile reads-per-window value
			
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
			# until optimum k is determined. Optimum k's clusters will be stored in k_loc_dict
			best_cluster_group_for_each_k = {}
			
			k = 0
			while k < 50:
				k += 1
				break_for_next_k = False
				initial_pop_attempt = 1
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
								dist_sum += sum( [(data_point - cluster_point)**2 * freq for data_point, freq in hist_dict[win][chr].items()] )
								data_count_dict[cluster_point] = sum( [freq for freq in hist_dict[win][chr].values()] )
						else:
							# Code for when >1 clusters are present. Data points must be placed into the closest cluster
							dist_sum, data_count_dict = Cluster_Group_Assign_Data_Points(cluster_group, hist_dict[win][chr], max_reads)
						cluster_group_list.append([cluster_group, dist_sum, data_count_dict])
					
					# Iteratively create new generations given the most successful individuals in the previous generation.
					# Do not let this iterative process continue beyond the total allowable generations. If this limit is
					# reached, a new initial population will be created and the iterative process will begin again.
					gen_number = 1
					prev_gen_cluster_group_list = cluster_group_list
					while gen_number < total_generations:
						gen_number += 1
						cluster_group_list = Cluster_Next_Generation(cluster_group_list, hist_dict[win][chr], breeding_individuals, total_individuals, small_mut_freq, large_mut_freq, k, max_reads)
						
						# Check several criteria to see if the genetic algorithm has reached a consensus on cluster location
						# If so, save the single highest scoring cluster for this k value, then check to see if previous k values
						# scored better. If so, set k-1 as the optimum clustering value. If not, continue increasing k.
						
						# 1) For the beginning few generations, check to see if the top breeding individuals have cluster groups that are
						#	 identical to each other and to the previous generation
						if gen_number < total_generations / 3:
							prev_gen_tmp = [cluster_group for cluster_group, dist_sum, data_point_dict in sorted(prev_gen_cluster_group_list, key = lambda prev_gen_cluster_group_list: prev_gen_cluster_group_list[1])[:breeding_individuals]]
							total_matches = 0
							for cluster_group, dist_sum, data_point_dict in cluster_group_list[:breeding_individuals]:
								if cluster_group in prev_gen_tmp and cluster_group == cluster_group_list[0][0]:
									total_matches += 1
							if total_matches == breeding_individuals:
								best_cluster_group = cluster_group_list[0]
								best_cluster_group_for_each_k[k] = best_cluster_group
								break_for_next_k = True
								print("Consensus reached for k=", str(k), " in ", gen_number-1, " generations due to identity among breeding individuals to themselves and to the previous generation.\n",
									  str(best_cluster_group), "\n",sep='') # Delete me
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
									avg_dist += fabs(cur_cluster - first_cluster)
							avg_dist /= ( (breeding_individuals - 1) * k )
							if avg_dist < 0.02 * max_reads_95th:
								best_cluster_group = cluster_group_list[0]
								best_cluster_group_for_each_k[k] = best_cluster_group
								break_for_next_k = True
								print("Consensus reached for k=", str(k), " in ", gen_number, " generations due to negligible cluster point changes.\n",
									  str(best_cluster_group), "\n",sep='') # Delete me
								break
						
						prev_gen_cluster_group_list = cluster_group_list
					
					if break_for_next_k == True:
						# If break_for_next_k == True:
						# gen_number < total_generations. The clusters converged.
						# If break_for_next_k == False:
						# gen_number == total_generations. The clusters did not converge, so create a new initial population and try again.
						break
					
					initial_pop_attempt += 1
					print("Failed to converge within ", str(total_generations), " generations. Starting Try #", str(initial_pop_attempt), " with a new initial population.")
				
				if k > 1 and log(best_cluster_group_for_each_k[k-1][1] * (k-1)) < log(best_cluster_group_for_each_k[k][1] * k):
					# Determine whether or not k-1 dist_sum (normalized for the number of clusters used) was the optimal value for k.
					# If so, halt genetic algorithm and return the k-1 cluster groups to be used for the creation of HMM probabilities.
					k_cluster_count_dict[win][chr] = k-1
					k_loc_dict[win][chr] = best_cluster_group_for_each_k[k-1]
					print("k=", str(k), " is optimal:\n", str(best_cluster_group_for_each_k[k-1]), "\n", sep='')
	
	return(k_cluster_count_dict, k_loc_dict)

def Prepare_HMM(sam, tp_positions, chr_len, chr, win_list, outdir, Rpath):
	# Prepare_HMM creates (1) per-window read count histograms and (2) observation files, as well as
	# (3) details for the HMM parameter file such as read count thresholds and copy number state names.
	#
	# At the end of this function, Dist_To_Params is called to output this data to files, and finally
	# HMM_Dup_Search is called to run the HMM on these files and output the duplication search results.
	
	# Create list of chromosomes to run HMM on
	chr_list = []
	if chr == "All":
		chr_list = [chr_entry for chr_entry in chr_len]
	else:
		chr_list = [chr_entry for chr_entry in chr.split(',')]
	
	# Create the HMM observation files (path_dict, an OrderedDict containings # of reads per window in sequential order)
	# and reads-per-window histograms (hist_dict, a dictionary whose values are chromosome-specific Counters representing reads-per-window)
	hist_dict, path_dict = SAM_to_wincov(sam, win_list, chr_list, outdir, "False")
	
	# Read in histograms if hist_dict is empty for a particular window. This happens when histograms
	# had already been generated for all chromosomes for that window, and is simply a time-saving feature
	# when performing re-runs on the same data set. Histogram and observation files are not re-created if
	# they already exist, as this is the most time-consuming part of the entire pipeline.
	for win in win_list:
		if not hist_dict[win]:
			# If this dictionary entry is empty for some window, grab those histograms from files
			for chr in chr_list:
				hist_i_fn = os.path.join(outdir, str(win) + "bp-window/" + str(chr) + "_hist.txt")
				with open(hist_i_fn) as infile:
					inlines = infile.readlines()
					hist = Counter()
					for line in inlines[1:]:
						cov, freq = line.split("\t")
						hist[int(cov)] = int(freq)
				hist_dict[win][chr] = hist
	
	# Pass histograms to k-means clustering function w/ genetic algorithm-based positioning of groups
	k_cluster_count_dict, k_loc_dict, k_data_count_dict = K_Means_Clustering(hist_dict)

	# Cycle over each window, chromosome, and threshold value combination provided and run HMM for each combination
	for win in win_list:
		for chr in chr_list:
			#for threshold_count in temp_t_count[win][chr]:
			# Create an HMM parameter file (poisson probabilities of given # of reads in a window based on duplication status of that window)
			threshold_count = Dist_To_Params(hist_dict[win][chr], k_cluster_count_dict[win][chr], k_loc_dict[win][chr], k_data_count_dict[win][chr], win, chr, sam, outdir)
			
			# Run HMM using generated parameter file and observation file
			HMM_Dup_Search(outdir, tp_positions, chr, win, threshold_count, k_loc_dict[win][chr], k_data_count_dict[win][chr])

def SAM_to_wincov(sam_file, win_list, chr_list, outdirpart, control_use):
	# Check to see if, for a given window size, all chromosomes have already had their
	# per-window read count histogram generated. If so, do not re-create the existing
	# data files necessary for the HMM to run.
	win_list_needed = win_list
	if control_use == "False":
		win_list_needed = []
		for win in win_list:
			not_needed = 0
			for chr in chr_list:
				# Check here to see if an observation file exists. If so, assume all files have
				# already been generated for this chromosome & window combination
				i_fn = os.path.join(outdirpart, str(win) + "bp-window/" + str(chr) + "_obs.txt")
				if os.path.exists(i_fn): not_needed += 1
			# If we did not see an observation file for all chromosomes to be run on, re-create all
			# the HMM observation and histogram files
			if not_needed != len(chr_list): win_list_needed.append(win)
	
	hist_dict = {win: {} for win in win_list}
	path_dict = {win: {} for win in win_list}
	chr_string = ','.join(chr_list)
	# If list is not empty, create paths for each chromosome/window pair
	# If list is empty, we don't need to re-create existing data files, and thus the following 'if' block is skipped
	if win_list_needed:
		# Get per-window read count histograms and paths for the given window sizes
		hist_dict, path_dict = SAM.Coverage_Window(sam_file, win_list_needed, chr_string, control_use)
		if control_use == "True":
			# Only need path_dict for a file being used only to find the exp:control ratio of % total reads
			# Thus, don't need to output paths or histograms, so skip past the following 'elif' block
			pass
		elif control_use == "False":
			# Check to see if folders have already been made for the provided windows. If not, make them.
			for chr in chr_list:
				for win in win_list:
					outdir = os.path.join(outdirpart,str(win) + "bp-window/")
					print(outdirpart)
					print(outdir)
					if not os.path.exists(outdir): os.makedirs(outdir)
			
			# Write out per-window read count histograms and paths
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

# Grab command line options
sam, gff, chr, win_list, outdir, Rpath = Command_line()

# Determine the names of chromosomes and their lengths
chr_len = SAM.Chr_Lengths(sam)

# If a GFF file is provided, determine the location of all transposons
# This is useful for quickly identifying false positive duplication calls
tp_positions = Grab_Transposons(gff, chr_len)

# Generate the necessary data to perform an HMM-based duplication search,
# then run the HMM on these parameters and observations.
# Detailed commentary is provided within sub-routines
Prepare_HMM(sam, tp_positions, chr_len, chr, win_list, outdir, Rpath)