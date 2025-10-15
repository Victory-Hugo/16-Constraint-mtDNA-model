import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import regional_constraint
import sys
from typing import Dict, IO, List, Tuple
sys.path.append('../mitochondrial_constraint/')
import calculate_oe.oe_functions


def calculate_fdr(rc_final: str, shuffle_final: str, shuffle_list: List[str], locus: str):
	"""Calculate the FDR for each area of regional constraint.

	:param rc_final: path to file with non-overlapping significant regions in real alignment
	:param shuffle_final: path to file with non-overlapping significant regions in all shuffles
	:param shuffle_list: list of unique shuffles to iterate across
	:param locus: locus name
	:return: dictionary with fdr estimates
	"""
	# extract a list of p-value thresholds and lengths to test, per the regions found in the real alignment
	# lengths in bp for protein and non-proteins
	p_list = []
	for row in csv.DictReader(open(rc_final), delimiter='\t'):
		if row["locus"] == locus:
			p_list.append((row["locus"], float(row["kmer_pvalue"]), int(row["length"])))
			
	# now, estimate the fdr for each candidate, using parallel processing
	manager = mp.Manager()
	shuffle_dict = manager.dict()
	pool = mp.Pool(mp.cpu_count())
	for item in p_list:  # list of conditions to count false positives for
		pool.apply_async(parallelize_counting, args=(item, shuffle_list, shuffle_final, shuffle_dict))
	pool.close()
	pool.join()
	fdr = {}
	for item in shuffle_dict:
		fdr[item] = shuffle_dict[item] / len(shuffle_list)
		print("FDR estimate for region with following p-value and length", item, "is", fdr[item])

	return fdr


def parallelize_counting(
		item: Tuple[str, float, int], shuffle_list: List[str], shuffle_final: str,
		shuffle_dict: Dict[Tuple[str, float, int], float]):
	"""Function to enable parallelization of counting of shuffles with a false positives.

	:param item: the locus, p-value threshold and length to assess
	:param shuffle_list: list of unique shuffles to iterate across
	:param shuffle_final: path to file with non-overlapping significant regions in all shuffles
	:param shuffle_dict: dictionary with item as key, and value as the count of shuffles with a false positive
	:return: shuffle dict
	"""
	shuffle_counts = []  # reset for each condition
	for shuffle_num in shuffle_list:  # list of all shuffles to iterate through
		count = 0
		shuffle_file_to_read = open(shuffle_final)
		for row in csv.DictReader(shuffle_file_to_read, delimiter='\t'):
			# if region is in the locus in the loop, and passes p-value and length filters, and in the shuffle in loop
			# length in bp here for all loci
			if (row["locus"] == item[0]) and (float(row["kmer_pvalue"]) <= item[1]) and (int(row["length"]) == item[2]) and (
					row["shuffle"] == shuffle_num.split('-')[0]) and (row["out_dir"].split('/')[-1] == shuffle_num.split('-')[1]):
				count = 1
		shuffle_file_to_read.close()
		shuffle_counts.append(count)
	# now write the total number of shuffles with a false positive for the condition tested to a dictionary
	shuffle_dict[item] = sum(shuffle_counts)
	
	return shuffle_dict


def filter_by_FDR(
		locus: str, significant_kmers: str, shuffle: str, dict: Dict[Tuple[int, int, str], Tuple[float, int, str]],
		pvalue_threshold: float, filters: Dict[int, float]):
	"""Remove candidates which don't pass the FDR threshold and then reapply the greedy algorithm.

	:param locus: locus name
	:param significant_kmers: path to file with all significant (overlapping) regions to filter across
	:param shuffle: shuffle number to filter across
	:param dict: dictionary with key as start and end coordinates of region, and values are a tuple of
	pvalue, length and shuffle number
	:param pvalue_threshold: highest accepted p-value for inclusion for entire analysis
	:param filters: dictionary where key is the length and value the p-value threshold for filtering
	:return: dictionary of output of greedy algorithm, with tuple key of kmer_start, kmer_end, and shuffle, with values a
	tuple of pvalue, length and shuffle number
	"""
	# exclude significant kmers that didn't pass FDR threshold, and repeat greedy algorithm - first read in as dictionary
	outlier_kmers = {}
	for row in csv.DictReader(open(significant_kmers), delimiter='\t'):
		if (row["locus"] == locus) and (float(row["pvalue"]) < pvalue_threshold):  # p-value threshold for entire analysis
			# note length in input file is in codons for protein genes, convert to bp
			length = (int(row["length"]) * 3) if any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) else int(row["length"])
			# if in real alignment or in the shuffle in the loop
			if (shuffle == "0") or (
					(row["shuffle"] == shuffle.split('-')[0]) and (row["out_dir"].split('/')[-1] == shuffle.split('-')[1])):
				include = "yes"
				for filter_key in sorted(filters):  # sort so goes in order of length
					if (length == filter_key) and (float(row["pvalue"]) >= filters[filter_key]):
						include = "no"  # then mark for filtering
				if include == "yes":  # then add to dictionary for input for greedy algorithm
					outlier_kmers[(int(row["kmer_start"]), int(row["kmer_end"]))] = (float(row["pvalue"]), length)
	
	# using parallel processing to apply greedy algorithm
	if len(outlier_kmers) > 0:
		# to more easily handle control region, shift bases 1-576 to following numbering after 16569
		new_dict = {}
		for key in sorted(outlier_kmers):
			key_start = (16569 + key[0]) if (locus == "CR" and key[0] < 577) else key[0]
			key_end = (16569 + key[1]) if (locus == "CR" and key[1] < 577) else key[1]
			new_dict[(key_start, key_end)] = (outlier_kmers[key][0], outlier_kmers[key][1])
		df = pd.DataFrame.from_dict(new_dict, orient='index', columns=['pvalue', 'length'])
		df.index = pd.MultiIndex.from_tuples(df.index, names=['kmer_start', 'kmer_end'])
		df = df.reset_index()
		pool = mp.Pool(mp.cpu_count())
		# this step will return the output of the greedy algorithm, as dict - note the dict is input of function
		for key in list(sorted(new_dict)):  # for every significant kmer
			pool.apply_async(reextract_nonoverlapping, args=(key, df, dict, shuffle))
		pool.close()
		pool.join()
		
	if dict is not None:  # to avoid None assignment
		return dict


def reextract_nonoverlapping(
		key: Tuple[int, int],
		df: pd.DataFrame,
		filtered_outlier_kmers: Dict[Tuple[int, int, str], Tuple[float, int, str]],
		shuffle: str):
	"""Extracts list of non-overlapping kmers by lowest p-value, and length if needed.

	:param key: a tuple key of kmer_start, kmer_end, and shuffle, control region shifted bases 1-576 to numbering >16569
	:param df: pandas dataframe list of all significant kmers
	:param filtered_outlier_kmers: dictionary with tuple key of kmer_start and kmer_end, with values are a tuple of
	pvalue, length and shuffle number
	:param shuffle: shuffle number to extract across
	:return: filtered_outlier_kmers dictionary
	"""
	kmer_start = key[0]
	kmer_end = key[1]
	# subset to all overlapping significant kmers, and sort so that the most significant, longest is top
	# the filter is a precondition for two intervals to overlap
	pre_filter_df = df[
		(df['kmer_start'] <= np.int64(kmer_end)) &
		(df['kmer_end'] >= np.int64(kmer_start))].sort_values(['pvalue', 'length'], ascending=(True, False))
	
	# extract kmer p-value and length
	kmer_pvalue = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))]['pvalue'].iloc[0]
	kmer_length = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))]['length'].iloc[0]
	# then remove the kmer under evaluation from df
	filter_df = pre_filter_df.drop(pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))].index)
		
	# extract non-overlapping using following three conditions
	# (1) if there are no overlapping kmers
	if filter_df.shape[0] == 0:
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)
		# note need shuffle in key to handle if same start, end observed in different shuffle by chance
	
	# (2) if the p-value of the kmer is lower than all overlapping kmers
	elif kmer_pvalue < filter_df['pvalue'].iloc[0]:
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)
	
	# (3) if the p-value is the same as the lowest overlapping p-value
	# and the kmer is longer, or same length, as the overlapping kmer with the same p-value
	elif (kmer_pvalue == filter_df['pvalue'].iloc[0]) and (kmer_length >= filter_df['length'].iloc[0]):
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)
	
	if filtered_outlier_kmers is not None:  # to avoid None assignment
		return filtered_outlier_kmers


def annotate_final(
		key: Tuple[int, int, str], dict_key: Tuple[int, int, str], locus: str, output_file: IO,
		input_file: str, obs_value: str, excluded_sites: List[int], fit_parameters: str,
		pos_to_res: Dict[Tuple[str, int], int],
		real_filtered_outlier_kmers: Dict[Tuple[int, int, str], Tuple[float, int, str]],
		fdr_dict: Dict[Tuple[str, float, int], float]):
	"""Annotate the final regional constraint areas and write to file.

	:param key: a tuple key of kmer_start, kmer_end, and shuffle ("0" for real) in mtDNA coordinates
	:param dict_key: a tuple key of kmer_start, kmer_end, and shuffle ("0" for real) for real_filtered_outlier_kmers
	:param locus: locus name
	:param output_file: output file to write to
	:param input_file: path to annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param excluded_sites: list of base positions to exclude from calculations
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param pos_to_res: dictionary to convert base coordinates to protein residues
	:param real_filtered_outlier_kmers: dictionary with tuple key of kmer_start, kmer_end, and shuffle, with values a
	tuple of pvalue, length and shuffle number
	:param fdr_dict: dictionary with fdr estimates for each regional constraint length and p-value
	"""
	# need these dictionaries and lists for when calculating metrics for final regions
	per_gene = {}
	for row in csv.DictReader(open(
			'required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		# make control region match input
		row["SYMBOL"] = "CR" if ((int(row["POS"]) <= 576) or (int(row["POS"]) >= 16024)) else row["SYMBOL"]
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]

	# list of positions to calculate metrics across
	kmer_coords = list(range(key[0], key[1] + 1)) if (key[0] < key[1]) else \
		list(range(key[0], 16570)) + list(range(1, key[1] + 1))  # handle across m.16569-1

	per_kmer = calculate_oe.oe_functions.initialize_sum_dict(identifier_list=[(key[0], key[1])])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in calculate_oe.oe_functions.ori_region) else 'ref_exc_ori'
			if int(row["POS"]) in kmer_coords:
				# include all bases if RNA gene/non-coding loci, but only missense (most severe) if protein gene
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
					x in row["consequence"] for x in calculate_oe.oe_functions.more_severe_than_missense)):
					per_kmer = calculate_oe.oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=(key[0], key[1]), region=region_to_use,
						observed=row[obs_value], likelihood=row["Likelihood"], dict=per_kmer)
	
	# now calculate observed, expected
	obs_max_het = calculate_oe.oe_functions.calculate_obs(sum_dict=per_kmer, identifier=(key[0], key[1]))
	obs_max_het = 0.0 if obs_max_het < 0.001 else obs_max_het  # to handle very small decimals from processing
	exp_max_het = calculate_oe.oe_functions.calculate_exp(
		sum_dict=per_kmer, identifier=(key[0], key[1]), fit_parameters=fit_parameters)
	total_all = calculate_oe.oe_functions.calculate_total(sum_dict=per_kmer, identifier=(key[0], key[1]))
	ratio_oe = obs_max_het / exp_max_het
	(lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het)
	# to key in residues for outlier_kmers dicts
	output_file.write(
		str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
		str(pos_to_res[(locus, key[0])] if (locus, key[0]) in pos_to_res else '') + '\t' +
		str(pos_to_res[(locus, key[1])] if (locus, key[1]) in pos_to_res else '') + '\t' +
		str(obs_max_het) + '\t' + str(exp_max_het) + '\t' + str(ratio_oe) + '\t' +
		str(lower_CI) + '\t' + str(upper_CI) + '\t' + str(total_all) + '\t' +
		str(real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][0]) + '\t' +
		str(real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][1]) + '\t' +
		str(fdr_dict[(
			locus, real_filtered_outlier_kmers[(
				dict_key[0], dict_key[1], "0")][0], real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][1])]) + '\n')
	output_file.flush()  # to help it be printed directly to file
	

def apply_FDR_filter(
		rc_final: str, shuffle_final: str, fdr_threshold: float, real_all: str, shuffle_all: str,
		shuffle_list: List[str], input_file: str, obs_value: str, excluded_sites: List[int], fit_parameters: str,
		pvalue_threshold: float, loci_type: str, file: IO):
	"""Return a list of high confidence areas of regional constraint.

	:param rc_final: file with non-overlapping significant candidates regional constraint areas in real alignment
	:param shuffle_final: file with non-overlapping significant candidates regional constraint areas in all shuffles
	:param fdr_threshold: fdr threshold
	:param real_all: file with overlapping significant candidates regional constraint areas in real alignment
	:param shuffle_all: file with overlapping significant candidates regional constraint areas in all shuffles
	:param shuffle_list: list of shuffles to iterate through
	:param input_file: path to annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param excluded_sites: list of base positions to exclude from calculations
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param pvalue_threshold: highest accepted p-value for inclusion for entire analysis
	:param loci_type: either protein or rn
	:param file: output file to write results to
	"""
	# process one locus at a time
	loci_list = [
		"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L",
		"MT-ND5", "MT-ND6"] if (loci_type == "protein") else ["MT-RNR1", "MT-RNR2"]
	
	# initiate required dictionaries
	pos_to_res = regional_constraint.pos_to_res()
	res_to_pos = regional_constraint.res_to_pos()
	# record the starting files
	initial_rc_file = rc_final
	initial_shuffle_file = shuffle_final
	manager = mp.Manager()
	
	for locus in loci_list:
		iteration = 0
		fdr_dict = calculate_fdr(
			rc_final=initial_rc_file, shuffle_final=initial_shuffle_file, shuffle_list=shuffle_list, locus=locus)
		
		result = "pass"
		for item in fdr_dict:
			result = "fail" if (fdr_dict[item] >= fdr_threshold) else result
			
		if result == "pass":  # if all candidates for the locus satisfy threshold on first pass
			print(locus, " passed!")
			real_filtered_outlier_kmers = {}
			for row in csv.DictReader(open(initial_rc_file), delimiter='\t'):
				if row["locus"] == locus:
					real_filtered_outlier_kmers[
						(int(row["kmer_start"]), int(row["kmer_end"]), "0")] = (float(row["kmer_pvalue"]), int(row["length"]), "0")
			# write to file
			for key in real_filtered_outlier_kmers:
				annotate_final(
					key=key, dict_key=key, locus=locus, output_file=file, input_file=input_file,
					obs_value=obs_value, excluded_sites=excluded_sites, fit_parameters=fit_parameters,
					pos_to_res=pos_to_res, real_filtered_outlier_kmers=real_filtered_outlier_kmers, fdr_dict=fdr_dict)
		
		# if at least one interval did not pass the FDR threshold, initiate the filters dictionary for locus
		filters = {}
		while result == "fail":  # if at least one candidate in locus failed FDR threshold
			iteration += 1
			print(datetime.datetime.now(), locus, "needs filtering")
			
			# define filters to apply
			for key in fdr_dict:  # key is (row["locus"], float(row["kmer_pvalue"]), int(row["length"]))
				if fdr_dict[key] >= fdr_threshold:
					# create a dictionary of filters where length is key and p-value threshold is the value
					if key[2] not in filters:
						filters[key[2]] = key[1]
					# if the length already has a filter and the p-value is lower than previous entry, replace it
					elif (key[2] in filters) and (key[1] < filters[key[2]]):
						filters[key[2]] = key[1]
						
			# first filter and reapply greedy algorithm for real alignment
			real_filtered_outlier_kmers = manager.dict()
			real_filtered_outlier_kmers = filter_by_FDR(
				locus=locus, significant_kmers=real_all, shuffle="0", dict=real_filtered_outlier_kmers,
				pvalue_threshold=pvalue_threshold, filters=filters)
			
			# then filter and reapply greedy algorithm to the shuffles
			shuffle_filtered_outlier_kmers = manager.dict()
			for shuffle in shuffle_list:
				shuffle_filtered_outlier_kmers = filter_by_FDR(
					locus=locus, significant_kmers=shuffle_all, shuffle=shuffle, dict=shuffle_filtered_outlier_kmers,
					pvalue_threshold=pvalue_threshold, filters=filters)
			
			# generate inputs for FDR estimates by writing output of greedy algorithm to file - length in bp
			rc_final = 'output_files/regional_constraint/temp_regional_constraint_intervals.txt'
			file2 = open(rc_final, "w")
			file2.write("kmer_start	kmer_end	locus	kmer_pvalue	length" + '\n')
			for key in sorted(real_filtered_outlier_kmers):  # sort given asynchronous parallelization
				file2.write(
					str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
					str(real_filtered_outlier_kmers[key][0]) + '\t' +
					str(real_filtered_outlier_kmers[key][1]) + '\t' + '\n')
				file2.flush()  # to help it be printed directly to file
			
			shuffle_final = 'output_files/regional_constraint/temp_shuffle_regional_constraint_intervals.txt'
			file3 = open(shuffle_final, "w")
			file3.write("kmer_start	kmer_end	locus	kmer_pvalue	length	shuffle	out_dir" + '\n')
			for key in sorted(shuffle_filtered_outlier_kmers):
				file3.write(
					str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][0]) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][1]) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][2].split('-')[0]) + '\t' +
					str('out_dir/' + shuffle_filtered_outlier_kmers[key][2].split('-')[1]) + '\n')
				file3.flush()  # to help it be printed directly to file
		
			# now calculate FDR again
			fdr_dict = calculate_fdr(
				rc_final=rc_final, shuffle_final=shuffle_final, shuffle_list=shuffle_list, locus=locus)
			
			result = "pass"
			for item in fdr_dict:
				result = "fail" if (fdr_dict[item] >= fdr_threshold) else result
				
			if result == "pass":  # then calculate obs, exp, etc and write to file
				for key in real_filtered_outlier_kmers:
					# before writing to file will need to do the following steps
					# need to handle control region (as shifted bases 1-576 to following numbering after 16569)
					key_start = (key[0] - 16569) if (key[0] > 16569) else key[0]
					key_end = (key[1] - 16569) if (key[1] > 16569) else key[1]
					key = (key_start, key_end, "0")
					# convert residues to mtDNA coordinates, corresponding to first base of first codon, last base of last codon
					# need this since the outlier dict uses residues not bp for start end
					residue_key = key
					if any(x in locus for x in ["MT-A", "MT-C", "MT-N"]):
						key = (min(res_to_pos[(locus, key[1])]), max(res_to_pos[(locus, key[0])])) if (locus == "MT-ND6") \
							else (min(res_to_pos[(locus, key[0])]), max(res_to_pos[(locus, key[1])]))
					annotate_final(
						key=(key[0], key[1], "0"), dict_key=residue_key, locus=locus, output_file=file, input_file=input_file,
						obs_value=obs_value, excluded_sites=excluded_sites, fit_parameters=fit_parameters,
						pos_to_res=pos_to_res, real_filtered_outlier_kmers=real_filtered_outlier_kmers, fdr_dict=fdr_dict)
				print("All done for locus", locus)
				continue  # to move to next locus
			# if at least one interval does not pass, then this loop will repeat
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	parser.add_argument(
		"-obs", type=str, help="Population dataset from which observed maximum heteroplasmy is obtained")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calculations")
	args = parser.parse_args()
	
	# set default, for gnomAD
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "gnomad_max_hl"
	if args.parameters is None:
		args.parameters = 'output_files/calibration/linear_model_fits.txt'
	if args.exc_sites is None:
		# exclude “artifact_prone_sites” in gnomAD positions - 301, 302, 310, 316, 3107, and 16182 (3107 already excluded)
		# variants at these sites were not called in gnomAD and therefore these positions removed from calculations
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	print(datetime.datetime.now(), "Starting to filter regional constraint by FDR!")
	print("This may take a while...")
	
	# write final output file
	file = open('output_files/regional_constraint/final_regional_constraint_intervals.txt', "w")
	file.write(
		"start	end	locus	protein_pos_start	protein_pos_end	obs_max_het	exp_max_het	ratio_oe	"
		"lower_CI	upper_CI	total	pvalue	length	fdr" + '\n')
	
	# assemble list of unique shuffles to use - first for protein
	number_shuffles = 70
	number_out_dir = 14
	shuffle_list = []
	for shuffle_number in range(1, number_shuffles + 1):
		for out_dir in range(1, number_out_dir + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# manually add last shuffle to equal 1000 (out directory 15)
	for shuffle_number in range(1, 20 + 1):
		shuffle_list.append(str(shuffle_number) + '-' + str(15))
	
	apply_FDR_filter(
		rc_final='output_files/regional_constraint/real_alignment/regional_constraint_intervals.txt',
		shuffle_final='output_files/regional_constraint/protein_shuffle/all_shuffles_nonoverlapping.txt',
		fdr_threshold=0.10,
		real_all='output_files/regional_constraint/real_alignment/all_significant_kmers.txt',
		shuffle_all='output_files/regional_constraint/protein_shuffle/all_sig_shuffles.txt',
		shuffle_list=shuffle_list,
		input_file=args.input, obs_value=args.obs, excluded_sites=args.exc_sites, fit_parameters=args.parameters,
		pvalue_threshold=0.005, loci_type="protein", file=file)
	# pvalue_threshold select based on highest p-value across candidates from real alignment
	
	# assemble list of unique shuffles to use - now for RNA
	# note this was done in three batches
	number_shuffles = 19
	number_out_dir = 20
	shuffle_list = []
	for shuffle_number in range(1, number_shuffles + 1):
		for out_dir in range(1, number_out_dir + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# batch 2
	for shuffle_number in range(1, 10 + 1):
		for out_dir in range(21, 32 + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# batch 3
	for shuffle_number in range(1, 10 + 1):
		for out_dir in range(33, 82 + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	
	print("Calculating for", len(shuffle_list), "shuffles")
			
	apply_FDR_filter(
		rc_final='output_files/regional_constraint/real_alignment/regional_constraint_intervals.txt',
		shuffle_final='output_files/regional_constraint/rn_shuffle/all_shuffles_nonoverlapping.txt',
		fdr_threshold=0.10,
		real_all='output_files/regional_constraint/real_alignment/all_significant_kmers.txt',
		shuffle_all='output_files/regional_constraint/rn_shuffle/all_sig_shuffles.txt',
		shuffle_list=shuffle_list,
		input_file=args.input, obs_value=args.obs, excluded_sites=args.exc_sites, fit_parameters=args.parameters,
		pvalue_threshold=0.005, loci_type="rn", file=file)
	
	# delete temporary file outputs	once complete
	os.remove('output_files/regional_constraint/temp_regional_constraint_intervals.txt')
	os.remove('output_files/regional_constraint/temp_shuffle_regional_constraint_intervals.txt')
	
	print(datetime.datetime.now(), "Script finished!")
	