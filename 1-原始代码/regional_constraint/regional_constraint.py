import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import sys
from typing import Dict, List, Tuple, Union
sys.path.append('../mitochondrial_constraint/')
import calculate_oe.oe_functions


def res_to_pos():
	"""Create dictionary to convert protein residue to list of base positions encoding said residue.

	:return: res_to_pos dictionary with tuple key of locus and residue number, with values as list of base positions
	"""
	dict = {}
	for row in csv.DictReader(
			open("required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		if row["Protein_position"]:
			if (row["SYMBOL"], int(row["Protein_position"])) not in dict:
				dict[(row["SYMBOL"], int(row["Protein_position"]))] = [int(row["POS"])]
			else:
				dict[(row["SYMBOL"], int(row["Protein_position"]))].append(int(row["POS"]))
	return dict
	
	
def pos_to_res():
	"""Create dictionary to convert a base position to protein residue number.

	:return: pos_to_res dictionary with tuple key of locus and base position, with residue number value
	"""
	dict = {}
	for row in csv.DictReader(
			open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		if row["Protein_position"] and ((row["SYMBOL"], int(row["POS"])) not in dict):
			dict[(row["SYMBOL"], int(row["POS"]))] = int(row["Protein_position"])
	return dict


def calculate_per_pos(
		locus: str,
		locus_coords: List[int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int]):
	"""Sum the observed maximum heteroplasmy, mutation likelihoods and count for each base position.

	:param locus: locus name
	:param locus_coords: list of base coordinates encoding the locus
	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param excluded_sites: list of base positions to exclude from calculations
	:return: per_pos dictionary with tuple key of the base position, mutation group, value to sum, and region with variant,
	and value corresponding to the value to sum in key
	"""
	# to handle variants in two genes
	per_gene = {}
	for row in csv.DictReader(open(
			'required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		# make control region match input
		row["SYMBOL"] = "CR" if ((int(row["POS"]) <= 576) or (int(row["POS"]) >= 16024)) else row["SYMBOL"]
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]
	
	# convert locus coordinates to list of str
	per_pos = calculate_oe.oe_functions.initialize_sum_dict(identifier_list=[str(x) for x in locus_coords])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in calculate_oe.oe_functions.ori_region) else 'ref_exc_ori'
			if row["POS"] in [str(x) for x in locus_coords]:
				# include all bases if RNA gene/non-coding loci, but only missense (most severe) if protein gene
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
					x in row["consequence"] for x in calculate_oe.oe_functions.more_severe_than_missense)):
					per_pos = calculate_oe.oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=row["POS"], region=region_to_use,
						observed=row[obs_value], likelihood=row["Likelihood"], dict=per_pos)
	return per_pos


def calculate_per_res(
		locus: str,
		locus_coords: List[int],
		codon_coords: List[int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int],
		pos_to_res: Dict[Tuple[str, int], int]):
	"""Sum the observed maximum heteroplasmy, mutation likelihoods and count for each protein residue.

	:param locus: locus name
	:param locus_coords: list of base coordinates encoding the locus
	:param codon_coords: list of protein residues, starting from 1
	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param excluded_sites: list of base positions to exclude from calculations
	:param pos_to_res: dictionary to convert base coordinates to protein residues
	:return: per_res dictionary with a tuple key of the residue, mutation group, value to sum, and region with variant,
	and value corresponding to the value to sum in key
	"""
	# to handle variants in two genes
	per_gene = {}
	for row in csv.DictReader(open(
			'required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]
	
	# convert coords to list of str
	per_res = calculate_oe.oe_functions.initialize_sum_dict(identifier_list=[str(x) for x in codon_coords])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in calculate_oe.oe_functions.ori_region) else 'ref_exc_ori'
			if row["POS"] in [str(x) for x in locus_coords]:
				# include all bases if RNA gene/non-coding loci, but only missense (most severe) if protein gene
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
					x in row["consequence"] for x in calculate_oe.oe_functions.more_severe_than_missense)):
					# note need to use pos_to_res to handle positions in two genes
					per_res = calculate_oe.oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=str(pos_to_res[(locus, int(row["POS"]))]),
						region=region_to_use, observed=row[obs_value], likelihood=row["Likelihood"], dict=per_res)
	return per_res


def shuffle(
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		coords: List[str]):
	"""Function that creates shuffles, or random permutations, of the real alignment.
	Randomly reassign the observed, mutation likelihoods and total variants of a position to another position.

	:param per_unit_dict: per_pos for non-protein, per_res for proteins
	:param coords: list of base coordinates or protein residue numbers for shuffle
	:return: shuffled dictionary with a tuple key of the base position, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	"""
	# first, change format of the per_unit_dict for easy shuffling
	shuffle_dict = {}
	col_names = []  # to preserve order of appendage, same for all positions
	for key in per_unit_dict:
		pos = key[0]
		if pos not in shuffle_dict:
			shuffle_dict[pos] = [per_unit_dict[key]]
			if (str(key[1]) + '-' + str(key[2]) + '-' + str(key[3])) not in col_names:
				col_names.append(str(key[1]) + '-' + str(key[2]) + '-' + str(key[3]))
		else:
			shuffle_dict[pos].append(per_unit_dict[key])
			if (str(key[1]) + '-' + str(key[2]) + '-' + str(key[3])) not in col_names:
				col_names.append(str(key[1]) + '-' + str(key[2]) + '-' + str(key[3]))
	print("after first reformat", shuffle_dict, col_names)

	df = pd.DataFrame.from_dict(shuffle_dict, orient='index', columns=col_names)
	pd.set_option('display.max_columns', None)
	pd.set_option('display.max_rows', None)
	print("before shuffle", df)

	df = df.sample(frac=1).reset_index()  # sample all rows without replacement and return all, reset index
	print("after shuffle", df)

	df.set_axis(coords, axis='index', inplace=True)
	print("renumber index", df)

	df = df.drop(columns=['index'])  # drop the column with the original index
	shuffle_dict = df.T.to_dict(orient='list')  # to set keys as index
	print("new shuffled", shuffle_dict)
	pd.set_option('display.max_rows', 5)

	# now convert back to input format
	per_unit_dict = {}
	# column is the identifiers in the key, index to extract order in values in the shuffle_dict
	for index, column in enumerate(col_names):
		for key in shuffle_dict:  # key is the position
			if column.split('-')[1] == "count":  # needs to be int
				per_unit_dict[(
					str(key), column.split('-')[0], column.split('-')[1], column.split('-')[2])] = int(shuffle_dict[key][index])
			else:  # these are floats
				per_unit_dict[(
					str(key), column.split('-')[0], column.split('-')[1], column.split('-')[2])] = float(shuffle_dict[key][index])
	print("after second reformat", per_unit_dict)

	return per_unit_dict

	
def calculate_kmers_step1(
		pos: int,
		kmer_length: int,
		coords: List[int],
		locus_oe: float,
		locus_lower_CI: float,
		fit_parameters: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		outlier_kmers: Dict[Tuple[int, int], Tuple[float, int]],
		sig_threshold: float):
	"""Calculate obs:exp and p-value for all possible kmers that start at position 1 in the locus.

	:param pos: start position of the kmer, either as base coordinate or residue
	:param kmer_length: the length of the kmer to evaluate
	:param coords: list of coordinates in locus, codon_coords for proteins and locus_coords for non-proteins
	:param locus_oe: the observed:expected ratio for the locus
	:param locus_lower_CI: the lower bound of the locus 90% CI around the observed:expected ratio
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param per_unit_dict: dictionary with a tuple key of the base position or residue, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param kmers_dict: dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param outlier_kmers: dictionary of significant kmers, keys are start, end coordinates, values are p-value and length
	:param sig_threshold: the p-value threshold to apply, only kmers with p<threshold are retained
	:return: kmers_dict dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:return: outlier_kmers dictionary of significant kmers, keys are start, end coordinates, values are p-value and length
	"""
	kmer_start = pos
	kmer_end = pos + kmer_length - 1  # -1 so hits end of gene
	
	if kmer_end > 16569:  # to handle circular genome
		kmer_range = list(range(kmer_start, 16570)) + list(range(1, (kmer_end - 16569 + 1)))
		kmer_end = int(kmer_range[-1])
	else:
		kmer_range = list(range(kmer_start, kmer_end + 1))  # all positions in the kmer

	if kmer_end in coords:  # to bound by end of locus
		# first initiate the entry in the dictionary
		for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to fit
			for value in ['obs_max_het', 'sum_LR', 'count']:  # the three values to sum
				for region in ['ref_exc_ori', 'ori']:  # the two mutational models that could be applied
					kmers_dict[(kmer_start, kmer_end), mut_group, value, region] = 0
		for pos in kmer_range:
			for key in per_unit_dict:
				if int(key[0]) == pos:
					kmers_dict[((kmer_start, kmer_end), key[1], key[2], key[3])] += per_unit_dict[key]
		
		# now calculate p-value, first calculate observed, expected
		exp_max_het = calculate_oe.oe_functions.calculate_exp(
			sum_dict=kmers_dict, identifier=(kmer_start, kmer_end), fit_parameters=fit_parameters)
		obs_max_het = calculate_oe.oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		obs_max_het = 0.0 if obs_max_het < 0.001 else obs_max_het  # to handle very small decimals from processing
		total_all = calculate_oe.oe_functions.calculate_total(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		# note expected here is calculated based on the locus obs:exp ratio
		less_than_pvalue = calculate_oe.oe_functions.calculate_pvalue(
			obs_max_het=obs_max_het, total=total_all, exp_max_het=(exp_max_het * locus_oe))
		
		# print((kmer_start, kmer_end), obs_max_het, exp_max_het, total_all, less_than_pvalue)
		
		if float(less_than_pvalue) < sig_threshold:  # only consider kmers below this p-value
			(lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_CI(
				obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=2.5)
			if upper_CI < locus_lower_CI:  # only consider if CI does not overlap loci CI
				outlier_kmers[(kmer_start, kmer_end)] = (less_than_pvalue, kmer_length)

	return kmers_dict, outlier_kmers


def calculate_kmers_step2(
		pos: int,
		kmer_length: int,
		locus: str,
		coords: List[int],
		fit_parameters: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		outlier_kmers: Dict[Tuple[int, int], Tuple[float, int]],
		locus_oe: float,
		locus_lower_CI: float,
		sig_threshold: float):
	"""Calculate obs:exp and p-value for all possible kmers that do not start at position 1 in the locus.

	:param pos: start position of the kmer, either as base coordinate or residue
	:param kmer_length: the length of the kmer to evaluate
	:param locus: name of locus
	:param coords: list of coordinates in locus, codon_coords for proteins and locus_coords for non-proteins
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param per_unit_dict: dictionary with a tuple key of the base position or residue, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param kmers_dict: dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param outlier_kmers: dictionary of significant kmers, keys are start, end coordinates, values are p-value and length
	:param locus_oe: the observed:expected ratio for the locus
	:param locus_lower_CI: the lower bound of the locus 90% CI around the observed:expected ratio
	:param sig_threshold: the p-value threshold to apply, only kmers with p<threshold are retained
	:return: kmers_dict dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:return: outlier_kmers dictionary of significant kmers, keys are start, end coordinates, values are p-value and length
	"""
	kmer_start = pos
	kmer_end = pos + kmer_length - 1
	# now use kmers_dict from step 1 to quickly calculate
	# for example, the values for kmer starting at 2 and ending at 10 will be equal to
	# the values for the kmer that starts at 1 and ends at 10 MINUS the values for position 1
	kmer_start_minus_one = 16569 if (locus == "CR" and pos == 1) else (pos - 1)
	
	if kmer_end > 16569:  # to handle circular genome
		k_range = list(range(kmer_start, 16570)) + list(range(1, (kmer_end - 16569 + 1)))
		kmer_end = int(k_range[-1])
		
	if kmer_end in coords:  # to bound by locus end
		# add entry to dictionary, at same time extract list of keys
		keys = []
		for identifier in [(kmer_start, kmer_end)]:
			for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to fit
				for value in ['obs_max_het', 'sum_LR', 'count']:  # the three values to sum
					for region in ['ref_exc_ori', 'ori']:  # the two mutational models that could be applied
						kmers_dict[identifier, mut_group, value, region] = 0
						keys.append(mut_group + '-' + value + '-' + region)
		
		# calculate the values for the kmer by taking the kmer values that are one position longer and then
		# subtracting the values for that position alone
		for key in keys:
			kmers_dict[((kmer_start, kmer_end), key.split('-')[0], key.split('-')[1], key.split('-')[2])] += \
				(kmers_dict[(
					(kmer_start_minus_one, kmer_end), key.split('-')[0], key.split('-')[1], key.split('-')[2])] - per_unit_dict[
					(str(kmer_start_minus_one), key.split('-')[0], key.split('-')[1], key.split('-')[2])])
		
		# now calculate p-value, first calculate observed, expected
		exp_max_het = calculate_oe.oe_functions.calculate_exp(
			sum_dict=kmers_dict, identifier=(kmer_start, kmer_end), fit_parameters=fit_parameters)
		obs_max_het = calculate_oe.oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		obs_max_het = 0.0 if obs_max_het < 0.001 else obs_max_het  # to handle very small decimals from processing
		total_all = calculate_oe.oe_functions.calculate_total(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		# note expected here is calculated based on the locus obs:exp ratio
		less_than_pvalue = calculate_oe.oe_functions.calculate_pvalue(
			obs_max_het=obs_max_het, total=total_all, exp_max_het=(exp_max_het * locus_oe))
		
		# print(locus, (kmer_start, kmer_end), obs_max_het, exp_max_het, total_all, less_than_pvalue)
		
		if float(less_than_pvalue) < sig_threshold:  # only consider kmers below this p-value
			(lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_CI(
				obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=2.5)
			if upper_CI < locus_lower_CI:  # only consider if CI does not overlap loci CI
				outlier_kmers[(kmer_start, kmer_end)] = (less_than_pvalue, kmer_length)

	return kmers_dict, outlier_kmers


def extract_nonoverlapping_kmers(
		key: Tuple[int, int],
		df: pd.DataFrame,
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		locus: str,
		fit_parameters: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		filtered_outlier_kmers: Dict[Tuple[int, int], Tuple[float, float, float, float, float, int, float, int]]):
	"""Extracts list of non-overlapping kmers by lowest p-value, and length if needed.

	:param key: a tuple key of kmer_start, kmer_end, with control region shifted bases 1-576 to numbering after 16569
	:param df: pandas dataframe list of all significant kmers
	:param kmers_dict: dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param locus: name of locus
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param res_to_pos: dictionary to convert protein residues to base coordinates
	:param filtered_outlier_kmers: dictionary with tuple key of kmer_start and kmer_end, with values of interest
	:return: filtered_outlier_kmers dictionary with tuple key of kmer_start and kmer_end, with values of interest
	"""
	kmer_start = key[0]
	kmer_end = key[1]
	# subset to all overlapping significant kmers, and sort so that the most significant, longest is top
	# the filter is a precondition for two intervals to overlap
	pre_filter_df = df[
		(df['kmer_start'] <= np.int64(kmer_end)) &
		(df['kmer_end'] >= np.int64(kmer_start))].sort_values(['less_than_pvalue', 'length'], ascending=(True, False))

	# extract kmer p-value and length
	kmer_pvalue = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
				pre_filter_df['kmer_end'] == np.int64(kmer_end))]['less_than_pvalue'].iloc[0]
	kmer_length = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
				pre_filter_df['kmer_end'] == np.int64(kmer_end))]['length'].iloc[0]
	# then remove the kmer under evaluation from df
	filter_df = pre_filter_df.drop(pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))].index)
	
	# extract non-overlapping using following three conditions
	# (1) if there are no overlapping kmers
	if filter_df.shape[0] == 0:
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 1", key, "to final dict")
		
	# (2) if the p-value of the kmer is lower than all overlapping kmers
	elif kmer_pvalue < filter_df['less_than_pvalue'].iloc[0]:
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 2", key, "to final dict")
	
	# (3) if the p-value is the same as the lowest overlapping p-value
	# and the kmer is longer, or same length, as the overlapping kmer with the same p-value
	elif (kmer_pvalue == filter_df['less_than_pvalue'].iloc[0]) and (kmer_length >= filter_df['length'].iloc[0]):
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 3", key, "to final dict")

	return filtered_outlier_kmers


def annotate_nonoverlapping(
		key: Tuple[int, int],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		fit_parameters: str,
		locus: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		kmer_pvalue: float,
		kmer_length: int,
		filtered_outlier_kmers: Dict[Tuple[int, int], Tuple[float, float, float, float, float, int, float, int]]):
	"""Annotate a list of non-overlapping kmers for writing to file.
	
	:param key: a tuple key of kmer_start, kmer_end, with control region shifted bases 1-576 to numbering after 16569
	:param kmers_dict: dictionary with a tuple key of the kmer start and end, mutation group, value to sum,
	and region with variant, and value corresponding to the value to sum in key
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param locus: name of locus
	:param res_to_pos: dictionary to convert protein residues to base coordinates
	:param kmer_pvalue: p-value for the kmer
	:param kmer_length: length of the kmer, residues for protein or base pairs for non-proteins
	:param filtered_outlier_kmers: dictionary with tuple key of kmer_start and kmer_end, with values of interest
	:return: filtered_outlier_kmers dictionary with tuple key of kmer_start and kmer_end, with values of interest
	"""
	# need to handle control region (as shifted bases 1-576 to following numbering after 16569)
	key_start = (key[0] - 16569) if (key[0] > 16569) else key[0]
	key_end = (key[1] - 16569) if (key[1] > 16569) else key[1]
	key = (key_start, key_end)
	
	# now calculate obs:exp ratio and CI
	obs_max_het = calculate_oe.oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=key)
	obs_max_het = 0 if obs_max_het < 0.001 else obs_max_het  # to handle very small decimals from processing
	exp_max_het = calculate_oe.oe_functions.calculate_exp(
		sum_dict=kmers_dict, identifier=key, fit_parameters=fit_parameters)
	ratio_oe = obs_max_het / exp_max_het
	total_all = calculate_oe.oe_functions.calculate_total(sum_dict=kmers_dict, identifier=key)
	(lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=2.5)
	
	# convert residues to mtDNA coordinates, corresponding to first base of first codon, last base of last codon
	# for MT-ND6 on reverse strand, this corresponds to last base of last codon, first base of first codon
	if any(x in locus for x in protein_genes):
		key = (min(res_to_pos[(locus, key[1])]), max(res_to_pos[(locus, key[0])])) if (locus == "MT-ND6") \
			else (min(res_to_pos[(locus, key[0])]), max(res_to_pos[(locus, key[1])]))
		kmer_length = kmer_length * 3  # to convert to base pair length
	
	# save to dictionary, now all in base coordinates
	filtered_outlier_kmers[key] = (
		obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI, total_all, kmer_pvalue, kmer_length)
	
	return filtered_outlier_kmers


def check_locus_oe(
		input: str,
		input_noncoding: str,
		locus: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		fit_parameters: str):
	"""Calculate the obs:exp ratio for each locus, this also serves as a code check.

	:param input: path to file with loci obs:exp values
	:param input_noncoding: path to file with noncoding loci obs:exp values
	:param locus: name of locus
	:param per_unit_dict: per_res for proteins and per_pos for non-proteins
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:return: locus_oe, the obs:exp ratio for the locus, and the lower CI of the gene 90% CI
	"""
	locus_oe = ''
	for file in (input, input_noncoding):
		for row in csv.DictReader(open(file), delimiter='\t'):
			row["locus"] = "CR" if (row["locus"] == "MT-CR") else row["locus"]  # for ease of parsing later
			if (row["consequence"] == "SNV" or row["consequence"] == "missense") and (row["locus"] == locus):
				locus_oe = float(row["obs:exp"])

	# calculate locus obs:exp here, as a check, using per_unit_dict
	loci_dict = calculate_oe.oe_functions.initialize_sum_dict(identifier_list=[locus])
	for key in per_unit_dict:
		loci_dict[(locus, key[1], key[2], key[3])] += per_unit_dict[key]

	# now calculate observed:expected ratio
	obs_max_het = calculate_oe.oe_functions.calculate_obs(sum_dict=loci_dict, identifier=locus)
	exp_max_het = calculate_oe.oe_functions.calculate_exp(
		sum_dict=loci_dict, identifier=locus, fit_parameters=fit_parameters)
	total_all = calculate_oe.oe_functions.calculate_total(sum_dict=loci_dict, identifier=locus)
	(lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=2.5)
	
	if round(locus_oe, 2) != round(obs_max_het / exp_max_het, 2):  # round as slight difference many decimal places down
		print("Error: gene obs:exp does not equal expected value", locus_oe, "vs", obs_max_het / exp_max_het)
		sys.exit()
	
	return locus_oe, lower_CI


def regional_constraint(
		out_dir: str,
		n_shuffles: int,
		loci_type: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		pos_to_res: Dict[Tuple[str, int], int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int],
		fit_parameters: str,
		loci_input: str,
		noncoding_input: str,
		sig_threshold: float):
	"""Identify a non-overlapping set of regions in a gene that are more significantly constrained than the gene.

	:param out_dir: the directory to save the output files to
	:param n_shuffles: the number of shuffles to complete - if 0, then don't shuffle (do real alignment)
	:param loci_type: options are 'p' for protein, 'rn' for RNA and non-coding, or 'both'
	:param res_to_pos: dictionary to convert protein residues to base coordinates
	:param pos_to_res: dictionary to convert base coordinates to protein residues
	:param input_file: path to annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param excluded_sites: list of base positions to exclude from calculations
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param loci_input: path to file with loci obs:exp values
	:param noncoding_input: path to file with noncoding loci obs:exp values
	:param sig_threshold: the p-value threshold to apply, only kmers with p < threshold are retained
	"""
	# make output files
	out_dir = "real_alignment" if n_shuffles == 0 else out_dir
	out_dir = 'output_files/regional_constraint/' + out_dir
	kmer_file = open('%s/all_significant_kmers.txt' % out_dir, "w")
	kmer_file.write("kmer_start	kmer_end	locus	pvalue	length	shuffle" + '\n')

	kmer_file2 = open('%s/regional_constraint_intervals.txt' % out_dir, "w")
	header = "kmer_start	kmer_end	locus	protein_pos_start	protein_pos_end	obs_max_het	exp_max_het	ratio_oe	lower_CI	upper_CI	total	kmer_pvalue	length	shuffle	out_dir"
	kmer_file2.write(header + '\n')

	# n_shuffles = 0 is REAL alignment
	number_shuffles = [0] if (int(n_shuffles) == 0) else list(range(1, n_shuffles + 1))
	
	for shuffle_num in number_shuffles:
		for row in csv.DictReader(open('required_files/databases/mitomap_genome_loci.txt'), delimiter='\t'):
			row["Map_Locus"] = row["Map_Locus"] if (row["Map_Locus"] != "MT-CR") else "CR"  # for ease of parsing later
			if any(x in row["Map_Locus"] for x in loci_list):
				if (loci_type == "p") and not any(x in row["Map_Locus"] for x in protein_genes):
					continue  # to skip this locus
				elif (loci_type == "rn") and any(x in row["Map_Locus"] for x in protein_genes):
					continue  # to skip this locus
	
				# SET PARAMETERS, GENERATE INPUTS #
				print("For locus", row["Map_Locus"])
				locus_start = int(row["Starting"])
				locus_end = int(row["Ending"])
				locus = row["Map_Locus"]
				locus_coords = list(range(locus_start, locus_end + 1)) if (locus_start < locus_end) else \
					list(range(locus_start, 16570)) + list(range(1, locus_end + 1))  # handle across m.16569-1
				max_kmer_length = len(locus_coords)
				# if gene is a protein gene, generate a list of protein residues, rewrite maximum length
				if any(x in locus for x in protein_genes):
					codon_coords = list(range(1, pos_to_res[(locus, locus_start)] + 1)) if (locus == "MT-ND6") else \
						list(range(1, pos_to_res[(locus, locus_end)] + 1))
					max_kmer_length = len(codon_coords)
				else:
					codon_coords = None
				# protein kmers of same size have about half as much expected, using same protein minimum as for local constraint
				# list of all possible kmer sizes to try, max is length of locus
				kmer_lengths = list(range(10, max_kmer_length + 1)) if any(x in locus for x in protein_genes) else \
					list(range(20, max_kmer_length + 1))
				
				# make per_pos dictionary
				per_pos = calculate_per_pos(
					locus=locus, locus_coords=locus_coords, input_file=input_file,
					obs_value=obs_value, excluded_sites=excluded_sites)
				# calculate per residue if protein gene
				per_res = calculate_per_res(
					locus=locus, locus_coords=locus_coords, codon_coords=codon_coords,
					input_file=input_file, obs_value=obs_value, excluded_sites=excluded_sites,
					pos_to_res=pos_to_res) if any(x in locus for x in protein_genes) else None
				# make per_unit_dictionary, will use per_res if protein and per_pos if not protein
				per_unit_dict = per_res if any(x in locus for x in protein_genes) else per_pos
				coords = codon_coords if any(x in locus for x in protein_genes) else locus_coords
				
				# create permuted sequence, if specified
				if shuffle_num > 0:
					per_unit_dict = shuffle(per_unit_dict=per_unit_dict, coords=coords)
					print(datetime.datetime.now(), "Shuffle number", shuffle_num, "for locus", locus, "has completed!")
				
				# calculate locus obs:exp as needed for regional constraint analysis, also as a check
				(locus_oe, locus_lower_CI) = check_locus_oe(
					input=loci_input, input_noncoding=noncoding_input, locus=locus, per_unit_dict=per_unit_dict,
					fit_parameters=fit_parameters)
	
				# REGIONAL CONSTRAINT #
				# using parallel processing
				manager = mp.Manager()
				kmers_dict = manager.dict()  # new each iteration/shuffle
				outlier_kmers = manager.dict()
	
				# calculate obs:exp and p-values for all possible kmers in gene/locus - use loop to parallelize
				# these functions return kmers_dict and outlier_kmers dictionaries
				for pos in coords:  # for every base or residue in the locus
					pool = mp.Pool(mp.cpu_count())
					if pos == coords[0]:  # if the first base position/residue in the locus
						for kmer_length in kmer_lengths:  # for every possible kmer length to be evaluated
							pool.apply_async(calculate_kmers_step1, args=(
								pos, kmer_length, coords, locus_oe, locus_lower_CI, fit_parameters, per_unit_dict, kmers_dict, outlier_kmers, sig_threshold))
					else:  # not first position of locus
						for kmer_length in kmer_lengths:
							pool.apply_async(calculate_kmers_step2, args=(
								pos, kmer_length, locus, coords, fit_parameters, per_unit_dict, kmers_dict, outlier_kmers, locus_oe, locus_lower_CI, sig_threshold))
					pool.close()
					pool.join()  # postpones the execution of next line of code until all processes in the queue are done.
	
				print(datetime.datetime.now(), "Shuffle number", shuffle_num, "for locus", locus, "has all kmers calculated!")
	
				# print the results for all significant kmers
				for key in sorted(outlier_kmers):  # sort given asynchronous parallelization
					kmer_file.write(
							str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' + str(outlier_kmers[key][0])
							+ '\t' + str(outlier_kmers[key][1]) + '\t' + str(shuffle_num) + '\n')
					kmer_file.flush()  # to help it be printed directly to file
		
				# now select for significant kmers which DO NOT have an overlapping kmer with a lower p-value
				if len(outlier_kmers) > 0:
					# to more easily handle control region, shift bases 1-576 to following numbering after 16569
					filtered_outlier_kmers = manager.dict()
					new_dict = {}
					for key in sorted(outlier_kmers):
						key_start = (16569 + key[0]) if (locus == "CR" and key[0] < 577) else key[0]
						key_end = (16569 + key[1]) if (locus == "CR" and key[1] < 577) else key[1]
						new_dict[(key_start, key_end)] = (outlier_kmers[key][0], outlier_kmers[key][1])
					# df is made from the base_counts made to handle the control region
					df = pd.DataFrame.from_dict(new_dict, orient='index', columns=['less_than_pvalue', 'length'])
					df.index = pd.MultiIndex.from_tuples(df.index, names=['kmer_start', 'kmer_end'])
					df = df.reset_index()
					pool = mp.Pool(mp.cpu_count())
					for key in list(sorted(new_dict)):  # for every significant kmer
						pool.apply_async(extract_nonoverlapping_kmers, args=(
							key, df, kmers_dict, locus, fit_parameters, res_to_pos, filtered_outlier_kmers))
					pool.close()
					pool.join()
					
					for key in sorted(filtered_outlier_kmers):  # sort given asynchronous parallelization
						kmer_file2.write(
							str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
							str(pos_to_res[(locus, key[0])] if (locus, key[0]) in pos_to_res else '') + '\t' +
							str(pos_to_res[(locus, key[1])] if (locus, key[1]) in pos_to_res else '') + '\t' +
							'\t'.join(str(x) for x in filtered_outlier_kmers[key]) + '\t' +
							str(shuffle_num) + '\t' + str(out_dir) + '\n')
						kmer_file2.flush()  # to help it be printed directly to file
					
					print("For loci", row["Map_Locus"], "the non-overlapping regions have been extracted!")
	

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
	parser.add_argument(
		"-n_shuffles", type=int, help="Number of shuffled alignments to generate, optional")
	parser.add_argument(
		"-loci_type", type=str, help="Options are 'p' for protein, 'rn' for RNA and non-coding, or 'both' for both")
	parser.add_argument(
		"-loci_input", type=str, help="Path to file with locus obs:exp values")
	parser.add_argument(
		"-noncoding_input", type=str, help="Path to file with noncoding loci obs:exp values")
	parser.add_argument(
		"-out_dir", type=str, help="Name of the output directory")
	parser.add_argument(
		"-sig_threshold", type=float, help="P-value threshold to apply")
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
	if args.n_shuffles is None:
		args.n_shuffles = 0  # real alignment
	if args.loci_input is None:
		args.loci_input = 'output_files/oe/genes_obs_exp.txt'
	if args.noncoding_input is None:
		args.noncoding_input = 'output_files/oe/noncoding_obs_exp.txt'
	
	for path in ['output_files/regional_constraint/%s' % args.out_dir]:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Starting to analyze regional constraint!")
	
	# loci to evaluate
	loci_list = ["MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-RNR1", "MT-RNR2"]
	protein_genes = ["MT-A", "MT-C", "MT-N"]
	
	res_to_pos = res_to_pos()
	pos_to_res = pos_to_res()

	regional_constraint(
		out_dir=args.out_dir, n_shuffles=args.n_shuffles, loci_type=args.loci_type, input_file=args.input,
		obs_value=args.obs, excluded_sites=args.exc_sites,
		fit_parameters=args.parameters, loci_input=args.loci_input, noncoding_input=args.noncoding_input,
		res_to_pos=res_to_pos, pos_to_res=pos_to_res, sig_threshold=args.sig_threshold)

	print(datetime.datetime.now(), "Script finished!")
