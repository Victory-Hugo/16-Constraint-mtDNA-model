import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
from numpy import mean
import os
import pandas as pd
import sys
from typing import Dict, List, Tuple
sys.path.append('../mitochondrial_constraint/')
import calculate_oe.oe_functions
import build_model.annotate_mutations
import build_model.compile_denovo


def lookup_dict(input_file: str, obs_value: str):
	"""Create a dictionary that can be used to look up observed maximum heteroplasmy and likelihood values.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:return: dictionary where the key is a tuple of position and mutation, and value is tuple of region with variant, and
	observed and likelihood values
	"""
	dict = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		# restrict to base and residue substitutions, ie missense where most severe, RNA and non-coding variants only
		if (any(x in row["consequence"] for x in ["missense", "non_coding_transcript", "intergenic"])) and not any(
				x in row["consequence"] for x in calculate_oe.oe_functions.more_severe_than_missense):
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in calculate_oe.oe_functions.ori_region) else 'ref_exc_ori'
			observed = float(row[obs_value])
			likelihood = float(row["Likelihood"])
			dict[(row["POS"], mutation)] = (region_to_use, observed, likelihood)
	return dict


def parallelize_kmers(
	position: int, 
	kmer_length: int,
	lookup_dictionary: Dict[Tuple[str, str], Tuple[str, float, float]],
	excluded_sites: List[int],
	fit_parameters: str,
	loci_dict: Dict[Tuple[int, int], Tuple[str, str]],
	kmer_dict: Dict[Tuple[int, int], Tuple[int, float, float, float, float, float, str, str]]):
	"""Function to parallize the calculation of the obs:exp ratio and 90% confidence interval of every kmer.

	:param position: the start position of the kmer
	:param kmer_length: the user-specified length of every kmer
	:param lookup_dictionary: dictionary where the key is a tuple of position and mutation, and value is
	tuple of region with variant, and observed and likelihood values
	:param excluded_sites: list of base positions to exclude from calculations
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param loci_dict: dictionary with detail annotations of loci per MITOMAP annotations
	:param kmer_dict: dictionary with all relevant annotations for every kmer
	:return: kmers dictionary with all relevant annotations for every kmer
	"""
	kmer_range = list(range(position, position + kmer_length))
	start = int(kmer_range[0])
	end = int(kmer_range[-1])
	if end > 16569:  # to handle circular genome
		kmer_range = list(range(start, 16570)) + list(range(1, end - 16569 + 1))
		end = int(kmer_range[-1])

	# initiate dictionary - can do this for each kmer since won't be added to later
	kmer_sum = calculate_oe.oe_functions.initialize_sum_dict(identifier_list=[(start, end)])
	
	# now sum the observed and likelihood values for the kmer, using the look-up dictionary
	for pos in kmer_range:
		if pos not in excluded_sites:
			for key in lookup_dictionary:
				if key[0] == str(pos):
					kmer_sum = calculate_oe.oe_functions.sum_obs_likelihood(
						mutation=key[1], identifier=(start, end), region=lookup_dictionary[key][0], dict=kmer_sum,
						observed=lookup_dictionary[key][1], likelihood=lookup_dictionary[key][2])
	
	# calculate observed and expected, and CI
	(total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI) = calculate_oe.oe_functions.calculate_oe(
		item=(start, end), sum_dict=kmer_sum, fit_parameters=fit_parameters, output="dict", file=None, max_parameter=2.5)

	# annotate with overlapping loci
	loci, loci_type = [], []  # list of overlapping loci
	for key in loci_dict:
		# condition for overlap, different for those spanning the artificial break
		if start < end:  # if start of kmer is less than end, ie not spanning artificial break
			if (start <= key[1]) and (end >= key[0]) and (loci_dict[key][0] not in loci):
				loci.append(loci_dict[key][0])
				loci_type.append(loci_dict[key][1])
		elif end < start:  # kmers overlapping the origin
			if ((end >= key[0]) and (end <= key[1])) or ((start >= key[0]) and (start <= key[1])) and (
					loci_dict[key][0] not in loci):
				loci.append(loci_dict[key][0])
				loci_type.append(loci_dict[key][1])
	loci = str(loci).strip('[]').replace("'", "")
	loci_type = str(loci_type).strip('[]').replace("'", "")
	
	# add entry to dictionary
	kmer_dict[(start, end)] = (total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI, loci, loci_type)
	return kmer_dict


def kmers(
	kmer_length: int,
	lookup_dictionary: Dict[Tuple[str, str], Tuple[str, float, float]],
	excluded_sites: List[int],
	fit_parameters: str,
	out_prefix: str):
	"""Function to calculate the obs:exp ratios and other annotations for kmers.

	:param kmer_length: the user-specified length of every kmer
	:param lookup_dictionary: dictionary where the key is a tuple of position and mutation, and value is
	tuple of region with variant, and observed and likelihood values
	:param excluded_sites: list of base positions to exclude from calculations
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param out_prefix: prefix to add to output file name
	:return: df, a pandas dataframe with details of all kmers and annotations
	:return: coding_intervals dictionary with detail annotations of loci per MITOMAP
	"""
	# generate dictionary of MITOMAP loci by their coordinates
	mitomap_loci = {}
	for row in csv.DictReader(open('required_files/databases/mitomap_genome_loci.txt'), delimiter='\t'):
		if int(row["Starting"]) > int(row["Ending"]):
			mitomap_loci[(1, int(row["Ending"]))] = (row["Map_Locus"], row["Description"])
			mitomap_loci[(int(row["Starting"]), 16569)] = (row["Map_Locus"], row["Description"])
		else:
			mitomap_loci[(int(row["Starting"]), int(row["Ending"]))] = (row["Map_Locus"], row["Description"])
	
	# initiate dictionaries, to be used in parallelized function
	manager = mp.Manager()
	kmers = manager.dict()
	pool = mp.Pool(mp.cpu_count())
	
	# calculate oe ratios and CIs for every possible kmer
	for position in list(range(1, 16570)):
		pool.apply_async(parallelize_kmers, args=(
			position, kmer_length, lookup_dictionary, excluded_sites, fit_parameters, mitomap_loci, kmers))
	pool.close()  # close pool and let all the processes complete
	pool.join()  # postpones the execution of next line of code until all processes in the queue are done
	
	# using pd to assign percentile ranks to each kmer by their upper CI value
	df = pd.DataFrame.from_dict(kmers, orient='index', columns=[
		'total_all', 'obs_max_het', 'exp_max_het', 'ratio_oe', 'lower_CI', 'upper_CI', 'loci', 'loci_type'])
	df.index = pd.MultiIndex.from_tuples(df.index, names=['start', 'end'])
	df = df.reset_index()
	print(df)
	df['pctRank_OEUF'] = df['upper_CI'].rank(pct=True, ascending=False)  # percentile rank by ratio obs:exp upper CI
	df = df.sort_values(['pctRank_OEUF'], ascending=False)  # smallest to largest
	
	# also rank separately for protein and RNA genes, and non-coding regions
	df['loci_type'] = df['loci_type'].map(str)
	df['pctRank_OEUF_RNA'] = df.upper_CI[
		df['loci_type'].str.contains('RNA')].rank(pct=True, ascending=False)  # if includes any RNA gene
	df['pctRank_OEUF_protein'] = df.upper_CI[
		df['loci_type'].str.contains('subunit|Cytochrome')].rank(pct=True, ascending=False)  # if includes any protein
	# only non-coding, excluding if contains RNA or protein sequence
	df['pctRank_OEUF_noncoding'] = df.upper_CI[
		~df['loci_type'].str.contains('RNA|subunit|Cytochrome')].rank(pct=True, ascending=False)

	# adding annotations, first create required dictionaries
	phylop_dict = build_model.annotate_mutations.phylop_annotate()
	uniprot_dict = build_model.annotate_mutations.uniprot_annotations()
	CI_dict = build_model.annotate_mutations.curated_func_sites()
	mod_dict = build_model.annotate_mutations.RNA_domains_mods()
	vep_dict = build_model.annotate_mutations.vep_annotate()
	prot_pos_dict = {}
	
	# create dictionary which converts DNA position to protein position
	for key in vep_dict:
		if key[3] == "codon":
			prot_pos_dict[int(key[1])] = str(vep_dict[key]).strip('[]').replace("'", "").replace(" ", "")
	# create dictionary of protein annotations
	sites_dict = {}
	for key in uniprot_dict:
		if "site" in str(uniprot_dict[key]):
			sites_dict[key] = [str(uniprot_dict[key]).strip('[]').replace("'", "").replace(" ", "")]
	# combine protein annotations into one dictionary
	for key in CI_dict:
		if key in sites_dict:
			sites_dict[key].append(CI_dict[key])
		else:
			sites_dict[key] = [CI_dict[key]]
	
	# initiate lists which will form basis of new column
	prot_start, prot_end, phylop_cons, prot_sites, rna_modified = [], [], [], [], []
	
	# iterate through each kmer
	for index in df.index:
		start = df['start'][index].astype(int)
		end = df['end'][index].astype(int)
		lookup_range = list(range(start, end + 1))
		if start > (16570 - kmer_length):  # to handle circular genome
			lookup_range = list(range(start, 16570)) + list(range(1, end + 1))
		
		# add protein position range - flip for MT-ND6, only protein on reverse strand
		if 'MT-ND6' in df['loci'][index]:
			prot_start.append(prot_pos_dict[end])
			prot_end.append(prot_pos_dict[start])
		else:
			prot_start.append(prot_pos_dict[start])
			prot_end.append(prot_pos_dict[end])
			
		# create new lists to handle multiple annotations per range
		p_cons, p_sites, r_mod = [], [], []
		# collect values across the kmer
		for i in lookup_range:
			if i != 3107:  # where ref=N
				# extract conservation metrics
				p_cons.append(float(phylop_dict[i]))
				# extract sites and modifications
				if (i in sites_dict) and (sites_dict[i] not in p_sites):
					p_sites.append(sites_dict[i])
				if ((str(i), "modified") in mod_dict) and (mod_dict[(str(i), "modified")] not in r_mod):
					r_mod.append(mod_dict[(str(i), "modified")])
				if (("\n" + str(i) + "\n") in open('required_files/other_annotations/rRNA_bridge_bases.txt').read()) and (
						"rRNA_bridge_base" not in r_mod):
					r_mod.append("rRNA_bridge_base")
		# take mean of conservation metric across kmer
		phylop_cons.append(mean(p_cons))
		# format sites and modifications
		prot_sites.append(str(p_sites).strip('[').strip(']').strip('\'').replace("', '", ",").replace(", ", "").strip())
		rna_modified.append(str(r_mod).strip('[').strip(']').strip('\'').replace("', '", ",").replace(", ", "").strip())
		
	# append annotations to df
	df['protein_position_start'] = prot_start
	df['protein_position_end'] = prot_end
	df['phylop_mean'] = phylop_cons
	df['protein_sites'] = prot_sites
	df['RNA_modifications'] = rna_modified
	print(df)

	# write to file
	np.savetxt(
		'output_files/local_constraint/%skmers_local_constraint.txt' % out_prefix, df, fmt='\t'.join(
			(['%i'] * 3) + (['%.6f'] * 5) + (['%s'] * 2) + (['%.6f'] * 4) + (['%s'] * 2) + (['%.6f'] * 1) + (['%s'] * 2)),
		header="start	end	variant_count	observed	expected	obs:exp	lower_CI	upper_CI	loci	loci_types	pctRank_OEUF	pctRank_OEUF_RNA	pctRank_OEUF_protein	pctRank_OEUF_noncoding	protein_pos_start	protein_pos_end	phylop_mean	protein_sites	RNA_sites")
	
	return df, mitomap_loci, sites_dict


def per_base(
	df: pd.DataFrame,
	kmer_length: int,
	mitomap_loci: Dict[Tuple[int, int], Tuple[str, str]]):
	"""Function to calculate the mean percentile ranking of overlapping kmers for every position in the mtDNA, and
	other annotations.

	:param df: a pandas dataframe with details of all kmers and associated annotations
	:param kmer_length: the user-specified length of every kmer
	:param mitomap_loci: dictionary with detail annotations of loci per MITOMAP annotations
	:return: df2, a pandas dataframe with details of statistics for every position in mtDNA and associated annotations
	"""
	
	dict = {}
	for i in list(range(1, 16570)):  # every base in the mtDNA
		# subset to all kmers overlapping position i
		if i > (16570 - kmer_length):  # to handle circular genome
			filtered_df = df[
				(
						(df['start'] >= np.int64(i - kmer_length)) &
						(df['end'] <= np.int64(i + kmer_length - 16570))
				) | (
						(df['start'] <= np.int64(i)) &
						(df['end'] >= np.int64(i))
				)]
		elif i <= kmer_length:  # to handle circular genome
			filtered_df = df[
				(
						(df['start'] >= np.int64(i - kmer_length + 16570)) &
						(df['end'] <= np.int64(i + kmer_length))
				) | (
						(df['start'] <= np.int64(i)) & (df['end'] >= np.int64(i))
				)]
		else:
			filtered_df = df[
				((df['start'] <= np.int64(i)) & (df['end'] >= np.int64(i))) |
				(df['start'] == np.int64(i)) |
				(df['end'] == np.int64(i))
			]
		# for each position, calculate the mean obs:exp ratio upper CI of overlapping kmers
		mean_OEUF = filtered_df['upper_CI'].mean()
		mean_exp = filtered_df['exp_max_het'].mean()

		# annotate regions to enable ranking separately for protein and RNA genes, and non-coding
		loci, loci_type = [], []
		for key in mitomap_loci:
			if ((i >= key[0]) and (i <= key[1])) or (i == key[0]) or (i == key[1]):
				loci.append(mitomap_loci[key][0])
				loci_type.append(mitomap_loci[key][1])
		dict[i] = (mean_OEUF, loci, loci_type, mean_exp)

	df2 = pd.DataFrame(dict).T.reset_index()
	df2.columns = ['pos', 'mean_OEUF', 'loci', 'loci_type', 'mean_exp']
	df2['pctRank_mean_OEUF'] = df2['mean_OEUF'].rank(pct=True, ascending=False)
	df2 = df2.sort_values(['pctRank_mean_OEUF'])
	
	# also rank RNA and protein separately, and non-coding
	df2['loci_type'] = df2['loci_type'].map(str)
	df2['pctRank_mean_OEUF_RNA'] = df2.mean_OEUF[
		df2['loci_type'].str.contains('RNA')].rank(pct=True, ascending=False)
	df2['pctRank_mean_OEUF_protein'] = df2.mean_OEUF[
		df2['loci_type'].str.contains('subunit|Cytochrome')].rank(pct=True, ascending=False)
	df2['pctRank_mean_OEUF_noncoding'] = df2.mean_OEUF[
		~df2['loci_type'].str.contains('RNA|subunit|Cytochrome')].rank(pct=True, ascending=False)
	
	return df2


def annotate(
	input_file: str,
	dataframe: pd.DataFrame,
	sites: Dict[int, List[str]],
	out_prefix: str):
	"""Function to provide per-base metrics for every position in the mtDNA.

	:param input_file: path to the file with this information, which is the output of annotate_output.py
	:param dataframe: a pandas dataframe with details of statistics for every position and annotations
	:param sites: dictionary where key is the position and value are protein site annotations
	:param out_prefix: prefix to add to output file name
	"""

	f = open('output_files/local_constraint/%sper_base_local_constraint.txt' % out_prefix, "w")
	header = "POS	REF	ALT	symbol	consequence	protein_position	MLC_pos_score	MLC_var_score	mean_OEUF	pctRank_mean_OEUF	pctRank_OEUF_protein	pctRank_OEUF_RNA	pctRank_OEUF_noncoding	mean_expected	in_phylotree	phyloP_score	protein_sites	RNA_modified	rRNA_bridge_base	apogee_class	mitotip_class	hmtvar_class	mitomap_status	mitomap_plasmy	mitomap_disease	clinvar_interp"
	f.write(header + '\n')
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mean_OEUF = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'mean_OEUF'].values).strip('[]')
		mean_exp = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'mean_exp'].values).strip('[]')
		pct_rank_OEUF = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF'].values).strip('[]')
		pct_rank_OEUF_RNA = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_RNA'].values).strip('[]')
		pct_rank_OEUF_protein = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_protein'].values).strip('[]')
		pct_rank_OEUF_noncoding = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_noncoding'].values).strip('[]')
		MLC_pos_score = pct_rank_OEUF
		
		# these scores are for local functional constraint, derived from RNA/noncoding and missense variants only
		# assign other functional variants within the range, determined from consideration of their oe distribution
		MLC_var_score = pct_rank_OEUF
		if ("stop_gain" in row["consequence"]) and ('stop_gained&start_lost' not in row["consequence"]):
			MLC_var_score = str(1)
			pct_rank_OEUF_protein = 'nan'
		if (("synonymous" in row["consequence"]) or ("stop_retained" in row["consequence"])) and not any(
				x in row["consequence"] for x in calculate_oe.oe_functions.more_severe_than_syn):
			MLC_var_score = str(0)
			pct_rank_OEUF_protein = 'nan'
		# these are start and stop loss - there are very few of these
		elif any(x in row["consequence"] for x in ["stop_lost", "start_lost", "incomplete_terminal"]):
			MLC_var_score = str(0.70)  # manually looked up equivalent
			pct_rank_OEUF_protein = 'nan'
			
		prot_sites = str(sites[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in sites else ''
		
		RNA_bridge = "Yes" if ("\n" + row["POS"] + "\n") in open(
			'required_files/other_annotations/rRNA_bridge_bases.txt').read() else "No"
		
		f.write(
			row["POS"] + '\t' + row["REF"] + '\t' + row["ALT"] + '\t' +
			row["symbol"] + '\t' + row["consequence"] + '\t' + row["protein_position"] + '\t' +
			str(MLC_pos_score) + '\t' + str(MLC_var_score) + '\t' +
			mean_OEUF + '\t' + pct_rank_OEUF + '\t' +
			pct_rank_OEUF_protein + '\t' + pct_rank_OEUF_RNA + '\t' + pct_rank_OEUF_noncoding + '\t' + mean_exp + '\t' +
			row["in_phylotree"] + '\t' + row["phyloP_score"] + '\t' +
			prot_sites + '\t' + row["RNA_modified"] + '\t' + RNA_bridge + '\t' +
			row["apogee_class"] + '\t' + row["mitotip_class"] + '\t' + row["hmtvar_class"] + '\t' +
			row["mitomap_status"] + '\t' + row["mitomap_plasmy"] + '\t' + row["mitomap_disease"] + '\t' +
			row["clinvar_interp"] + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	parser.add_argument(
		"-obs", type=str, help="Population dataset from which observed maximum heteroplasmy is obtained")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-prefix", type=str, help="Output files prefix")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calculations")
	parser.add_argument(
		"-kmer_length", type=int, help="Window size to use for analysis in bp")
	args = parser.parse_args()
	
	# set default, for gnomAD
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "gnomad_max_hl"
	if args.parameters is None:
		args.parameters = 'output_files/calibration/linear_model_fits.txt'
	if args.prefix is None:
		args.prefix = ""
	if args.exc_sites is None:
		# exclude “artifact_prone_sites” in gnomAD positions - 301, 302, 310, 316, 3107, and 16182 (3107 already excluded)
		# variants at these sites were not called in gnomAD and therefore these positions removed from calculations
		args.exc_sites = [301, 302, 310, 316, 16182]
	if args.kmer_length is None:
		args.kmer_length = 30
	
	for path in ['output_files/local_constraint']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Analyze local constraint in the mtDNA!")
	print("This may take a little time depending on available computing power...")
	
	(df, mitomap_loci, sites_dict) = kmers(
		kmer_length=args.kmer_length, lookup_dictionary=lookup_dict(input_file=args.input, obs_value=args.obs),
		excluded_sites=args.exc_sites, fit_parameters=args.parameters, out_prefix=args.prefix)
	
	df2 = per_base(df=df, kmer_length=args.kmer_length, mitomap_loci=mitomap_loci)
	
	annotate(input_file=args.input, dataframe=df2, sites=sites_dict, out_prefix=args.prefix)

	print(datetime.datetime.now(), "Script finished!")
