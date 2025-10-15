import argparse
import datetime
from oe_functions import *
import os
from typing import List


def consequences_oe(
		input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for functional classes of variation.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%sconsequences_obs_exp.txt' % output_prefix, "w")
	header = "consequence	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	consequences = [
			"synonymous", "missense", "stop_gain", "start_lost", "stop_lost", "rRNA", "tRNA", "intergenic"]
	
	conseq_sum = initialize_sum_dict(identifier_list=consequences)
	# for each locus, sum observed and likelihood for each variant type
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			name = ''
			# to get as most severe
			if ("synonymous" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_syn):
				name = 'synonymous'
			elif ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
				name = 'missense'
			elif ("stop_gain" in row["consequence"]) and ('stop_gained&start_lost' not in row["consequence"]):
				name = 'stop_gain'
			elif "start_lost" in row["consequence"]:
				name = 'start_lost'
			# note incomplete_terminal_codon_variant are equivalent to stop lost
			elif ("stop_lost" in row["consequence"]) or ("incomplete_terminal" in row["consequence"]):
				name = 'stop_lost'
			elif row["symbol"].startswith('MT-R'):
				name = 'rRNA'
			elif row["symbol"].startswith('MT-T'):
				name = 'tRNA'
			elif "intergenic" in row["consequence"]:
				name = 'intergenic'
			if name != '':
				conseq_sum = sum_obs_likelihood(
					mutation=mutation, identifier=name, region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=conseq_sum)
	
	for consequence in consequences:
		calculate_oe(item=consequence, sum_dict=conseq_sum, fit_parameters=fit_parameters, file=file)


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
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calibration")
	args = parser.parse_args()
	
	# set defaults, for gnomAD
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
	
	for path in ['output_files/oe']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for each functional class of variation")
	
	consequences_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
	