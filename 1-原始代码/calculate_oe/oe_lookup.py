import argparse
import datetime
from oe_functions import *
import os


def oe_lookup(
		start: int, end: int,
		input_file: str, obs_value: str, fit_parameters: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval between a pair of given coordinates.

	:param start: start coordinate
	:param end: end coordinate
	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	# initialize dictionary so all values are 0
	loci_sum = initialize_sum_dict(identifier_list=['all SNVs'])
	# for each locus, sum observed and likelihood for each variant class
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			# this includes ALL SNVs, and start has to be less than end
			if (int(row["POS"]) >= start) and (int(row["POS"]) <= end):
				loci_sum = sum_obs_likelihood(
					mutation=mutation, identifier='all SNVs', region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=loci_sum)
	# calculate ratio and confidence interval
	for variant_type in ['all SNVs']:
		exp_max_het = calculate_exp(sum_dict=loci_sum, identifier=variant_type, fit_parameters=fit_parameters)
		obs_max_het = calculate_obs(identifier=variant_type, sum_dict=loci_sum)
		ratio_oe = obs_max_het / exp_max_het
		total_all = calculate_total(identifier=variant_type, sum_dict=loci_sum)
		(lower_CI, upper_CI) = calculate_CI(
			obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het)
	print(
		"For interval", start, "to", end, "the observed and expected values are", obs_max_het, "and", exp_max_het,
		". The ratio is ", ratio_oe, "and the OEUF is", upper_CI, ".")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	parser.add_argument(
		"-obs", type=str, help="Population dataset from which observed maximum heteroplasmy is obtained")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calibration")
	parser.add_argument(
		"-start", type=int, help="Start coordinate for the interval to calculate within")
	parser.add_argument(
		"-end", type=int, help="Start coordinate for the interval to calculate within")
	args = parser.parse_args()
	
	# set defaults, for gnomAD
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
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for a given interval")
	
	oe_lookup(
		start=args.start, end=args.end, input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
