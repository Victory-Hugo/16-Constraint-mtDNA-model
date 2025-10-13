import csv
import numpy as np
from scipy.stats import beta
from typing import Dict, List, TextIO, Tuple, Union


def initialize_sum_dict(identifier_list: List[Union[str, Tuple[int, int]]]):
	"""Generate a dictionary that will be used to sum the observed maximum heteroplasmy, mutation likelihood scores and
	count for a group of variants.

	:param identifier_list: identifiers for each of the variant categories for which the values are being summed
	:return: a dictionary with a tuple key of the identifier, mutation group, value to sum, and region with variant,
	and value of 0 for each key
	"""
	# initialize dictionary so all values are 0
	dict = {}
	for identifier in identifier_list:
		for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to fit
			for value in ['obs_max_het', 'sum_LR', 'count']:  # the three values to sum
				for region in ['ref_exc_ori', 'ori']:  # the two mutational models that could be applied
					dict[identifier, mut_group, value, region] = 0
	return dict


def sum_obs_likelihood(
		mutation: str, identifier: Union[str, Tuple[int, int]], region: str,
		observed: Union[str, float], likelihood: Union[str, float],
		dict: Dict[Tuple[str, str, str, str], Union[float, int]]):
	"""Sum the observed maximum heteroplasmy, mutation likelihood scores and count for a category of variants.

	:param mutation: the mutation type, in Ref>Alt format
	:param identifier: the identifier of the variant category for which the values are being summed
	:param region: should be either 'ref_exc_ori' or 'ori', to indicate which mutational model to apply
	:param observed: the value to use for the observed maximum heteroplasmy
	:param likelihood: the value to use for the mutation likelihood score
	:param dict: the name of the dictionary to append to
	:return: a dictionary with a tuple key of the identifier, mutation group, value to sum, and region with variant,
	and value corresponding to the value to sum in key
	"""
	# highly mutable G>A and T>C are fit separately
	mut_group = 'G>A_and_T>C' if (mutation == 'G>A' or mutation == 'T>C') else 'other'
	# note the dictionary should already be initialized (to start at 0 for all anticipated keys)
	dict[(identifier, mut_group, 'obs_max_het', region)] += float(observed)
	dict[(identifier, mut_group, 'sum_LR', region)] += float(likelihood)
	dict[(identifier, mut_group, 'count', region)] += 1
	return dict


def linear_model(count: int, intercept: float, coefficient: float, sum_lr: float):
	"""Apply a linear model fit on neutral variation to calculate the expected sum maximum heteroplasmy.
	A clamp is applied such that the expected value cannot be <0, or greater than the maximum possible value.

	:param count: the number of variants to calculate the expected for
	:param intercept: linear equation intercept
	:param coefficient: linear equation coefficient
	:param sum_lr: the sum of the mutation likelihood scores for the variants to calculate the expected for
	:return: expected sum maximum heteroplasmy
	"""
	return max(float(0), min(float(count), intercept + (coefficient * sum_lr)))


def calibrate(
		sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]], mut_group: str, identifier: str, region: str,
		parameters: Dict[Tuple[str, str, str], float]):
	"""Apply the linear model, which returns the expected sum maximum heteroplasmy.

	:param sum_dict: dictionary with the possible variant count and sum mutation scores needed for calculation
	:param mut_group: mutation group to calibrate, either 'G>A_and_T>C' or 'other'
	:param identifier: the identifier of the variant category to calibrate
	:param region: should be either 'ref_exc_ori' or 'ori', to indicate which calibration to apply
	:param parameters: dictionary with the linear equation coefficients and intercepts to use for calibration
	:return: the output of linear model, the expected sum maximum heteroplasmy for the mutation group and region
	"""
	return linear_model(
		count=sum_dict[(identifier, mut_group, 'count', region)],
		intercept=parameters[(mut_group, region, 'intercept')],
		coefficient=parameters[(mut_group, region, 'coefficient')],
		sum_lr=sum_dict[(identifier, mut_group, 'sum_LR', region)])


def calculate_exp(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]], fit_parameters: str):
	"""Calculate the expected sum maximum heteroplasmy for a category of variants.

	:param sum_dict: dictionary with the possible variant count and sum mutation scores needed for calculation
	:param identifier: the identifier of the variant category to calculate the expected for
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:return: total expected sum maximum heteroplasmy
	"""
	# create dictionary with the linear equation coefficients and intercepts to use for calculating expected
	parameters = {}
	for row in csv.DictReader(open(fit_parameters), delimiter='\t'):
		parameters[(row["mutation_group"], row["region"], row["item"])] = float(row["value"])
	
	# calculate expected separately for 'G>A_and_T>C' and 'other' mutation groups, and 'ref_exc_ori' and 'ori' regions
	# if the group has 0 variants, the calibrate function will return 0 due to the clamp
	exp = calibrate(
			sum_dict=sum_dict, mut_group='G>A_and_T>C', identifier=identifier, region='ref_exc_ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='other', identifier=identifier, region='ref_exc_ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='G>A_and_T>C', identifier=identifier, region='ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='other', identifier=identifier, region='ori', parameters=parameters)
	return exp


def calculate_obs(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]]):
	"""Calculate the observed sum maximum heteroplasmy for a category of variants.

	:param identifier: the identifier of the variant category to calculate the observed for
	:param sum_dict: dictionary with the observed maximum heteroplasmy
	:return: total observed sum maximum heteroplasmy, across both mutation groups and region types
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'obs_max_het', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'obs_max_het', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'obs_max_het', 'ori')] + sum_dict[
		(identifier, 'other', 'obs_max_het', 'ori')]


def calculate_total(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]]):
	"""Calculate the total count of possible variants for a group of variants.

	:param identifier: the identifier of the variant category to calculate the observed for
	:param sum_dict: dictionary with the observed maximum heteroplasmy
	:return: total count of possible variants, across both mutation groups and region types
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'count', 'ori')] + sum_dict[
		(identifier, 'other', 'count', 'ori')]


def calculate_CI(obs_max_het: float, total: int, exp_max_het: float, max_parameter: float = 2.5):
	"""Calculate the 90% confidence interval around the observed:expected ratio using a beta distribution.

	:param obs_max_het: observed sum maximum heteroplasmy
	:param total: number of possible variants
	:param exp_max_het: expected sum maximum heteroplasmy
	:param max_parameter: maximum varying parameter value, default = 2.5
	:return: lower_CI and upper_CI, the lower and upper bounds of the 90% confidence interval
	"""
	# compute the density of the beta distribution for a given pair of observed and expected values
	# use expected times a varying parameter, define range of varying parameters to assess
	varying_parameters = np.arange(0.001, max_parameter, 0.001).tolist()
	
	# express observed as proportion of possible to bound between 0-1
	prop_poss_maxhet = obs_max_het / total
	# if 0 or 1 set within range to enable assessment, use minimum possible observed to be conservative
	if prop_poss_maxhet == 0:
		prop_poss_maxhet = 0.1 / total
	elif prop_poss_maxhet == 1:
		prop_poss_maxhet = 1 - (0.1 / total)
	
	beta_distrib = []
	for vp in varying_parameters:
		a = (exp_max_het * vp) + 1  # number 'successes'
		b = total - (exp_max_het * vp) + 1  # number 'failures'
		beta_distrib.append(beta.pdf(prop_poss_maxhet, a, b))

	# compute cdf of this function
	cumulative_beta = []
	for value in np.cumsum(beta_distrib):
		cumulative_beta.append(value / max(np.cumsum(beta_distrib)))

	# extract the value of the varying parameter at points corresponding to 5% and 95%, for bounds of 90% CI
	lower_list, upper_list = [], []
	for value in cumulative_beta:
		if value < 0.05:
			lower_list.append(cumulative_beta.index(value))
		if value > 0.95:
			upper_list.append(cumulative_beta.index(value))
	
	if (obs_max_het > 0) and (len(lower_list) > 0):
		lower_idx = max(lower_list)  # what is the position of the maximum value that has probability <0.05
		lower_CI = varying_parameters[lower_idx]  # 5th
	else:
		lower_CI = 0
	upper_idx = min(upper_list)  # what is the position of the minimum value that has probability >0.95
	upper_CI = varying_parameters[upper_idx]  # 95th
	
	return lower_CI, upper_CI


def calculate_pvalue(obs_max_het: float, total: int, exp_max_het: float):
	"""Calculate if the observed value is significantly lower than the expected, using a beta distribution.

	:param obs_max_het: observed sum maximum heteroplasmy
	:param total: number of possible variants
	:param exp_max_het: expected sum maximum heteroplasmy
	:return: less_than_pvalue, the probability of observing a value that is less than or equal to the observed,
	given the expected value
	"""
	# express observed as proportion of possible to bound between 0-1
	prop_poss_maxhet = obs_max_het / total
	# if 0 or 1 set within range to enable assessment, use minimum possible observed to be conservative
	if prop_poss_maxhet == 0:
		prop_poss_maxhet = 0.1 / total
	elif prop_poss_maxhet == 1:
		prop_poss_maxhet = 1 - (0.1 / total)
	a = exp_max_het + 1
	b = total - exp_max_het + 1
	less_than_pvalue = beta.cdf(prop_poss_maxhet, a, b)
	return less_than_pvalue


def calculate_oe(
		item: Union[str, Tuple[int, int]], sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		fit_parameters: str, file: Union[TextIO, None], output: str = "write", max_parameter: float = 2.5):
	"""Calculate observed, expected, obs:exp ratio, and 90% confidence interval.

	:param item: the identifier of the variant category to calculate the expected for
	:param sum_dict: dictionary with the possible variant count and sum mutation scores needed for calculation
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param file: the file to write the results to, if output is "write"
	:param output: whether to write to file, or save as dictionary, default is "write", other "dict"
	:param max_parameter: maximum varying parameter value for calculate_CI, default = 2.5
	"""
	exp_max_het = calculate_exp(sum_dict=sum_dict, identifier=item, fit_parameters=fit_parameters)
	obs_max_het = calculate_obs(sum_dict=sum_dict, identifier=item)
	ratio_oe = obs_max_het / exp_max_het
	total_all = calculate_total(sum_dict=sum_dict, identifier=item)

	(lower_CI, upper_CI) = calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=max_parameter)
	
	print("Calculating values for", item)
	if output == "write":
		file.write(
			item + '\t' + str(total_all) + '\t' + str(obs_max_het) + '\t' + str(exp_max_het)
			+ '\t' + str(ratio_oe) + '\t' + str(lower_CI) + '\t' + str(upper_CI) + '\n')
	elif output == "dict":
		return total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI


# variables to use

# ori refers to OriB-OriH region with known difference in mutational signature, m.16197-191, across artificial break
ori_region = list(range(16197, 16570)) + list(range(1, 191 + 1))

# Possible consequences in protein genes, ranked in order of severity:
# stop_gained
# stop_gained&start_lost, start_lost, stop_lost, incomplete_terminal_codon_variant&coding_sequence_variant (=start lost)
# missense_variant
# synonymous_variant, stop_retained_variant (equivalent to synonymous, very few)
more_severe_than_missense = ["stop_gain", "start_lost", "stop_lost", "incomplete_terminal"]
more_severe_than_syn = ["missense", "stop_gain", "start_lost", "stop_lost", "incomplete_terminal"]
