import argparse
import datetime
from oe_functions import *
import os
from typing import List


def calibrate(input_file: str, obs_value: str, output_prefix: str, excluded_sites: List[int]):
	"""Sum the observed maximum heteroplasmy, mutation likelihood scores and count for neutral variants,
	in the reference sequence excluding the OriB-OriH region.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	# first, extract list of genes/loci and their length
	gene_length = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		for gene in row["symbol"].split(','):  # split is to handle positions in two genes
			# mark control region separately to other non-coding
			gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
			gene = 'other_non-coding' if gene == '' else gene
			# don't include positions in the ori region or in excluded sites
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in excluded_sites):
				if gene not in gene_length:
					gene_length[gene] = 1
				else:
					gene_length[gene] += 1
	
	# now initialize a dictionary that will be used to sum values for each gene/loci, and set all values to 0
	calibration = initialize_sum_dict(identifier_list=list(gene_length.keys()))
	
	# determine the phyloP threshold that will be used to identify neutral variants - this is the bottom decile
	phylop = []
	catch_list = []  # use to get to unique positions, since this file has three rows per positions (ie 3 possible SNVs)
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # second element [1] is bottom decile

	# write a file of all the neutral variants used
	f = open('output_files/calibration/neutral_variants_used.txt', "w")
	f.write('variant' + '\t' + 'position' + '\t' + 'source' + '\n')
	
	# now build dictionary
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for gene in row["symbol"].split(','):  # to handle variants within two genes
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in excluded_sites):
				# criteria for neutral = haplogroup variant in phylotree or in bottom decile phyloP (threshold)
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# mark control region separately to other non-coding
					gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
					gene = 'other_non-coding' if gene == '' else gene
					# sum values for each gene/loci
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=gene, region='ref_exc_ori', dict=calibration,
						observed=row[obs_value], likelihood=row["Likelihood"])
					# write to file, note variants in two genes are written twice
					reason = 'haplogroup_variant' if (int(row["in_phylotree"]) == 1) else 'lowest_decile_phyloP'
					f.write('m.' + str(row["POS"]) + row["REF"] + '>' + row["ALT"] + '\t' + str(row["POS"]) + '\t' + reason + '\n')
	
	# write to file for plotting
	f = open('output_files/calibration/%sloci_obs_vs_scores.txt' % output_prefix, "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to calibrate
		for gene in gene_length:
			f.write(
				mut_group + '\t' + gene + '\t' +
				str(calibration[(gene, mut_group, 'obs_max_het', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'sum_LR', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'count', 'ref_exc_ori')]) + '\t' +
				str(gene_length[gene] / 3) + '\n')  # need to divide by 3 as count is number positions x 3 (possible SNVs)


def calibrate_ori(input_file: str, obs_value: str, output_prefix: str, excluded_sites: List[int]):
	"""Sum the observed maximum heteroplasmy, mutation likelihood scores and count for neutral variants,
	in the OriB-OriH region.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	# first, segment the ori region into equally sized blocks of 70 bp, using approximate size of tRNA genes
	ori_blocks = {}
	block = 1
	for i in range(0, len(ori_region), 70):
		if len(ori_region[i:i + 70]) == 70:
			ori_blocks["block_" + str(block)] = ori_region[i:i + 70]
		else:  # to ensure all ori bases included, make the last block larger if needed
			ori_blocks["block_" + str(block - 1)] = ori_blocks["block_" + str(block - 1)] + \
													ori_region[i:i + ori_region[-1]]
		block += 1
	
	# now initialize a dictionary that will be used to sum values for each block, and set all values to 0
	calibration = initialize_sum_dict(identifier_list=list(ori_blocks.keys()))
	
	# determine the phyloP threshold that will be used to identify neutral variants - this is the bottom decile
	phylop = []
	catch_list = []  # use to get to unique positions, since this file has three rows per positions (ie 3 possible SNVs)
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # second element [1] is bottom decile
	
	# write a file of all the neutral variants used
	f = open('output_files/calibration/neutral_variants_used.txt', "a")
	
	# now build dictionary
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for n_block in ori_blocks:
			# if the variant is in the ori block in the loop
			if (int(row["POS"]) in ori_blocks[n_block]) and (int(row["POS"]) not in excluded_sites):
				# criteria for neutral = haplogroup variant in phylotree or in bottom decile phyloP (threshold)
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# sum values for each block
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=n_block, region='ori', dict=calibration,
						observed=row[obs_value], likelihood=row["Likelihood"])
					# write to file
					reason = 'haplogroup_variant' if (int(row["in_phylotree"]) == 1) else 'lowest_decile_phyloP'
					f.write('m.' + str(row["POS"]) + row["REF"] + '>' + row["ALT"] + '\t' + str(row["POS"]) + '\t' + reason + '\n')
	
	# write to file for plotting
	f = open('output_files/calibration/%sloci_obs_vs_scores_ori.txt' % output_prefix, "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to calibrate
		for n_block in ori_blocks:
			f.write(
				mut_group + '\t' + n_block + '\t' +
				str(calibration[(n_block, mut_group, 'obs_max_het', 'ori')]) + '\t' +
				str(calibration[(n_block, mut_group, 'sum_LR', 'ori')]) + '\t' +
				str(calibration[(n_block, mut_group, 'count', 'ori')]) + '\t' +
				str(len(ori_blocks[n_block])) + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	parser.add_argument(
		"-obs", type=str, help="Population dataset from which observed maximum heteroplasmy is obtained")
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
	if args.prefix is None:
		args.prefix = ""
	if args.exc_sites is None:
		# exclude “artifact_prone_sites” in gnomAD positions - 301, 302, 310, 316, 3107, and 16182 (3107 already excluded)
		# variants at these sites were not called in gnomAD and therefore these positions removed from calculations
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	for path in ['output_files/calibration']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Preparing files for model calibration!" + '\n')
	
	calibrate(input_file=args.input, obs_value=args.obs, output_prefix=args.prefix, excluded_sites=args.exc_sites)
	calibrate_ori(input_file=args.input, obs_value=args.obs, output_prefix=args.prefix, excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
