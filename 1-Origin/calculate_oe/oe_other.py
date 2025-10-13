import argparse
import datetime
from oe_functions import *
import os


def protein_annot_oe(
		input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for missense by UniProt and other annotations.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%sprotein_annotation_obs_exp.txt' % output_prefix, "w")
	header = "protein_annotation	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# assess functional sites, and transmembrane which is the only topological domain available for all proteins
	prot_sum = initialize_sum_dict(identifier_list=["site", "transmembrane", "not-transmembrane"])
	
	# for each annotation, sum observed and likelihood scores for missense variants
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			# only include missense where it is the most severe consequence
			if ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
				if ('site' in row["uniprot_annotation"]) or row["other_prot_annotation"]:
					prot_sum = sum_obs_likelihood(
						mutation=mutation, identifier="site", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"], dict=prot_sum)
						
	for annot in ["site"]:
		calculate_oe(item=annot, sum_dict=prot_sum, fit_parameters=fit_parameters, file=file)


def tRNA_pos_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval at each tRNA position.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%stRNA_position_obs_exp.txt' % output_prefix, "w")
	header = "tRNA_position	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# create list of tRNA positions
	tRNA_positions = []
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if (row["tRNA_position"] not in tRNA_positions) and (',' not in row["tRNA_position"]):
			tRNA_positions.append(row["tRNA_position"])
	tpos_sum = initialize_sum_dict(identifier_list=tRNA_positions)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			if row["symbol"].startswith('MT-T'):
				for tpos in row["tRNA_position"].split(','):  # some positions have two since in two genes
					tpos_sum = sum_obs_likelihood(
						mutation=mutation, identifier=tpos, region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"], dict=tpos_sum)
	# so this is currently double counting tRNA positions in two genes
	
	for tRNA_pos in tRNA_positions:
		calculate_oe(item=tRNA_pos, sum_dict=tpos_sum, fit_parameters=fit_parameters, file=file)


def RNA_base_types_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for different base types in RNA genes.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%sRNA_base_types_obs_exp.txt' % output_prefix, "w")
	header = "RNA_base_type	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# list of RNA base types
	RNA_types = ["WC_tRNA", "non-WC_tRNA", "loop-or-other_tRNA", "WC_rRNA", "non-WC_rRNA", "loop-or-other_rRNA"]
	RNA_type_sum = initialize_sum_dict(identifier_list=RNA_types)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			for RNA_type in row["RNA_base_type"].split(','):  # to iterate through
				if row["symbol"].startswith('MT-T'):
					RNA_type_sum = sum_obs_likelihood(
						mutation=mutation, identifier=RNA_type + "_tRNA", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"], dict=RNA_type_sum)
				if row["symbol"].startswith('MT-R'):
					RNA_type_sum = sum_obs_likelihood(
						mutation=mutation, identifier=RNA_type + "_rRNA", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"], dict=RNA_type_sum)
	
	for RNA_type in RNA_types:
		calculate_oe(item=RNA_type, sum_dict=RNA_type_sum, fit_parameters=fit_parameters, file=file)
	
	# also capture RNA modified bases
	RNA_mods = ["modified_tRNA", "modified_rRNA", "non-modified_tRNA", "non-modified_rRNA"]
	RNA_mod_sum = initialize_sum_dict(identifier_list=RNA_mods)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		identifier = ""
		if int(row["POS"]) not in excluded_sites:
			if row["symbol"].startswith('MT-T'):
				identifier = "modified_tRNA" if row["RNA_modified"] else "non-modified_tRNA"
			if row["symbol"].startswith('MT-R'):
				identifier = "modified_rRNA" if row["RNA_modified"] else "non-modified_rRNA"
			if identifier != "":
				RNA_mod_sum = sum_obs_likelihood(
					mutation=mutation, identifier=identifier, region='ref_exc_ori',
					observed=row[obs_value], likelihood=row["Likelihood"], dict=RNA_mod_sum)
	
	for mod in RNA_mods:
		calculate_oe(item=mod, sum_dict=RNA_mod_sum, fit_parameters=fit_parameters, file=file)

	
def tRNA_domain_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for different tRNA domains.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%stRNA_domains_obs_exp.txt' % output_prefix, "w")
	header = "tRNA_domain	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# create list of tRNA domains
	tRNA_domains = []
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if (row["tRNA_domain"] not in tRNA_domains) and (',' not in row["tRNA_domain"]):
			tRNA_domains.append(row["tRNA_domain"])
	tRNA_dom_sum = initialize_sum_dict(identifier_list=tRNA_domains)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			if row["symbol"].startswith('MT-T'):
				for tRNA_dom in row["tRNA_domain"].split(','):  # to iterate through
					tRNA_dom_sum = sum_obs_likelihood(
						mutation=mutation, identifier=tRNA_dom, region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"], dict=tRNA_dom_sum)
	
	for tRNA_dom in tRNA_domains:
		calculate_oe(item=tRNA_dom, sum_dict=tRNA_dom_sum, fit_parameters=fit_parameters, file=file)
	

def insilico_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for different in silico predictors.
	APOGEE is for missense, MitoTip and HmtVar are for tRNA. All three are recommended in ACMG modifications.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%sinsilicos_obs_exp.txt' % output_prefix, "w")
	header = "prediction	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# start with APOGEE, create list of classifications
	apogee = ["Pathogenic-APOGEE", "Neutral-APOGEE"]
	apogee_sum = initialize_sum_dict(identifier_list=apogee)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			# now only include missense where it is the most severe consequence
			if ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
				for apogee_class in row["apogee_class"].split(','):  # some positions have two since in two genes
					if apogee_class:
						apogee_sum = sum_obs_likelihood(
							mutation=mutation, identifier=apogee_class + "-APOGEE", region='ref_exc_ori',
							observed=row[obs_value], likelihood=row["Likelihood"], dict=apogee_sum)
	
	for apogee_class in apogee:
		calculate_oe(item=apogee_class, sum_dict=apogee_sum, fit_parameters=fit_parameters, file=file)
		
	# now MitoTip
	mitotip = [
		'likely pathogenic-MitoTip', 'possibly pathogenic-MitoTip', 'possibly benign-MitoTip', 'likely benign-MitoTip']
	mitotip_sum = initialize_sum_dict(identifier_list=mitotip)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			if row["symbol"].startswith('MT-T'):
				mitotip_sum = sum_obs_likelihood(
					mutation=mutation, identifier=row["mitotip_class"] + "-MitoTip", region='ref_exc_ori',
					observed=row[obs_value], likelihood=row["Likelihood"], dict=mitotip_sum)
	
	for mitotip_class in mitotip:
		calculate_oe(item=mitotip_class, sum_dict=mitotip_sum, fit_parameters=fit_parameters, file=file)
	
	# now HmtVar
	hmtvar = ['pathogenic-HmtVar', 'likely_pathogenic-HmtVar', 'likely_polymorphic-HmtVar', 'polymorphic-HmtVar']
	hmtvar_sum = initialize_sum_dict(identifier_list=hmtvar)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			if row["symbol"].startswith('MT-T') and row["hmtvar_class"] != "None":
				hmtvar_sum = sum_obs_likelihood(
					mutation=mutation, identifier=row["hmtvar_class"] + "-HmtVar", region='ref_exc_ori',
					observed=row[obs_value], likelihood=row["Likelihood"], dict=hmtvar_sum)
	
	for hmtvar_class in hmtvar:
		calculate_oe(item=hmtvar_class, sum_dict=hmtvar_sum, fit_parameters=fit_parameters, file=file)
			

def disease_vars_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for different disease-associated variants.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%sdisease_variants_obs_exp.txt' % output_prefix, "w")
	header = "classification	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# first for MITOMAP
	mitomap_status = ["Cfrm-MITOMAP", "Reported-MITOMAP"]
	mitomap_sum = initialize_sum_dict(identifier_list=mitomap_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			status = ''
			if "Cfrm" in row["mitomap_status"]:
				status = "Cfrm-MITOMAP"
			elif "Reported" in row["mitomap_status"]:
				status = "Reported-MITOMAP"
			if status != '':
				mitomap_sum = sum_obs_likelihood(
					mutation=mutation, identifier=status, region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=mitomap_sum)
	
	for status in mitomap_status:
		calculate_oe(item=status, sum_dict=mitomap_sum, fit_parameters=fit_parameters, file=file)
	
	# next clinvar
	clinvar_status = [
		"Benign-ClinVar", "Likely benign-ClinVar", "Uncertain significance-ClinVar", "Likely pathogenic-ClinVar", "Pathogenic-ClinVar"]
	clinvar_sum = initialize_sum_dict(identifier_list=clinvar_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			if (row["clinvar_interp"] + "-ClinVar") in clinvar_status:
				region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
				clinvar_sum = sum_obs_likelihood(
					mutation=mutation, identifier=row["clinvar_interp"] + "-ClinVar", region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=clinvar_sum)
	
	for status in clinvar_status:
		calculate_oe(item=status, sum_dict=clinvar_sum, fit_parameters=fit_parameters, file=file)


def vus_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""Calculate the observed:expected ratio and 90% confidence interval for different subsets of VUS.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	:param obs_value: the column header of observed value
	:param fit_parameters: the path to the file with the linear equation coefficients and intercepts to use
	:param output_prefix: string added to start of output file name
	:param excluded_sites: list of base positions to exclude from calculations
	"""
	file = open('output_files/oe/%svus_obs_exp.txt' % output_prefix, "w")
	header = "classification	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# first for MITOMAP
	mitomap_status = ["MITOMAP-Reported", "MITOMAP-Reported-PP3-PM2s", "MITOMAP-Reported-BP4-BS1"]
	mitomap_sum = initialize_sum_dict(identifier_list=mitomap_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			status = ''
			if "Reported" in row["mitomap_status"]:
				status = "MITOMAP-Reported"
				if row["consequence"] == "missense_variant":  # just missense
					if (row["apogee_class"] == "Pathogenic") and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif (row["apogee_class"] == "Neutral") and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
				if row["symbol"].startswith('MT-T'):  # subset to tRNA
					if ("pathogenic" in row["mitotip_class"]) and ("pathogenic" in row["hmtvar_class"])\
							and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif ("benign" in row["mitotip_class"]) and ("polymorphic" in row["hmtvar_class"])\
							and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
			if status != '':
				mitomap_sum = sum_obs_likelihood(
					mutation=mutation, identifier=status, region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=mitomap_sum)
	
	for status in mitomap_status:
		calculate_oe(item=status, sum_dict=mitomap_sum, fit_parameters=fit_parameters, file=file)
	
	# next clinvar
	clinvar_status = ["ClinVar-Uncertain significance", "ClinVar-Uncertain significance-PP3-PM2s"]
	# manually removing "Uncertain significance-ClinVar-BP4-BS1" as none satisfy this, quick fix to avoid error
	clinvar_sum = initialize_sum_dict(identifier_list=clinvar_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			status = ''
			if row["clinvar_interp"] == "Uncertain significance":
				status = "ClinVar-Uncertain significance"
				if row["consequence"] == "missense_variant":  # just missense
					if (row["apogee_class"] == "Pathogenic") and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif (row["apogee_class"] == "Neutral") and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
				if row["symbol"].startswith('MT-T'):  # subset to tRNA
					if ("pathogenic" in row["mitotip_class"]) and ("pathogenic" in row["hmtvar_class"]) and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif ("benign" in row["mitotip_class"]) and ("polymorphic" in row["hmtvar_class"]) and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
			if status != '':
				clinvar_sum = sum_obs_likelihood(
					mutation=mutation, identifier=status, region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=clinvar_sum)
	
	for status in clinvar_status:
		calculate_oe(item=status, sum_dict=clinvar_sum, fit_parameters=fit_parameters, file=file)


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
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for protein annotations")

	protein_annot_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)

	print(datetime.datetime.now(), "Calculate the observed:expected ratio for tRNA position")

	tRNA_pos_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)

	print(datetime.datetime.now(), "Calculate the observed:expected ratio for RNA base types")

	RNA_base_types_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)

	print(datetime.datetime.now(), "Calculate the observed:expected ratio for tRNA domains")

	tRNA_domain_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)

	print(datetime.datetime.now(), "Calculate the observed:expected ratio for in silicos")

	insilico_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)

	print(datetime.datetime.now(), "Calculate the observed:expected ratio for disease variants")

	disease_vars_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for subsets of VUS")
	
	vus_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
