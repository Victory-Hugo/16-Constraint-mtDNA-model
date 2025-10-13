import argparse
import csv
import datetime
import os
from statistics import median


def gene_conservation(input_file: str):
	"""Calculate gene conservation metrics, to compare to observed:expected ratio values.

	:param input_file: annotated file with mutation scores and observed maximum heteroplasmy
	"""
	# first create dictionary of the phylop conservation scores for each mtDNA base
	dict = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		dict[row["POS"]] = float(row["phyloP_score"])
	
	# to handle bases in two genes, create a list of base coordinates in each gene
	gene_dict = {}
	for row in csv.DictReader(open(
			'required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		if row["SYMBOL"] not in gene_dict:
			gene_dict[row["SYMBOL"]] = []
		elif row["POS"] not in gene_dict[row["SYMBOL"]]:  # since each position will be listed multiple times in this file
			gene_dict[row["SYMBOL"]].append(row["POS"])
	
	# output file
	file = open('output_files/other/gene_phylop.txt', "w")
	header = "locus	median_phylop_score"
	file.write(header + '\n')
	
	for key in gene_dict:
		phylop_cons_values = []
		if key:  # to skip non-coding
			for base in gene_dict[key]:
				if base != '3107':  # this position is not listed in the input file (as Ref=N), so skip
					phylop_cons_values.append(dict[base])
			file.write(key + '\t' + str(median(phylop_cons_values)) + '\n')


def codon_usage(input_file: str):
	"""Calculate tRNA codon usage, to compare to observed:expected ratio values.

	:param input_file: annotated file with mutation scores and observed maximum heteroplasmy
	"""
	# count the codons for each tRNA in the mtDNA
	dict = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		for index, ref_aa in enumerate(row["amino_acids"].split(',')):  # some positions have two since in two genes
			ref_aa = 'MT-T' + str(list(row["amino_acids"].split(","))[index].split("/")[0])
			if ref_aa == "MT-TL" or ref_aa == "MT-TS":  # these amino acids have two tRNAs
				ref_codon = str(list(row["codon_change"].split(","))[index].split("/")[0]).upper()
				if ref_aa == "MT-TL" and ref_codon.startswith('T'):
					ref_aa = "MT-TL1"
				elif ref_aa == "MT-TL" and ref_codon.startswith('C'):
					ref_aa = "MT-TL2"
				elif ref_aa == "MT-TS" and ref_codon.startswith('T'):
					ref_aa = "MT-TS1"
				elif ref_aa == "MT-TS" and ref_codon.startswith('A'):
					ref_aa = "MT-TS2"
			if ref_aa in dict:
				dict[ref_aa] += 1
			else:
				dict[ref_aa] = 1
	
	# output file
	file = open('output_files/other/tRNA_codon_counts.txt', "w")
	header = "tRNA_codon	count"
	file.write(header + '\n')
	for key in dict:
		if key != "MT-T" and key != "MT-T*" and key != "MT-TX":  # these are non-relevant fields
			file.write(key + '\t' + str(dict[key] / 9) + '\n')  # divide by 9 since synthetic vcf has each codon 9x


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	args = parser.parse_args()
	
	# set default
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	
	for path in ['output_files/other']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	gene_conservation(input_file=args.input)
	codon_usage(input_file=args.input)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
