import argparse
import csv
import datetime
import json
import sys
sys.path.append('../mitochondrial_constraint/')
from build_model.compile_denovo import rcrs_pos_to_ref


def rcrs_pos_to_trinucleotide():
	"""Generate dictionary linking each position to its reference trinucleotide in the rCRS.

	:return: dictionary where the key is the position in rCRS, and the value is its reference trinucleotide
	"""
	# first, generate dictionary to convert coordinate to reference nucleotide
	rcrs_pos2ref = rcrs_pos_to_ref()
	# now, generate dictionary of coordinate to trinucleotide
	dict = {}
	for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'), delimiter='\t'):
		pos = int(row["POS"])
		ref = row["REF"]

		if pos == 16569:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(1)]  # dealing with circular genome
		elif pos == 1:
			trinucleotide = rcrs_pos2ref[str(16569)] + ref + rcrs_pos2ref[str(pos + 1)]  # dealing with circular genome
		elif ref == "N":
			continue  # ie skip, to handle the 'N' spacer expected at m.3107
		elif rcrs_pos2ref[str(pos + 1)] == "N":
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 2)]
		elif rcrs_pos2ref[str(pos - 1)] == "N":
			trinucleotide = rcrs_pos2ref[str(pos - 2)] + ref + rcrs_pos2ref[str(pos + 1)]
		else:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 1)]
		dict[str(pos)] = trinucleotide
	return dict


def gnomad_annotate():
	"""Generate dictionary with the maximum observed heteroplasmy (max_hl) and other annotations for each variant in gnomAD.

	:return: a dictionary, where the key is a tuple of the position, ref and alt, and the values are annotations
	"""
	dict = {}
	for row in csv.DictReader(
			open('required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv'), delimiter='\t'):
		if row["filters"] == "PASS":
			dict[(row["ref"], row["position"], row["alt"])] = (row["max_observed_heteroplasmy"], row["AF_hom"], row["AF_het"], row["AC_hom"], row["AC_het"])
	return dict


def phylop_annotate():
	"""Generate dictionary with the phyloP scores for conservation, from 100 vertebrates.

	:return: a dictionary, where the key is the position, and the value is the phyloP conservation score
	"""
	dict = {}
	pos = 0  # to handle a header row
	for row in open('required_files/insilicos/chrM.phyloP100way.wigFix'):
		dict[pos] = row.replace('\n', '')
		pos += 1
	return dict


def vep_annotate():
	"""Create a dictionary of the VEP annotations for every possible single nucleotide variant in the mtDNA.

	:return: dictionary where tuple of the variant and value identifier is key, and value is list of annotations
	"""
	vep = {}
	# use vcf where variants in two genes are split across two rows, for easy parsing
	for row in csv.DictReader(
			open("required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		# gene or locus
		if (row["REF"], row["POS"], row["ALT"], "symbol") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")].append(row["SYMBOL"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")] = [row["SYMBOL"]]
		# consequences
		if (row["REF"], row["POS"], row["ALT"], "consequence") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")].append(row["Consequence"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")] = [row["Consequence"]]
		# amino acids
		if (row["REF"], row["POS"], row["ALT"], "aa") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "aa")].append(row["Amino_acids"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "aa")] = [row["Amino_acids"]]
		# protein position
		if (row["REF"], row["POS"], row["ALT"], "codon") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "codon")].append(row["Protein_position"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "codon")] = [row["Protein_position"]]
		# codon change
		if (row["REF"], row["POS"], row["ALT"], "codon_change") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "codon_change")].append(row["Codons"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "codon_change")] = [row["Codons"]]
	return vep


def tRNA_positions():
	"""Create a dictionary of the tRNA position numbers (ranging from 1-73) encoded by mtDNA.

	:return: dictionary where position in the mtDNA is key and the value is the tRNA position
	"""
	dict = {}
	for row in csv.DictReader(open('required_files/other_annotations/tRNA_positions.txt'), delimiter='\t'):
		if row["m.pos"] in dict:  # if the mtDNA position is in two tRNA genes
			dict[row["m.pos"]].append(row["tRNA_position"])
		else:
			dict[row["m.pos"]] = [row["tRNA_position"]]
	return dict


def RNA_domains_mods():
	"""Create a dictionary of modified RNA bases and tRNA domains encoded by mtDNA.

	:return: dictionary where position in the mtDNA is key and values are whether the base is modified and tRNA domain
	"""
	dict = {}
	for row in csv.DictReader(
			open('required_files/other_annotations/mito_RNA_modifications_and_domains.txt'), delimiter='\t'):
		# handle the fact few bases are in two RNA genes
		if row["MODIFICATION"] != "#N/A":
			if "," in row["GENE"]:  # bases in two RNA genes coded together
				dict[row["POS"], "modified"] = [row["MODIFICATION"], row["MODIFICATION"]]
			else:
				dict[row["POS"], "modified"] = [row["MODIFICATION"]]
		if row["DOMAIN"]:
			if "," in row["GENE"]:  # bases in two RNA genes coded together
				dict[row["POS"], "domain"] = [row["DOMAIN"], row["DOMAIN"]]
			else:
				dict[row["POS"], "domain"] = [row["DOMAIN"]]
	# Note, there are 4 discrepancies from these domain annotations and the pair type
	# at pair m.8307+8311 and m.10039+10046, in mamit-tRNAdb they are drawn as pairs, but in the pdf are in loop domain
	return dict


def RNA_base_type():
	"""Create a dictionary of RNA base types encoded by mtDNA (Watson-Crick (WC), non-WC, or loop/other).

	:return: dictionary where position in the mtDNA is key and the value is the base type
	"""
	# first, create a dictionary for bases in pairs so can look-up if WC or not
	RNA_dict = {}
	for row in csv.DictReader(open('required_files/other_annotations/all_RNA_bases.tsv'), delimiter='\t'):
		if (row["Type"] == "b") and row["Pair_coordinate"]:  # type b excludes base m.3107N, filter to bases in pairs
			# note there are four bases in two tRNAs, hence using gene as second key
			RNA_dict[(row["Genomic_coordinate"], row["file"])] = row["RNA.base"]
	# now iterate through RNA bases
	dict = {}
	for row in csv.DictReader(open('required_files/other_annotations/all_RNA_bases.tsv'), delimiter='\t'):
		if row["Type"] == "b":  # only excludes base m.3107N
			# first determine base type
			if row["Pair_coordinate"]:  # if in pair
				base1 = RNA_dict[(row["Genomic_coordinate"], row["file"])]
				base2 = RNA_dict[(row["Pair_coordinate"], row["file"])]  # pairing base
				if (base1 == 'A' and base2 == 'T') or (base1 == 'T' and base2 == 'A'):
					base_type = "WC"
				elif (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
					base_type = "WC"
				else:
					base_type = "non-WC"
			else:
				base_type = "loop-or-other"
			# now build dictionary, appending to handle bases in two genes
			if row["Genomic_coordinate"] not in dict:
				dict[row["Genomic_coordinate"]] = [base_type]
			else:
				dict[row["Genomic_coordinate"]].append(base_type)
	return dict


def uniprot_annotations():
	"""Create a dictionary of the UniProt annotations for mtDNA-encoded proteins.

	:return: dictionary where position in the mtDNA is key and the value is the UniProt annotation
	"""
	# on the UniProt website, they display annotations for 'sites' and 'topology' for each protein, if available
	# but no family or domains for mtDNA-encoded proteins in UniProt
	dict = {}
	for row in csv.reader(open('required_files/other_annotations/uniprot_beds_chrMT_combined.txt'), delimiter='\t'):
		# restrict to annotations of interest
		if any(x in row[14] for x in ["binding", "metal"]):  # metal binding or binding site
			annotation = "site:" + row[14].split("UP000005640_9606_")[1].split(".bed")[0] + "-" + row[13]
			# per UniProt: start_coord = row[1] and end_coord = row[2], but there can be intervals/blocks between these
			# number of blocks representing the annotation is row[9]
			# row[10] are the block sizes, a comma separated list
			# row[11] are block starts, a comma separated list of block offsets relative to the annotation start
			nblock = list(range(1, int(row[9]) + 1))
			for block in nblock:
				start = int(row[1]) + int(row[11].split(",")[(block - 1)]) + 1  # need the plus 1 offset
				end = start + int(row[10].split(",")[(block - 1)]) - 1  # need the minus 1 offset
				if end > int(row[2]):
					sys.exit('Problem: the predicted end coordinate is greater than the provided')
				for pos in list(range(1, 16570)):
					if (pos >= start) and (pos <= end):
						if pos in dict:
							if annotation not in dict[pos]:
								dict[pos].append(annotation)
						else:
							dict[pos] = [annotation]
	return dict


def curated_func_sites():
	"""Residues involved in complex I proton transfer, curated from PMID:32972993.

	:return: dictionary where position in the mtDNA is key and the value is the annotation label
	"""
	# file with CI proton residues are labelled by their protein position/residue number
	# so first create a lookup to go from protein and protein position to mtDNA coordinate
	res_to_pos = {}
	for row in csv.DictReader(
			open("required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		if (row["SYMBOL"], row["Protein_position"]) not in res_to_pos:
			res_to_pos[(row["SYMBOL"], row["Protein_position"])] = [int(row["POS"])]
		else:
			res_to_pos[(row["SYMBOL"], row["Protein_position"])].append(int(row["POS"]))
	
	dict = {}
	for row in csv.DictReader(
			open('required_files/other_annotations/CI_proton_residues_PMID32972993.txt'), delimiter='\t'):
		for pos in list(range(1, 16570)):
			if pos in res_to_pos[(row["locus"], row["residue"])]:
				dict[pos] = "proton-transfer-PMID32972993"
	return dict


def apogee():
	"""Create a dictionary of the APOGEE scores for missense.

	:return: dictionary where position in the mtDNA is key and the value is the APOGEE score
	"""
	dict = {}
	for row in csv.DictReader(open('required_files/insilicos/MitImpact_db_3.0.6.txt'), delimiter='\t'):
		# note, for variants in two genes this won't distangle which gene has which score
		if (row["Ref"], row["Start"], row["Alt"]) in dict:  # i.e. variant is in two genes
			dict[(row["Ref"], row["Start"], row["Alt"])].append(row["APOGEE"])
		else:
			dict[(row["Ref"], row["Start"], row["Alt"])] = [row["APOGEE"]]
	return dict


def mitotip():
	"""Create a dictionary of the MitoTip in silico scores for tRNA variants.

	:return: dictionary where position in the mtDNA is key and the value is the MitoTip classification
	"""
	# note variants in two genes have same score
	dict = {}
	for row in csv.DictReader(open('required_files/insilicos/mitotip_scores.txt'), delimiter='\t'):
		prediction = ''
		if row["Quartile"] == "Q1":
			prediction = "likely pathogenic"
		elif row["Quartile"] == "Q2":
			prediction = "possibly pathogenic"
		elif row["Quartile"] == "Q3":
			prediction = "possibly benign"
		elif row["Quartile"] == "Q4":
			prediction = "likely benign"
		dict[(row["rCRS"], row["Position"], row["Alt"])] = prediction
	return dict


def hmtvar():
	"""Create a dictionary of the HmtVar in silico scores for tRNA variants.

	:return: dictionary where the key is a tuple of the ref, position and alt, and the value is the HmtVar classification
	"""
	dict = {}
	for row in csv.DictReader(open("required_files/insilicos/hmtvar_annotations.txt"), delimiter="\t"):
		insilico = ""
		# extract the in silico prediction from the annotation
		if len(row["HmtVar"]) > 3:
			annotation = json.loads(row["HmtVar"])
			insilico = str(annotation["pathogenicity"])
		dict[(row["REF"], row["POS"], row["ALT"])] = insilico
	return dict


def in_helix():
	"""Generate dictionary with the maximum observed heteroplasmy (max_hl) for each variant in HelixMTdb.

	:return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the max_hl
	"""
	dict = {}
	for row in csv.DictReader(open('required_files/databases/HelixMTdb_20200327.tsv'), delimiter='\t'):
		if row["alleles"].count(",") == 1:  # skip multiallelic variant calls
			pos = row["locus"].split("chrM:")[1]
			ref = row["alleles"].split("\"")[1]
			alt = row["alleles"].split("\"")[3]
			
			# calculate maximum heteroplasmy, as only provided for variants observed at heteroplasmy
			if float(row["AF_hom"]) == 0:
				max_het = float(row["max_ARF"])
			else:  # then homoplasmic are observed
				max_het = 1
			dict[(ref, pos, alt)] = (max_het, row["AF_hom"], row["AF_het"])
	return dict


def mitomap():
	"""Generate dictionary with the disease association statys for each variant in MITOMAP.

	:return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the status
	"""
	dict1, dict2 = {}, {}
	for row in csv.DictReader(open('required_files/databases/MITOMAP_disease_2022-05-25.txt'), delimiter='\t'):
		if ("Cfrm" in str(row["status"])) or ("Reported" in str(row["status"])):
			if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # if SNVs
				dict1[(row["ref"], row["pos"], row["alt"])] = (
					row["status"], row["homoplasmy"], row["heteroplasmy"], row["disease"])
	
	for row in csv.DictReader(open('required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt'), delimiter='\t'):
		if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # if SNVs
			# 56910 is total number of gb sequences to convert to allele freq
			dict2[(row["ref"], row["pos"], row["alt"])] = (int(row["gbcnt"]), (int(row["gbcnt"]) / 56910))
	
	return dict1, dict2


def clinvar():
	"""Generate dictionary with the clinical significance interpretation for each variant in ClinVar.

	:return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the interpretation
	"""
	dict = {}
	for row in csv.DictReader(open('required_files/databases/clinvar_result_chrMT_SNVonly_05252022.txt'), delimiter='\t'):
		pos = row["GRCh38Location"]
		alt = row["Canonical SPDI"].split(':')[-1]
		ref = row["Canonical SPDI"].split(':')[2]
		interp = row["Clinical significance (Last reviewed)"].split('(')[0]
		if (len(ref) == 1) and (len(alt) == 1) and (ref != alt):  # if SNVs
			#if "no assertion criteria" not in row["Review status"]:
				# exclude those only listed for cancer, some are cancer and mito diseases keep those
				# there shouldn't be any however that are not 'no assertion criteria provided for release used
			# identified through manual inspection of conditions, those only annotated for cancer excluded
			if interp != "Conflicting interpretations of pathogenicity":
				if row["Condition(s)"] != "Familial colorectal cancer" and \
						row["Condition(s)"] != "Familial cancer of breast" and \
						row["Condition(s)"] != "Acute megakaryoblastic leukemia|Mediastinal germ cell tumor" and \
						row["Condition(s)"] != "Neoplasm of ovary":
					dict[(ref, pos, alt)] = interp
	return dict


def chimp_ref_lookup():
	"""Parse an alignment of the human and chimpanzee reference mtDNA sequences and determine the ancestral chimp allele.
	Note this uses a shifted version of the human reference sequence.
	
	:return: a dictionary, where the key is a tuple of the ref and position, and the value is the chimp allele
	"""
	# read in alignment file generated by the msa package in R
	# create dictionary to parse results
	dict = {}
	with open('required_files/other_annotations/human-shifted_chimp_mt_aln.txt') as file:
		# mark position in the alignment file
		h_aln_pos = 1
		c_aln_pos = 1
		while True:
			line = file.readline()
			if not line:
				break
			if line.startswith('[1]'):  # these are the human sequence in the alignment
				sequence = line.strip().split(' ')[1]
				for base in sequence:
					dict[('human', h_aln_pos)] = base
					h_aln_pos += 1
			if line.startswith('[2]'):  # these are the chimp sequence in the alignment
				sequence = line.strip().split(' ')[1]
				for base in sequence:
					dict[('chimp', c_aln_pos)] = base
					c_aln_pos += 1
	
	# now parse dict to return chimp reference allele at each position
	# this is shifted human mtDNA, so starts at m.577 - this is to match the start position of the chimpanzee reference
	new_dict = {}
	pos = 577
	for key in dict:
		if (key[0] == 'human') and (dict[key] != '-'):  # these are gaps in human alignment
			aln_pos = key[1]
			human_ref = dict[key]
			chimp_ref = dict[('chimp', aln_pos)]
			new_dict[(human_ref, pos)] = chimp_ref
			pos += 1
			if pos == 16570:  # will happen to 16569
				pos = 1  # to renumber
	return new_dict


def annotate(input_file: str):
	"""Annotate the file with all possible mitochondrial mutations and their likelihood scores.

	:param input_file: the file with mutation likelihood scores, output of composite_likelihood_mito.py
	"""
	f = open('%s_annotated.txt' % input_file.split('.txt')[0], "w")
	header = "POS	REF	ALT	Likelihood	trinucleotide	symbol	consequence	amino_acids	protein_position	codon_change	gnomad_max_hl	gnomad_af_hom	gnomad_af_het	gnomad_ac_hom	gnomad_ac_het	in_phylotree	phyloP_score	tRNA_position	tRNA_domain	RNA_base_type	RNA_modified	rRNA_bridge_base	uniprot_annotation	other_prot_annotation	apogee_class	mitotip_class	hmtvar_class	helix_max_hl	helix_af_hom	helix_af_het	mitomap_gbcnt	mitomap_af	mitomap_status	mitomap_plasmy	mitomap_disease	clinvar_interp	chimp_ref"
	f.write(header + '\n')

	# generate required dictionaries
	rcrs_pos2trinuc = rcrs_pos_to_trinucleotide()
	gnomad = gnomad_annotate()
	phylop = phylop_annotate()
	vep = vep_annotate()
	tRNA_position = tRNA_positions()
	RNA_dom_mod = RNA_domains_mods()
	RNA_type = RNA_base_type()
	uniprot = uniprot_annotations()
	other_prot = curated_func_sites()
	apogee_scores = apogee()
	mitotip_scores = mitotip()
	hmtvar_scores = hmtvar()
	helix = in_helix()
	mitomap_vars1, mitomap_vars2 = mitomap()
	clinvar_vars = clinvar()
	chimp_dict = chimp_ref_lookup()
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		variant = row["REF"] + row["POS"] + row["ALT"]
		var_tuple = (row["REF"], row["POS"], row["ALT"])
		# annotate haplogroup variants, use '\n' in search to ensure unique matches
		in_phylo = 1 if ("\n" + variant + "\n") in open('required_files/databases/phylotree_variants.txt').read() else 0
		max_hl = gnomad[var_tuple][0] if var_tuple in gnomad else 0
		gnomad_af_hom = gnomad[var_tuple][1] if var_tuple in gnomad else 0
		gnomad_af_het = gnomad[var_tuple][2] if var_tuple in gnomad else 0
		gnomad_ac_hom = gnomad[var_tuple][3] if var_tuple in gnomad else 0
		gnomad_ac_het = gnomad[var_tuple][4] if var_tuple in gnomad else 0
		tRNA_pos = tRNA_position[row["POS"]] if row["POS"] in tRNA_position else ''
		tRNA_dom = RNA_dom_mod[(row["POS"], "domain")] if (row["POS"], "domain") in RNA_dom_mod else ''
		RNA_mod = RNA_dom_mod[(row["POS"], "modified")] if (row["POS"], "modified") in RNA_dom_mod else ''
		RNA_base = str(RNA_type[row["POS"]]).strip('[]').replace("'", "").replace(" ", "") if row["POS"] in RNA_type else ''
		RNA_bridge = "Yes" if ("\n" + row["POS"] + "\n") in open('required_files/other_annotations/rRNA_bridge_bases.txt').read() else "No"
		uniprot_annot = str(uniprot[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in uniprot else ''
		other_prot_annot = str(other_prot[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in other_prot else ''
		apogee_score = str(apogee_scores[var_tuple]).strip('[]').replace("'", "").replace(" ", "") \
			if var_tuple in apogee_scores else ''
		mitotip_score = mitotip_scores[var_tuple] if var_tuple in mitotip_scores else ''
		helix_max_hl = helix[var_tuple][0] if var_tuple in helix else 0
		helix_af_hom = helix[var_tuple][1] if var_tuple in helix else 0
		helix_af_het = helix[var_tuple][2] if var_tuple in helix else 0
		mitomap_ac = mitomap_vars2[var_tuple][0] if var_tuple in mitomap_vars2 else 0
		mitomap_af = mitomap_vars2[var_tuple][1] if var_tuple in mitomap_vars2 else 0
		mitomap_status = mitomap_vars1[var_tuple][0] if var_tuple in mitomap_vars1 else ''
		mitomap_plasmy = (mitomap_vars1[var_tuple][1] + '/' + mitomap_vars1[var_tuple][2]) if var_tuple in mitomap_vars1 else ''
		mitomap_dz = mitomap_vars1[var_tuple][3] if var_tuple in mitomap_vars1 else ''
		clinvar_int = clinvar_vars[var_tuple] if var_tuple in clinvar_vars else ''
		
		f.write(
			row["POS"] + '\t' + row["REF"] + '\t' + row["ALT"] + '\t' +
			row["Likelihood"] + '\t' +
			rcrs_pos2trinuc[row["POS"]] + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "symbol")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "consequence")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "aa")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "codon")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "codon_change")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(max_hl) + '\t' + str(gnomad_af_hom) + '\t' + str(gnomad_af_het) + '\t' + str(gnomad_ac_hom) + '\t' + str(gnomad_ac_het) + '\t' +
			str(in_phylo) + '\t' +
			phylop[int(row["POS"])] + '\t' +
			str(tRNA_pos).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(tRNA_dom).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			RNA_base + '\t' +
			str(RNA_mod).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(RNA_bridge) + '\t' +
			uniprot_annot + '\t' +
			other_prot_annot + '\t' +
			apogee_score + '\t' +
			mitotip_score + '\t' +
			str(hmtvar_scores[(row["REF"], row["POS"], row["ALT"])]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(helix_max_hl) + '\t' + str(helix_af_hom) + '\t' + str(helix_af_het) + '\t' +
			str(mitomap_ac) + '\t' + str(mitomap_af) + '\t' +
			mitomap_status + '\t' +
			mitomap_plasmy + '\t' +
			mitomap_dz + '\t' +
			clinvar_int + '\t' +
			chimp_dict[(row["REF"], int(row["POS"]))] + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-input", type=str, help="File with mutation likelihood scores")
	args = parser.parse_args()

	# set default
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mito_mutation_likelihoods.txt'

	print(datetime.datetime.now(), "Annotating mitochondrial mutation likelihoods!" + '\n')

	annotate(input_file=args.input)

	print(datetime.datetime.now(), "Script complete!" + '\n')
