import argparse
from Bio.PDB import PDBParser, MMCIFParser
import csv
import datetime
import sys
import os

# 添加模块搜索路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '1-构建模型'))

# 现在可以正常导入模块
import annotate_mutations


def calculate_distance(input_file: str, rc_file: str):
	"""计算蛋白和rRNA三维结构中距区域约束最近的位置。

	:param input_file: 包含突变似然评分和观测携带者计数的注释文件
	:param rc_file: 含区域约束最终区间的文件
	:return: 字典，记录蛋白或rRNA每个位点到区域约束的最小距离
	"""
	# mtDNA编码亚基及其在PDB文件中的链ID
	chain_id = {
		'MT-ND1': 's', 'MT-ND2': 'i', 'MT-ND3': 'j', 'MT-ND4': 'r', 'MT-ND4L': 'k', 'MT-ND5': 'l', 'MT-ND6': 'm',
		'MT-CYB': 'J', 'MT-CO1': 'A', 'MT-CO2': 'B', 'MT-CO3': 'C', 'MT-ATP6': 'A'}
	# 初始化区域约束内的残基字典
	rc_residues = {}
	for key in chain_id:
		rc_residues[key] = []
	# 根据区间生成残基编号列表
	for row in csv.DictReader(open(rc_file), delimiter='\t'):
		if row["locus"] in rc_residues:
			if row["locus"] == 'MT-ND6':  # 反义链
				for residue in list(range(int(row["protein_pos_end"]), int(row["protein_pos_start"]) + 1)):
					rc_residues[row["locus"]].append(residue)
			else:
				for residue in list(range(int(row["protein_pos_start"]), int(row["protein_pos_end"]) + 1)):
					rc_residues[row["locus"]].append(residue)
				
	# 创建解析器
	parser = PDBParser()
	# 从复合体I开始，加载结构
	structure = parser.get_structure('complexI', '0-required_files/pdb/5xtd.pdb')
	# 对于每个残基，计算其到区域约束的最小距离
	min_distance_dict = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if row["symbol"].startswith('MT-N'):
			min_distance = 100000  # 设定一个较大的初始值
			pair = ''
			for protein in ['MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6']:  # 复合体I全部蛋白
				for residue in rc_residues[protein]:  # 遍历其中所有受限残基
					try:
						distance = \
							structure[0][chain_id[row["symbol"]]][int(row["protein_position"])]['CA'] - \
							structure[0][chain_id[protein]][residue]['CA']  # 区域约束残基
					except KeyError:
						continue
					if distance < min_distance:
						min_distance = distance
						pair = 'same' if (row["symbol"] == protein) else 'different'
			min_distance = 'NA' if (min_distance == 100000) else min_distance
			min_distance_dict[row["POS"]] = (min_distance, pair)
	
	# 复合体III
	structure = parser.get_structure('complexIII', '0-required_files/pdb/5xte.pdb')
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if row["symbol"] == "MT-CYB":
			min_distance = 100000  # 设定一个较大的初始值
			pair = 'same'
			for protein in ['MT-CYB']:
				for residue in rc_residues[protein]:
					try:
						distance = \
							structure[0][chain_id[row["symbol"]]][int(row["protein_position"])]['CA'] - \
							structure[0][chain_id[row["symbol"]]][residue]['CA']  # 区域约束残基
					except KeyError:
						continue
					if distance < min_distance:
						min_distance = distance
			min_distance = 'NA' if (min_distance == 100000) else min_distance
			min_distance_dict[row["POS"]] = (min_distance, pair)
	
	# 复合体IV
	structure = parser.get_structure('complexIV', '0-required_files/pdb/5z62.pdb')
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if row["symbol"].startswith('MT-CO'):
			min_distance = 100000  # 设定一个较大的初始值
			pair = ''
			for protein in ['MT-CO1', 'MT-CO2', 'MT-CO3']:  # mtDNA编码的复合体IV蛋白
				for residue in rc_residues[protein]:  # 遍历其中受限残基
					try:
						distance = \
							structure[0][chain_id[row["symbol"]]][int(row["protein_position"])]['CA'] - \
							structure[0][chain_id[protein]][residue]['CA']  # 区域约束残基
					except KeyError:
						continue
					if distance < min_distance:
						min_distance = distance
						pair = 'same' if (row["symbol"] == protein) else 'different'
			min_distance = 'NA' if (min_distance == 100000) else min_distance
			min_distance_dict[row["POS"]] = (min_distance, pair)
	
	# 复合体V——因缺少人类结构，使用AlphaFold模型
	structure = parser.get_structure('complexV', '0-required_files/pdb/AF-P00846-F1-model_v2.pdb')
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if row["symbol"] == "MT-ATP6":
			min_distance = 100000  # 设定一个较大的初始值
			pair = 'same'
			for protein in ['MT-ATP6']:
				for residue in rc_residues[protein]:
					try:
						distance = \
							structure[0][chain_id[row["symbol"]]][int(row["protein_position"])]['CA'] - \
							structure[0][chain_id[row["symbol"]]][residue]['CA']  # 区域约束残基
					except KeyError:
						continue
					if distance < min_distance:
						min_distance = distance
			min_distance = 'NA' if (min_distance == 100000) else min_distance
			min_distance_dict[row["POS"]] = (min_distance, pair)
	
	# rRNA部分
	parser = MMCIFParser()
	structure = parser.get_structure('ribosome', '0-required_files/pdb/6zse.cif')
	chain_id = {'MT-RNR1': 'AA', 'MT-RNR2': 'XA'}
	# 初始化区域约束内的核苷酸字典
	rc_bases = {}
	for key in chain_id:
		rc_bases[key] = []
	# 根据区间生成核苷酸编号列表
	for row in csv.DictReader(open(rc_file), delimiter='\t'):
		if row["locus"] in rc_bases:
			for base in list(range(int(row["start"]), int(row["end"]) + 1)):
				rc_bases[row["locus"]].append(base)
	
	# 先构建参考RNA碱基字典
	rna_ref = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if "MT-R" in row["symbol"]:
			rna_ref[int(row["POS"])] = row["REF"]
	# 计算到区域约束的最小距离
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if row["symbol"].startswith('MT-R'):
			min_distance = 100000  # 设定一个较大的初始值
			pair = ''
			for rna in rc_bases:
				for base in rc_bases[rna]:  # 遍历其中受限核苷酸
					try:
						# 测量参与碱基配对的氮原子
						# 参考：https://chem.libretexts.org/Courses/California_Polytechnic_State_University_San_Luis_Obispo/Survey_of_Biochemistry_and_Biotechnology/10%3A_Supplemental_Modules_(Biochemistry)/10.01%3A_DNA/10.1.02%3A_Structure_of_DNA_and_RNA
						atom1 = 'N3' if (row["REF"] == "C" or row["REF"] == "T") else 'N1'
						atom2 = 'N3' if (rna_ref[base] == "C" or rna_ref[base] == "T") else 'N1'
						distance = \
							structure[0][chain_id[row["symbol"]]][int(row["POS"])][atom1] - \
							structure[0][chain_id[rna]][base][atom2]  # 原文为'N1'
						# 同时计算磷酸骨架之间的距离，以覆盖毗邻碱基
						distance_alt = \
							structure[0][chain_id[row["symbol"]]][int(row["POS"])]['P'] - \
							structure[0][chain_id[rna]][base]['P']
						# 取两者较小值，从而兼顾碱基配对和邻近位点
						if distance_alt < distance:
							distance = distance_alt
					except KeyError:
						continue
					if distance < min_distance:
						min_distance = distance
						pair = 'same' if (row["symbol"] == rna) else 'different'
			min_distance = 'NA' if (min_distance == 100000) else min_distance
			min_distance_dict[row["POS"]] = (min_distance, pair)

	return min_distance_dict


def annotate_rc(input_file: str, rc_file: str):
	"""为所有可能的mtDNA SNV添加区域约束注释。

	:param input_file: 包含突变似然评分和观测携带者计数的注释文件
	:param rc_file: 含区域约束最终区间的文件
	"""
	f = open('output/regional_constraint/mito_regional_constraint_annotation.txt', "w")
	header = "POS	REF	ALT	symbol	in_rc	min_distance_to_rc	distance_pair	consequence	amino_acids	protein_position	codon_change	carrier_count	callable_samples	allele_frequency	gnomad_max_hl	in_phylotree	phyloP_score	tRNA_position	tRNA_domain	RNA_base_type	RNA_modified	rRNA_bridge_base	uniprot_annotation	other_prot_annotation	apogee_class	mitotip_class	hmtvar_class	helix_max_hl	mitomap_gbcnt	mitomap_af	mitomap_status	mitomap_plasmy	mitomap_disease	clinvar_interp"
	f.write(header + '\n')
	
	# 生成所需的辅助字典
	distance_dict = calculate_distance(input_file=input_file, rc_file=rc_file)
	gnomad = annotate_mutations.gnomad_annotate()
	phylop = annotate_mutations.phylop_annotate()
	vep = annotate_mutations.vep_annotate()
	tRNA_position = annotate_mutations.tRNA_positions()
	RNA_dom_mod = annotate_mutations.RNA_domains_mods()
	RNA_type = annotate_mutations.RNA_base_type()
	uniprot = annotate_mutations.uniprot_annotations()
	other_prot = annotate_mutations.curated_func_sites()
	apogee_scores = annotate_mutations.apogee()
	mitotip_scores = annotate_mutations.mitotip()
	hmtvar_scores = annotate_mutations.hmtvar()
	helix = annotate_mutations.in_helix()
	mitomap_vars1, mitomap_vars2 = annotate_mutations.mitomap()
	clinvar_vars = annotate_mutations.clinvar()
	
	# 构建区域约束范围的字典
	dict = {}
	for row in csv.DictReader(open(rc_file), delimiter='\t'):
		dict[(int(row["start"]), int(row["end"]), row["locus"])] = row["upper_CI"]

	for row in csv.DictReader(open(input_file), delimiter='\t'):
		# 判断该位点是否处于区域约束范围内
		rc = "no"
		symbol = str(vep[(row["REF"], row["POS"], row["ALT"], "symbol")]).strip('[]').replace("'", "").replace(" ", "")
		for key in dict:
			if (int(row["POS"]) >= key[0]) and (int(row["POS"]) <= key[1]) and (symbol == key[2]):
				rc = "yes"
			
		# 复用其他模块的注释信息
		variant = row["REF"] + row["POS"] + row["ALT"]
		var_tuple = (row["REF"], row["POS"], row["ALT"])
		in_phylo = 1 if ("\n" + variant + "\n") in open('0-required_files/databases/phylotree_variants.txt').read() else 0
		max_hl = gnomad[var_tuple][0] if var_tuple in gnomad else 0
		tRNA_pos = tRNA_position[row["POS"]] if row["POS"] in tRNA_position else ''
		tRNA_dom = RNA_dom_mod[(row["POS"], "domain")] if (row["POS"], "domain") in RNA_dom_mod else ''
		RNA_mod = RNA_dom_mod[(row["POS"], "modified")] if (row["POS"], "modified") in RNA_dom_mod else ''
		RNA_base = str(RNA_type[row["POS"]]).strip('[]').replace("'", "").replace(" ", "") if row["POS"] in RNA_type else ''
		RNA_bridge = "Yes" if ("\n" + row["POS"] + "\n") in open('0-required_files/other_annotations/rRNA_bridge_bases.txt').read() else "No"
		uniprot_annot = str(uniprot[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in uniprot else ''
		other_prot_annot = str(other_prot[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in other_prot else ''
		apogee_score = str(apogee_scores[var_tuple]).strip('[]').replace("'", "").replace(" ", "") \
			if var_tuple in apogee_scores else ''
		mitotip_score = mitotip_scores[var_tuple] if var_tuple in mitotip_scores else ''
		helix_max_hl = helix[var_tuple][0] if var_tuple in helix else 0
		mitomap_ac = mitomap_vars2[var_tuple][0] if var_tuple in mitomap_vars2 else 0
		mitomap_af = mitomap_vars2[var_tuple][1] if var_tuple in mitomap_vars2 else 0
		mitomap_status = mitomap_vars1[var_tuple][0] if var_tuple in mitomap_vars1 else ''
		mitomap_plasmy = (mitomap_vars1[var_tuple][1] + '/' + mitomap_vars1[var_tuple][2]) if var_tuple in mitomap_vars1 else ''
		mitomap_dz = mitomap_vars1[var_tuple][3] if var_tuple in mitomap_vars1 else ''
		clinvar_int = clinvar_vars[var_tuple] if var_tuple in clinvar_vars else ''
		
		f.write(
			row["POS"] + '\t' + row["REF"] + '\t' + row["ALT"] + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "symbol")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(rc) + '\t' +
			str(distance_dict[row["POS"]][0] if (row["POS"]) in distance_dict else '') + '\t' +
			str(distance_dict[row["POS"]][1] if (row["POS"]) in distance_dict else '') + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "consequence")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "aa")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "codon")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "codon_change")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(row["carrier_count"]) + '\t' + str(row["callable_samples"]) + '\t' + str(row["allele_frequency"]) + '\t' +
			str(max_hl) + '\t' +
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
			str(helix_max_hl) + '\t' +
			str(mitomap_ac) + '\t' + str(mitomap_af) + '\t' + mitomap_status + '\t' + mitomap_plasmy + '\t' + mitomap_dz + '\t' +
			clinvar_int + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-input", type=str, help="Annotated file with mutation likelihood scores")
	parser.add_argument("-rc_input", type=str, help="File with areas of regional constraint")
	parser.add_argument(
		"-obs", type=str, help="Population dataset providing observed carrier counts")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calibration")
	args = parser.parse_args()

	# 设置默认参数
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.rc_input is None:
		args.rc_input = 'output/regional_constraint/real_alignment/regional_constraint_intervals.txt'
	if args.obs is None:
		args.obs = "carrier_count"
	if args.parameters is None:
		args.parameters = 'output/calibration/linear_model_fits.txt'
	if args.exc_sites is None:
		# 排除gnomAD中易受技术伪影影响的位置：301、302、310、316、3107和16182（3107已在前面排除）
		# 这些位点未在gnomAD中检测到，因此需在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]

	print(datetime.datetime.now(), "Annotating mitochondrial mutations with regional constraint details!" + '\n')
	
	annotate_rc(input_file=args.input, rc_file=args.rc_input)

	print(datetime.datetime.now(), "Script complete!" + '\n')
	
