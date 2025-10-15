import argparse
import csv
import datetime
import json
import os
import sys
from compile_denovo import rcrs_pos_to_ref


def rcrs_pos_to_trinucleotide():
	"""生成一个字典，将每个位点映射到 rCRS 中的参考三核苷酸。

	:return: 字典，键为 rCRS 中的位点，值为该位点的参考三核苷酸
	"""
	# 首先生成一个字典，用于将坐标转换为参考碱基
	rcrs_pos2ref = rcrs_pos_to_ref()
	# 接着构建坐标到三核苷酸的映射字典
	dict = {}
	for row in csv.DictReader(open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'), delimiter='\t'):
		pos = int(row["POS"])
		ref = row["REF"]

		if pos == 16569:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(1)]  # 处理环状基因组
		elif pos == 1:
			trinucleotide = rcrs_pos2ref[str(16569)] + ref + rcrs_pos2ref[str(pos + 1)]  # 处理环状基因组
		elif ref == "N":
			continue  # 跳过该位点，用于处理 m.3107 处预期的“N”间隔
		elif rcrs_pos2ref[str(pos + 1)] == "N":
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 2)]
		elif rcrs_pos2ref[str(pos - 1)] == "N":
			trinucleotide = rcrs_pos2ref[str(pos - 2)] + ref + rcrs_pos2ref[str(pos + 1)]
		else:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 1)]
		dict[str(pos)] = trinucleotide
	return dict


def gnomad_annotate():
	"""生成一个字典，记录 gnomAD 中每个变异的最大观测杂合度（max_hl）及其他注释信息。

	:return: 字典，键为包含位点、参考等位基因和替换等位基因的元组，值为对应的注释
	"""
	dict = {}
	for row in csv.DictReader(
			open('0-required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv'), delimiter='\t'):
		if row["filters"] == "PASS":
			dict[(row["ref"], row["position"], row["alt"])] = (row["max_observed_heteroplasmy"], row["AF_hom"], row["AF_het"], row["AC_hom"], row["AC_het"])
	return dict


def load_carrier_stats(stats_path: str):
	"""读取 VCF 汇总得到的携带者统计信息。

	:param stats_path: 由 compile_denovo.py 生成的变异统计文件路径
	:return: 字典，键为 (REF, POS, ALT)，值为 (carrier_count, callable_samples, allele_frequency)
	"""
	dict = {}
	if (stats_path is None) or (not os.path.exists(stats_path)):
		return dict
	for row in csv.DictReader(open(stats_path), delimiter='\t'):
		dict[(row["REF"], row["POS"], row["ALT"])] = (
			int(row["carrier_count"]),
			int(row["callable_samples"]),
			float(row["allele_frequency"]),
		)
	return dict


def phylop_annotate():
	"""生成一个字典，记录来自 100 种脊椎动物的 phyloP 保守性评分。

	:return: 字典，键为位点，值为对应的 phyloP 保守性评分
	"""
	dict = {}
	pos = 0  # 用于处理包含表头的首行
	for row in open('0-required_files/insilicos/chrM.phyloP100way.wigFix'):
		dict[pos] = row.replace('\n', '')
		pos += 1
	return dict


def vep_annotate():
	"""创建一个字典，存储 mtDNA 中所有可能单核苷酸变体的 VEP 注释。

	:return: 字典，键为包含变体和字段标识的元组，值为注释列表
	"""
	vep = {}
	# 使用将落在两个基因中的变体拆分成两行的 VCF，便于解析
	for row in csv.DictReader(
			open("0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		# 基因或基因座
		if (row["REF"], row["POS"], row["ALT"], "symbol") in vep:  # 即该变体位于两个基因中
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")].append(row["SYMBOL"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")] = [row["SYMBOL"]]
		# 影响后果
		if (row["REF"], row["POS"], row["ALT"], "consequence") in vep:  # 即该变体位于两个基因中
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")].append(row["Consequence"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")] = [row["Consequence"]]
		# 氨基酸
		if (row["REF"], row["POS"], row["ALT"], "aa") in vep:  # 即该变体位于两个基因中
			vep[(row["REF"], row["POS"], row["ALT"], "aa")].append(row["Amino_acids"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "aa")] = [row["Amino_acids"]]
		# 蛋白质位置
		if (row["REF"], row["POS"], row["ALT"], "codon") in vep:  # 即该变体位于两个基因中
			vep[(row["REF"], row["POS"], row["ALT"], "codon")].append(row["Protein_position"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "codon")] = [row["Protein_position"]]
		# 密码子变化
		if (row["REF"], row["POS"], row["ALT"], "codon_change") in vep:  # 即该变体位于两个基因中
			vep[(row["REF"], row["POS"], row["ALT"], "codon_change")].append(row["Codons"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "codon_change")] = [row["Codons"]]
	return vep


def tRNA_positions():
	"""创建一个字典，记录 mtDNA 编码的 tRNA 位置编号（范围 1-73）。

	:return: 字典，键为 mtDNA 中的位点，值为对应的 tRNA 位置
	"""
	dict = {}
	for row in csv.DictReader(open('0-required_files/other_annotations/tRNA_positions.txt'), delimiter='\t'):
		if row["m.pos"] in dict:  # 如果该 mtDNA 位点属于两个 tRNA 基因
			dict[row["m.pos"]].append(row["tRNA_position"])
		else:
			dict[row["m.pos"]] = [row["tRNA_position"]]
	return dict


def RNA_domains_mods():
	"""创建一个字典，记录 mtDNA 编码的 RNA 修饰碱基与 tRNA 结构域。

	:return: 字典，键为 mtDNA 位点，值表示碱基是否被修饰以及所属的 tRNA 结构域
	"""
	dict = {}
	for row in csv.DictReader(
			open('0-required_files/other_annotations/mito_RNA_modifications_and_domains.txt'), delimiter='\t'):
		# 处理少数位于两个 RNA 基因中的碱基
		if row["MODIFICATION"] != "#N/A":
			if "," in row["GENE"]:  # 同时属于两个 RNA 基因的碱基
				dict[row["POS"], "modified"] = [row["MODIFICATION"], row["MODIFICATION"]]
			else:
				dict[row["POS"], "modified"] = [row["MODIFICATION"]]
		if row["DOMAIN"]:
			if "," in row["GENE"]:  # 同时属于两个 RNA 基因的碱基
				dict[row["POS"], "domain"] = [row["DOMAIN"], row["DOMAIN"]]
			else:
				dict[row["POS"], "domain"] = [row["DOMAIN"]]
	# 注意，这些结构域注释与配对类型存在 4 处差异：
	# 在 m.8307+8311 和 m.10039+10046 配对位点，mamit-tRNAdb 将其绘制为配对结构，但 PDF 中标注为环区
	return dict


def RNA_base_type():
	"""创建一个字典，记录 mtDNA 编码的 RNA 碱基类型（Watson-Crick (WC)、非 WC 或环/其他）。

	:return: 字典，键为 mtDNA 中的位点，值为对应的碱基类型
	"""
	# 首先为成对碱基创建查找表，用于判断是否为 WC 配对
	RNA_dict = {}
	for row in csv.DictReader(open('0-required_files/other_annotations/all_RNA_bases.tsv'), delimiter='\t'):
		if (row["Type"] == "b") and row["Pair_coordinate"]:  # 类型 b 排除 m.3107N，仅保留成对碱基
			# 注意有四个碱基落在两个 tRNA 中，因此使用基因名作为第二个键
			RNA_dict[(row["Genomic_coordinate"], row["file"])] = row["RNA.base"]
	# 然后遍历 RNA 碱基
	dict = {}
	for row in csv.DictReader(open('0-required_files/other_annotations/all_RNA_bases.tsv'), delimiter='\t'):
		if row["Type"] == "b":  # 仅排除 m.3107N
			# 先确定碱基类型
			if row["Pair_coordinate"]:  # 如果存在配对
				base1 = RNA_dict[(row["Genomic_coordinate"], row["file"])]
				base2 = RNA_dict[(row["Pair_coordinate"], row["file"])]  # 配对碱基
				if (base1 == 'A' and base2 == 'T') or (base1 == 'T' and base2 == 'A'):
					base_type = "WC"
				elif (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
					base_type = "WC"
				else:
					base_type = "non-WC"
			else:
				base_type = "loop-or-other"
			# 构建字典，如遇位于两个基因中的碱基则追加
			if row["Genomic_coordinate"] not in dict:
				dict[row["Genomic_coordinate"]] = [base_type]
			else:
				dict[row["Genomic_coordinate"]].append(base_type)
	return dict


def uniprot_annotations():
	"""创建一个字典，记录 mtDNA 编码蛋白的 UniProt 注释。

	:return: 字典，键为 mtDNA 位点，值为对应的 UniProt 注释
	"""
	# 在 UniProt 网站上，如果有可用信息，会展示蛋白的“sites”和“topology”注释
	# 但 UniProt 中不提供 mtDNA 编码蛋白的家族或结构域信息
	dict = {}
	for row in csv.reader(open('0-required_files/other_annotations/uniprot_beds_chrMT_combined.txt'), delimiter='\t'):
		# 仅保留感兴趣的注释
		if any(x in row[14] for x in ["binding", "metal"]):  # 金属结合或结合位点
			annotation = "site:" + row[14].split("UP000005640_9606_")[1].split(".bed")[0] + "-" + row[13]
			# 按 UniProt 说明：start_coord = row[1]，end_coord = row[2]，但二者之间可能包含多个区间/区块
			# 表示注释的区块数量位于 row[9]
			# row[10] 为区块长度，使用逗号分隔
			# row[11] 为区块起点，相对于注释起点的偏移量，逗号分隔
			nblock = list(range(1, int(row[9]) + 1))
			for block in nblock:
				start = int(row[1]) + int(row[11].split(",")[(block - 1)]) + 1  # 需加 1 的偏移量
				end = start + int(row[10].split(",")[(block - 1)]) - 1  # 需减 1 的偏移量
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
	"""整理自 PMID:32972993 的复合物 I 质子转运相关残基。

	:return: 字典，键为 mtDNA 位点，值为对应的注释标签
	"""
	# 质子转运残基文件以蛋白质位置/残基编号标注
	# 因此先创建一个查找表，将蛋白及其残基位置映射到 mtDNA 坐标
	res_to_pos = {}
	for row in csv.DictReader(
			open("0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		if (row["SYMBOL"], row["Protein_position"]) not in res_to_pos:
			res_to_pos[(row["SYMBOL"], row["Protein_position"])] = [int(row["POS"])]
		else:
			res_to_pos[(row["SYMBOL"], row["Protein_position"])].append(int(row["POS"]))
	
	dict = {}
	for row in csv.DictReader(
			open('0-required_files/other_annotations/CI_proton_residues_PMID32972993.txt'), delimiter='\t'):
		for pos in list(range(1, 16570)):
			if pos in res_to_pos[(row["locus"], row["residue"])]:
				dict[pos] = "proton-transfer-PMID32972993"
	return dict


def apogee():
	"""创建一个字典，记录错义突变的 APOGEE 评分。

	:return: 字典，键为 mtDNA 位点，值为对应的 APOGEE 评分
	"""
	dict = {}
	for row in csv.DictReader(open('0-required_files/insilicos/MitImpact_db_3.0.6.txt'), delimiter='\t'):
		# 注意：对于位于两个基因中的变体，无法区分每个基因对应的评分
		if (row["Ref"], row["Start"], row["Alt"]) in dict:  # 即变体位于两个基因中
			dict[(row["Ref"], row["Start"], row["Alt"])].append(row["APOGEE"])
		else:
			dict[(row["Ref"], row["Start"], row["Alt"])] = [row["APOGEE"]]
	return dict


def mitotip():
	"""创建一个字典，记录 tRNA 变体的 MitoTip 预测评分。

	:return: 字典，键为 mtDNA 位点，值为对应的 MitoTip 分类
	"""
	# 注意：位于两个基因中的变体拥有相同的评分
	dict = {}
	for row in csv.DictReader(open('0-required_files/insilicos/mitotip_scores.txt'), delimiter='\t'):
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
	"""创建一个字典，记录 tRNA 变体的 HmtVar 预测评分。

	:return: 字典，键为包含 ref、位点和 alt 的元组，值为对应的 HmtVar 分类
	"""
	dict = {}
	for row in csv.DictReader(open("0-required_files/insilicos/hmtvar_annotations.txt"), delimiter="\t"):
		insilico = ""
		# 从注释中提取体外预测结果
		if len(row["HmtVar"]) > 3:
			annotation = json.loads(row["HmtVar"])
			insilico = str(annotation["pathogenicity"])
		dict[(row["REF"], row["POS"], row["ALT"])] = insilico
	return dict


def in_helix():
	"""生成一个字典，记录 HelixMTdb 中每个变体的最大观测杂合度（max_hl）。

	:return: 字典，键为包含 ref、位点和 alt 的元组，值为对应的 max_hl
	"""
	dict = {}
	for row in csv.DictReader(open('0-required_files/databases/HelixMTdb_20200327.tsv'), delimiter='\t'):
		if row["alleles"].count(",") == 1:  # 跳过多等位位点
			pos = row["locus"].split("chrM:")[1]
			ref = row["alleles"].split("\"")[1]
			alt = row["alleles"].split("\"")[3]
			
			# 计算最大杂合度，仅对存在杂合度的数据提供
			if float(row["AF_hom"]) == 0:
				max_het = float(row["max_ARF"])
			else:  # 出现纯合观察
				max_het = 1
			dict[(ref, pos, alt)] = (max_het, row["AF_hom"], row["AF_het"])
	return dict


def mitomap():
	"""生成一个字典，记录 MITOMAP 中每个变体的疾病关联状态。

	:return: 字典，键为包含 ref、位点和 alt 的元组，值为对应的注释状态
	"""
	dict1, dict2 = {}, {}
	for row in csv.DictReader(open('0-required_files/databases/MITOMAP_disease_2022-05-25.txt'), delimiter='\t'):
		if ("Cfrm" in str(row["status"])) or ("Reported" in str(row["status"])):
			if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # 若为 SNV
				dict1[(row["ref"], row["pos"], row["alt"])] = (
					row["status"], row["homoplasmy"], row["heteroplasmy"], row["disease"])
	
	for row in csv.DictReader(open('0-required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt'), delimiter='\t'):
		if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # 若为 SNV
			# 56910 为总的 gb 序列数，用于换算等位基因频率
			dict2[(row["ref"], row["pos"], row["alt"])] = (int(row["gbcnt"]), (int(row["gbcnt"]) / 56910))
	
	return dict1, dict2


def clinvar():
	"""生成一个字典，记录 ClinVar 中每个变体的临床意义解读。

	:return: 字典，键为包含 ref、位点和 alt 的元组，值为对应的解读
	"""
	dict = {}
	for row in csv.DictReader(open('0-required_files/databases/clinvar_result_chrMT_SNVonly_05252022.txt'), delimiter='\t'):
		pos = row["GRCh38Location"]
		alt = row["Canonical SPDI"].split(':')[-1]
		ref = row["Canonical SPDI"].split(':')[2]
		interp = row["Clinical significance (Last reviewed)"].split('(')[0]
		if (len(ref) == 1) and (len(alt) == 1) and (ref != alt):  # 若为 SNV
			#if "no assertion criteria" not in row["Review status"]:
				# 排除仅标注为癌症的条目；若兼有癌症与线粒体疾病则保留
				# 正常情况下不会出现未提供“no assertion criteria”的条目
			# 通过人工检查条件，排除仅注释为癌症的条目
			if interp != "Conflicting interpretations of pathogenicity":
				if row["Condition(s)"] != "Familial colorectal cancer" and \
						row["Condition(s)"] != "Familial cancer of breast" and \
						row["Condition(s)"] != "Acute megakaryoblastic leukemia|Mediastinal germ cell tumor" and \
						row["Condition(s)"] != "Neoplasm of ovary":
					dict[(ref, pos, alt)] = interp
	return dict


def chimp_ref_lookup():
	"""解析人类与黑猩猩参考 mtDNA 序列的比对，并确定祖先黑猩猩等位基因。
	注意，这里使用的是经过偏移的人类参考序列。
	
	:return: 字典，键为包含参考碱基与位点的元组，值为黑猩猩等位基因
	"""
	# 读取由 R 中 msa 包生成的比对文件
	# 创建字典以便解析结果
	dict = {}
	with open('0-required_files/other_annotations/human-shifted_chimp_mt_aln.txt') as file:
		# 在比对文件中记录位置
		h_aln_pos = 1
		c_aln_pos = 1
		while True:
			line = file.readline()
			if not line:
				break
			if line.startswith('[1]'):  # 比对中对应人类序列的行
				sequence = line.strip().split(' ')[1]
				for base in sequence:
					dict[('human', h_aln_pos)] = base
					h_aln_pos += 1
			if line.startswith('[2]'):  # 比对中对应黑猩猩序列的行
				sequence = line.strip().split(' ')[1]
				for base in sequence:
					dict[('chimp', c_aln_pos)] = base
					c_aln_pos += 1
	
	# 解析字典以返回每个位点的黑猩猩参考等位基因
	# 由于采用移位后的人类 mtDNA，起始位置为 m.577，以匹配黑猩猩参考序列的起点
	new_dict = {}
	pos = 577
	for key in dict:
		if (key[0] == 'human') and (dict[key] != '-'):  # 排除人类序列中的缺口
			aln_pos = key[1]
			human_ref = dict[key]
			chimp_ref = dict[('chimp', aln_pos)]
			new_dict[(human_ref, pos)] = chimp_ref
			pos += 1
			if pos == 16570:  # 处理 16569 之后的循环
				pos = 1  # 重新编号
	return new_dict


def annotate(input_file: str, variant_stats: str):
	"""为文件添加所有可能的线粒体突变及其可能性评分注释。

	:param input_file: 突变可能性评分文件，即 composite_likelihood_mito.py 的输出
	:param variant_stats: 携带者统计文件路径，由 compile_denovo.py 生成
	"""
	f = open('%s_annotated.txt' % input_file.split('.txt')[0], "w")
	header = "POS	REF	ALT	Likelihood	trinucleotide	symbol	consequence	amino_acids	protein_position	codon_change	carrier_count	callable_samples	allele_frequency	gnomad_max_hl	gnomad_af_hom	gnomad_af_het	gnomad_ac_hom	gnomad_ac_het	in_phylotree	phyloP_score	tRNA_position	tRNA_domain	RNA_base_type	RNA_modified	rRNA_bridge_base	uniprot_annotation	other_prot_annotation	apogee_class	mitotip_class	hmtvar_class	helix_max_hl	helix_af_hom	helix_af_het	mitomap_gbcnt	mitomap_af	mitomap_status	mitomap_plasmy	mitomap_disease	clinvar_interp	chimp_ref"
	f.write(header + '\n')

	# 生成所需的字典
	rcrs_pos2trinuc = rcrs_pos_to_trinucleotide()
	gnomad = gnomad_annotate()
	carriers = load_carrier_stats(variant_stats)
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
		# 注释单倍群变异，使用 '\n' 确保匹配唯一
		in_phylo = 1 if ("\n" + variant + "\n") in open('0-required_files/databases/phylotree_variants.txt').read() else 0
		max_hl = gnomad[var_tuple][0] if var_tuple in gnomad else 0
		gnomad_af_hom = gnomad[var_tuple][1] if var_tuple in gnomad else 0
		gnomad_af_het = gnomad[var_tuple][2] if var_tuple in gnomad else 0
		gnomad_ac_hom = gnomad[var_tuple][3] if var_tuple in gnomad else 0
		gnomad_ac_het = gnomad[var_tuple][4] if var_tuple in gnomad else 0
		carrier_count, callable_samples, allele_frequency = carriers.get(var_tuple, (0, 0, 0.0))
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
			str(carrier_count) + '\t' + str(callable_samples) + '\t' + str(allele_frequency) + '\t' +
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
	parser.add_argument(
		"--variant-stats",
		dest="variant_stats",
		type=str,
		help="Carrier statistics generated by compile_denovo.py",
	)
	args = parser.parse_args()

	# 设置默认输入文件
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods.txt'
	if args.variant_stats is None:
		args.variant_stats = 'output/denovo/variant_carrier_stats.txt'

	print(datetime.datetime.now(), "Annotating mitochondrial mutation likelihoods!" + '\n')

	annotate(input_file=args.input, variant_stats=args.variant_stats)

	print(datetime.datetime.now(), "Script complete!" + '\n')
