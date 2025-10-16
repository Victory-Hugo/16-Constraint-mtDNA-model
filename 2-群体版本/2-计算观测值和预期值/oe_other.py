import argparse
import datetime
from oe_functions import *
import os
'''
该脚本用于计算线粒体 DNA（mtDNA）不同功能注释、结构域、预测工具和疾病相关变异的观测值与预期值之比（obs/exp），并输出每类的90%置信区间。主要功能包括：
1. protein_annot_oe: 基于 UniProt 等蛋白注释，统计错义变异的 obs/exp。
2. tRNA_pos_oe: 按 tRNA 位置统计 obs/exp。
3. RNA_base_types_oe: 按 RNA 基因中不同碱基类型及修饰情况统计 obs/exp。
4. tRNA_domain_oe: 按 tRNA 结构域统计 obs/exp。
5. insilico_oe: 按 in silico 预测工具（APOGEE、MitoTip、HmtVar）分类统计 obs/exp。
6. disease_vars_oe: 按疾病相关变异（MITOMAP、ClinVar）分类统计 obs/exp。
7. vus_oe: 按 VUS 子集（MITOMAP、ClinVar）统计 obs/exp。
脚本通过命令行参数指定输入文件、观测值列名、线性模型参数文件、输出前缀及需排除的碱基位点。各函数均会输出对应类别的 obs/exp 结果至指定文件夹。适用于线粒体变异功能注释与 ACMG 证据整合分析。

'''

def protein_annot_oe(
		input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算基于 UniProt 等注释的错义变异的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%sprotein_annotation_obs_exp.txt' % output_prefix, "w")
	header = "protein_annotation	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 评估功能位点，以及所有蛋白均提供的跨膜拓扑域
	prot_sum = initialize_sum_dict(identifier_list=["site", "transmembrane", "not-transmembrane"])
	
	# 针对每个注释，汇总错义变异的观测值与似然得分
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			# 仅保留错义作为最严重后果的情况
			if ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
				if ('site' in row["uniprot_annotation"]) or row["other_prot_annotation"]:
					prot_sum = sum_obs_likelihood(
						mutation=mutation, identifier="site", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=prot_sum)
						
	for annot in ["site"]:
		calculate_oe(item=annot, sum_dict=prot_sum, fit_parameters=fit_parameters, file=file)


def tRNA_pos_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算每个 tRNA 位置的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%stRNA_position_obs_exp.txt' % output_prefix, "w")
	header = "tRNA_position	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 创建 tRNA 位置列表
	tRNA_positions = []
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if (row["tRNA_position"] not in tRNA_positions) and (',' not in row["tRNA_position"]):
			tRNA_positions.append(row["tRNA_position"])
	tpos_sum = initialize_sum_dict(identifier_list=tRNA_positions)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
	if int(row["POS"]) not in excluded_sites:
		if row["symbol"].startswith('MT-T'):
			for tpos in row["tRNA_position"].split(','):  # 某些位置因位于两个基因而有两个条目
				tpos_sum = sum_obs_likelihood(
					mutation=mutation, identifier=tpos, region='ref_exc_ori',
					observed=row[obs_value], likelihood=row["Likelihood"],
					callable_samples=row["callable_samples"], dict=tpos_sum)
	# 因此此处会对位于两个基因中的 tRNA 位置进行双重计数
	
	for tRNA_pos in tRNA_positions:
		calculate_oe(item=tRNA_pos, sum_dict=tpos_sum, fit_parameters=fit_parameters, file=file)


def RNA_base_types_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算 RNA 基因中不同碱基类型的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%sRNA_base_types_obs_exp.txt' % output_prefix, "w")
	header = "RNA_base_type	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 构建 RNA 碱基类型列表
	RNA_types = ["WC_tRNA", "non-WC_tRNA", "loop-or-other_tRNA", "WC_rRNA", "non-WC_rRNA", "loop-or-other_rRNA"]
	RNA_type_sum = initialize_sum_dict(identifier_list=RNA_types)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			for RNA_type in row["RNA_base_type"].split(','):  # 遍历所有类型
				if row["symbol"].startswith('MT-T'):
					RNA_type_sum = sum_obs_likelihood(
						mutation=mutation, identifier=RNA_type + "_tRNA", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=RNA_type_sum)
				if row["symbol"].startswith('MT-R'):
					RNA_type_sum = sum_obs_likelihood(
						mutation=mutation, identifier=RNA_type + "_rRNA", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=RNA_type_sum)
	
	for RNA_type in RNA_types:
		calculate_oe(item=RNA_type, sum_dict=RNA_type_sum, fit_parameters=fit_parameters, file=file)
	
	# 同时统计 RNA 修饰碱基
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
				observed=row[obs_value], likelihood=row["Likelihood"],
				callable_samples=row["callable_samples"], dict=RNA_mod_sum)
	
	for mod in RNA_mods:
		calculate_oe(item=mod, sum_dict=RNA_mod_sum, fit_parameters=fit_parameters, file=file)

	
def tRNA_domain_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算不同 tRNA 结构域的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%stRNA_domains_obs_exp.txt' % output_prefix, "w")
	header = "tRNA_domain	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 创建 tRNA 结构域列表
	tRNA_domains = []
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if (row["tRNA_domain"] not in tRNA_domains) and (',' not in row["tRNA_domain"]):
			tRNA_domains.append(row["tRNA_domain"])
	tRNA_dom_sum = initialize_sum_dict(identifier_list=tRNA_domains)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
	if int(row["POS"]) not in excluded_sites:
		if row["symbol"].startswith('MT-T'):
			for tRNA_dom in row["tRNA_domain"].split(','):  # 遍历所有结构域
				tRNA_dom_sum = sum_obs_likelihood(
					mutation=mutation, identifier=tRNA_dom, region='ref_exc_ori',
					observed=row[obs_value], likelihood=row["Likelihood"],
					callable_samples=row["callable_samples"], dict=tRNA_dom_sum)
	
	for tRNA_dom in tRNA_domains:
		calculate_oe(item=tRNA_dom, sum_dict=tRNA_dom_sum, fit_parameters=fit_parameters, file=file)
	

def insilico_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算不同 in silico 预测工具的观测值与预期值之比及其90%置信区间。
	APOGEE 针对错义变异，MitoTip 与 HmtVar 针对 tRNA，这三者均为 ACMG 调整建议中的推荐工具。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%sinsilicos_obs_exp.txt' % output_prefix, "w")
	header = "prediction	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 先处理 APOGEE，构建分类列表
	apogee = ["Pathogenic-APOGEE", "Neutral-APOGEE"]
	apogee_sum = initialize_sum_dict(identifier_list=apogee)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
	if int(row["POS"]) not in excluded_sites:
		# 仅保留错义作为最严重后果的情况
		if ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
			for apogee_class in row["apogee_class"].split(','):  # 某些位置因位于两个基因而出现两次
				if apogee_class:
					apogee_sum = sum_obs_likelihood(
						mutation=mutation, identifier=apogee_class + "-APOGEE", region='ref_exc_ori',
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=apogee_sum)
	
	for apogee_class in apogee:
		calculate_oe(item=apogee_class, sum_dict=apogee_sum, fit_parameters=fit_parameters, file=file)
		
	# 然后处理 MitoTip
	mitotip = [
		'likely pathogenic-MitoTip', 'possibly pathogenic-MitoTip', 'possibly benign-MitoTip', 'likely benign-MitoTip']
	mitotip_sum = initialize_sum_dict(identifier_list=mitotip)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
	if int(row["POS"]) not in excluded_sites:
		if row["symbol"].startswith('MT-T'):
			mitotip_sum = sum_obs_likelihood(
				mutation=mutation, identifier=row["mitotip_class"] + "-MitoTip", region='ref_exc_ori',
				observed=row[obs_value], likelihood=row["Likelihood"],
				callable_samples=row["callable_samples"], dict=mitotip_sum)
	
	for mitotip_class in mitotip:
		calculate_oe(item=mitotip_class, sum_dict=mitotip_sum, fit_parameters=fit_parameters, file=file)
	
	# 接着处理 HmtVar
	hmtvar = ['pathogenic-HmtVar', 'likely_pathogenic-HmtVar', 'likely_polymorphic-HmtVar', 'polymorphic-HmtVar']
	hmtvar_sum = initialize_sum_dict(identifier_list=hmtvar)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
	if int(row["POS"]) not in excluded_sites:
		if row["symbol"].startswith('MT-T') and row["hmtvar_class"] != "None":
			hmtvar_sum = sum_obs_likelihood(
				mutation=mutation, identifier=row["hmtvar_class"] + "-HmtVar", region='ref_exc_ori',
				observed=row[obs_value], likelihood=row["Likelihood"],
				callable_samples=row["callable_samples"], dict=hmtvar_sum)
	
	for hmtvar_class in hmtvar:
		calculate_oe(item=hmtvar_class, sum_dict=hmtvar_sum, fit_parameters=fit_parameters, file=file)
			

def disease_vars_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算不同疾病相关变异的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%sdisease_variants_obs_exp.txt' % output_prefix, "w")
	header = "classification	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 先处理 MITOMAP
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
				observed=row[obs_value], likelihood=row["Likelihood"],
				callable_samples=row["callable_samples"], dict=mitomap_sum)
	
	for status in mitomap_status:
		calculate_oe(item=status, sum_dict=mitomap_sum, fit_parameters=fit_parameters, file=file)
	
	# 接着处理 ClinVar
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
				observed=row[obs_value], likelihood=row["Likelihood"],
				callable_samples=row["callable_samples"], dict=clinvar_sum)
	
	for status in clinvar_status:
		calculate_oe(item=status, sum_dict=clinvar_sum, fit_parameters=fit_parameters, file=file)


def vus_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算不同 VUS 子集的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%svus_obs_exp.txt' % output_prefix, "w")
	header = "classification	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	# 先处理 MITOMAP
	mitomap_status = ["MITOMAP-Reported", "MITOMAP-Reported-PP3-PM2s", "MITOMAP-Reported-BP4-BS1"]
	mitomap_sum = initialize_sum_dict(identifier_list=mitomap_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			status = ''
			if "Reported" in row["mitomap_status"]:
				status = "MITOMAP-Reported"
				if row["consequence"] == "missense_variant":  # 仅考虑错义
					if (row["apogee_class"] == "Pathogenic") and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif (row["apogee_class"] == "Neutral") and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
				if row["symbol"].startswith('MT-T'):  # 仅针对 tRNA 子集
					if ("pathogenic" in row["mitotip_class"]) and ("pathogenic" in row["hmtvar_class"])\
							and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif ("benign" in row["mitotip_class"]) and ("polymorphic" in row["hmtvar_class"])\
							and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
	if status != '':
		mitomap_sum = sum_obs_likelihood(
			mutation=mutation, identifier=status, region=region_to_use,
			observed=row[obs_value], likelihood=row["Likelihood"],
			callable_samples=row["callable_samples"], dict=mitomap_sum)
	
	for status in mitomap_status:
		calculate_oe(item=status, sum_dict=mitomap_sum, fit_parameters=fit_parameters, file=file)
	
	# 接着处理 ClinVar
	clinvar_status = ["ClinVar-Uncertain significance", "ClinVar-Uncertain significance-PP3-PM2s"]
	# 手动移除 "Uncertain significance-ClinVar-BP4-BS1"，因为没有记录满足该条件，此举用于快速规避报错
	clinvar_sum = initialize_sum_dict(identifier_list=clinvar_status)
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		if int(row["POS"]) not in excluded_sites:
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			status = ''
			if row["clinvar_interp"] == "Uncertain significance":
				status = "ClinVar-Uncertain significance"
				if row["consequence"] == "missense_variant":  # 仅考虑错义
					if (row["apogee_class"] == "Pathogenic") and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif (row["apogee_class"] == "Neutral") and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
				if row["symbol"].startswith('MT-T'):  # 仅针对 tRNA 子集
					if ("pathogenic" in row["mitotip_class"]) and ("pathogenic" in row["hmtvar_class"]) and (float(row["mitomap_af"]) < 0.00002):
						status = status + "-PP3-PM2s"
					elif ("benign" in row["mitotip_class"]) and ("polymorphic" in row["hmtvar_class"]) and (float(row["mitomap_af"]) > 0.005):
						status = status + "-BP4-BS1"
	if status != '':
		clinvar_sum = sum_obs_likelihood(
			mutation=mutation, identifier=status, region=region_to_use,
			observed=row[obs_value], likelihood=row["Likelihood"],
			callable_samples=row["callable_samples"], dict=clinvar_sum)
	
	for status in clinvar_status:
		calculate_oe(item=status, sum_dict=clinvar_sum, fit_parameters=fit_parameters, file=file)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed carrier counts")
	parser.add_argument(
		"-obs", type=str, help="Population dataset providing observed carrier counts")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-prefix", type=str, help="Output files prefix")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calibration")
	args = parser.parse_args()
	
	# 设置 gnomAD 的默认值
if args.input is None:
	args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
if args.obs is None:
	args.obs = "carrier_count"
if args.parameters is None:
	args.parameters = 'output/calibration/linear_model_fits.txt'
if args.prefix is None:
	args.prefix = ""
if args.exc_sites is None:
	# 排除 gnomAD 中的“artifact_prone_sites”：301、302、310、316、3107 和 16182（3107 已排除）
	# 这些位点在 gnomAD 中未调用，因此在计算中剔除
	args.exc_sites = [301, 302, 310, 316, 16182]
	
	for path in ['output/oe']:
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
