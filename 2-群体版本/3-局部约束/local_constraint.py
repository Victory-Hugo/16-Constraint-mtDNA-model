import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
from numpy import mean
import os
import pandas as pd
import sys
from typing import Dict, List, Tuple
'''
该脚本用于分析线粒体DNA（mtDNA）上的局部功能约束，计算每个窗口（kmer）和每个位点的观测/期望（obs/exp）比值及其置信区间，并进行功能注释。主要流程包括：
1. 构建查找字典：从注释文件中提取每个位点和突变的观测最大杂合度及似然值。
2. 并行计算窗口（kmer）：对每个窗口计算obs/exp比值及置信区间，并注释重叠的MITOMAP功能区。
3. 汇总窗口结果：生成所有窗口的详细注释表，并按不同功能区（蛋白、RNA、非编码）进行分组排名。
4. 基于窗口结果，计算每个位点的平均obs/exp置信区间及相关注释，输出每个位点的功能约束评分。
5. 最终输出窗口和每个位点的局部约束结果文件，便于后续分析和可视化。
主要函数说明：
- lookup_dict: 构建用于查找观测值和似然值的字典。
- parallelize_kmers: 并行计算每个窗口的obs/exp比值及置信区间，并进行功能区注释。
- kmers: 计算所有窗口的obs/exp比值及相关注释，输出窗口级结果。
- per_base: 基于窗口结果，计算每个位点的平均obs/exp置信区间及相关注释，输出位点级结果。
- annotate: 整合所有注释信息，输出每个位点的详细功能约束评分和注释。
命令行参数：
- -input: 注释文件路径，包含突变似然分数和观测最大杂合度。
- -obs: 观测最大杂合度的列名。
- -parameters: 线性模型参数文件路径，用于计算期望值。
- -prefix: 输出文件前缀。
- -exc_sites: 需排除的位点列表。
- -kmer_length: 分析窗口大小（bp）。
输出文件：
- output/local_constraint/[prefix]kmers_local_constraint.txt：窗口级局部约束结果。
- output/local_constraint/[prefix]per_base_local_constraint.txt：位点级局部约束结果。
依赖模块：
- oe_functions, annotate_mutations, compile_denovo
- pandas, numpy, csv, multiprocessing, argparse, datetime, os, sys
注意事项：
- 需提前准备MITOMAP功能区注释文件和相关功能注释文件。
- 计算过程涉及多进程并行，运行时间与计算资源相关。
'''
# 添加模块搜索路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '2-计算观测值和预期值'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '1-构建模型'))

# 现在可以正常导入模块
import oe_functions
import annotate_mutations
import compile_denovo


def lookup_dict(input_file: str, obs_value: str):
	"""创建一个字典，用于查找观测到的最大杂合度和似然值。

	:param input_file: 包含突变似然评分和观测最大杂合度的注释文件
	:param obs_value: 观测值的列标题
	:return: 字典，键为（位置、突变）元组，值为包含变异所在区域、观测值和似然值的元组
	"""
	dict = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		# 仅保留碱基和氨基酸替换，即以错义为最严重影响、RNA及非编码突变
		if (any(x in row["consequence"] for x in ["missense", "non_coding_transcript", "intergenic"])) and not any(
				x in row["consequence"] for x in oe_functions.more_severe_than_missense):
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in oe_functions.ori_region) else 'ref_exc_ori'
			observed = float(row[obs_value])
			likelihood = float(row["Likelihood"])
			callable_samples = float(row["callable_samples"])
			dict[(row["POS"], mutation)] = (region_to_use, observed, likelihood, callable_samples)
	return dict


def parallelize_kmers(
	position: int, 
	kmer_length: int,
	lookup_dictionary: Dict[Tuple[str, str], Tuple[str, float, float, float]],
	excluded_sites: List[int],
	fit_parameters: str,
	loci_dict: Dict[Tuple[int, int], Tuple[str, str]],
	kmer_dict: Dict[Tuple[int, int], Tuple[int, float, float, float, float, float, str, str]]):
	"""并行计算每个kmer的obs:exp比值及其90%置信区间。

	:param position: kmer的起始位置
	:param kmer_length: 用户指定的kmer长度
	:param lookup_dictionary: 字典，键为（位置、突变）元组，值为包含区域、观测值、似然值和可判读样本数的元组
	:param excluded_sites: 需要在计算中排除的碱基位置列表
	:param fit_parameters: 包含线性方程系数和截距的文件路径
	:param loci_dict: 基于MITOMAP注释的位点详细注释字典
	:param kmer_dict: 保存每个kmer相关注释信息的字典
	:return: 更新后的kmer字典，包含所有相关注释信息
	"""
	kmer_range = list(range(position, position + kmer_length))
	start = int(kmer_range[0])
	end = int(kmer_range[-1])
	if end > 16569:  # 处理环状基因组
		kmer_range = list(range(start, 16570)) + list(range(1, end - 16569 + 1))
		end = int(kmer_range[-1])

	# 初始化字典——每个kmer独立初始化，后续不会追加
	kmer_sum = oe_functions.initialize_sum_dict(identifier_list=[(start, end)])
	
	# 使用查找字典汇总该kmer的观测值和似然值
	for pos in kmer_range:
		if pos not in excluded_sites:
			for key in lookup_dictionary:
				if key[0] == str(pos):
					kmer_sum = oe_functions.sum_obs_likelihood(
						mutation=key[1], identifier=(start, end), region=lookup_dictionary[key][0], dict=kmer_sum,
						observed=lookup_dictionary[key][1], likelihood=lookup_dictionary[key][2],
						callable_samples=lookup_dictionary[key][3])

	# 计算观测值、期望值及置信区间
	(total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI) = oe_functions.calculate_oe(
		item=(start, end), sum_dict=kmer_sum, fit_parameters=fit_parameters, output="dict", file=None, alpha=0.1)

	# 注释重叠的基因区
	loci, loci_type = [], []  # 重叠位点列表
	for key in loci_dict:
		# 重叠条件，对于跨越人工切点的片段需要单独处理
		if start < end:  # 起始小于终止，说明未跨越人工切点
			if (start <= key[1]) and (end >= key[0]) and (loci_dict[key][0] not in loci):
				loci.append(loci_dict[key][0])
				loci_type.append(loci_dict[key][1])
		elif end < start:  # kmer跨越原点
			if ((end >= key[0]) and (end <= key[1])) or ((start >= key[0]) and (start <= key[1])) and (
					loci_dict[key][0] not in loci):
				loci.append(loci_dict[key][0])
				loci_type.append(loci_dict[key][1])
	loci = str(loci).strip('[]').replace("'", "")
	loci_type = str(loci_type).strip('[]').replace("'", "")
	
	# 将结果写入字典
	kmer_dict[(start, end)] = (total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI, loci, loci_type)
	return kmer_dict


def kmers(
	kmer_length: int,
	lookup_dictionary: Dict[Tuple[str, str], Tuple[str, float, float, float]],
	excluded_sites: List[int],
	fit_parameters: str,
	out_prefix: str):
	"""计算所有kmer的obs:exp比值并添加相关注释。

	:param kmer_length: 用户指定的每个kmer长度
	:param lookup_dictionary: 字典，键为（位置、突变）元组，值为包含区域、观测值和似然值的元组
	:param excluded_sites: 需要排除的碱基位置列表
	:param fit_parameters: 线性方程系数和截距文件的路径
	:param out_prefix: 输出文件名前缀
	:return: df，为包含全部kmer及注释细节的pandas数据框
	:return: coding_intervals，为基于MITOMAP的位点注释字典
	"""
	# 按坐标生成MITOMAP位点字典
	mitomap_loci = {}
	for row in csv.DictReader(open('0-required_files/databases/mitomap_genome_loci.txt'), delimiter='\t'):
		if int(row["Starting"]) > int(row["Ending"]):
			mitomap_loci[(1, int(row["Ending"]))] = (row["Map_Locus"], row["Description"])
			mitomap_loci[(int(row["Starting"]), 16569)] = (row["Map_Locus"], row["Description"])
		else:
			mitomap_loci[(int(row["Starting"]), int(row["Ending"]))] = (row["Map_Locus"], row["Description"])
	
	# 初始化在并行函数中使用的字典
	manager = mp.Manager()
	kmers = manager.dict()
	pool = mp.Pool(mp.cpu_count())
	
	# 计算所有可能kmer的oe比值与置信区间
	for position in list(range(1, 16570)):
		pool.apply_async(parallelize_kmers, args=(
			position, kmer_length, lookup_dictionary, excluded_sites, fit_parameters, mitomap_loci, kmers))
	pool.close()  # 关闭进程池，让所有任务执行完毕
	pool.join()  # 阻塞后续代码，直到队列中的所有进程完成
	
	# 使用pandas依据上限置信区间为每个kmer计算百分位排名
	df = pd.DataFrame.from_dict(kmers, orient='index', columns=[
		'total_all', 'obs_max_het', 'exp_max_het', 'ratio_oe', 'lower_CI', 'upper_CI', 'loci', 'loci_type'])
	df.index = pd.MultiIndex.from_tuples(df.index, names=['start', 'end'])
	df = df.reset_index()
	print(df)
	df['pctRank_OEUF'] = df['upper_CI'].rank(pct=True, ascending=False)  # 按obs:exp上限CI的百分位排名
	df = df.sort_values(['pctRank_OEUF'], ascending=False)  # 从小到大排序
	
	# 同时对蛋白基因、RNA基因和非编码区域分别计算排名
	df['loci_type'] = df['loci_type'].map(str)
	df['pctRank_OEUF_RNA'] = df.upper_CI[
		df['loci_type'].str.contains('RNA')].rank(pct=True, ascending=False)  # 若包含任何RNA基因
	df['pctRank_OEUF_protein'] = df.upper_CI[
		df['loci_type'].str.contains('subunit|Cytochrome')].rank(pct=True, ascending=False)  # 若包含任何蛋白编码区
	# 仅限非编码区域，排除含RNA或蛋白序列的片段
	df['pctRank_OEUF_noncoding'] = df.upper_CI[
		~df['loci_type'].str.contains('RNA|subunit|Cytochrome')].rank(pct=True, ascending=False)

	# 添加注释，首先创建所需字典
	phylop_dict = annotate_mutations.phylop_annotate()
	uniprot_dict = annotate_mutations.uniprot_annotations()
	CI_dict = annotate_mutations.curated_func_sites()
	mod_dict = annotate_mutations.RNA_domains_mods()
	vep_dict = annotate_mutations.vep_annotate()
	prot_pos_dict = {}
	
	# 创建DNA位置到蛋白位置的映射字典
	for key in vep_dict:
		if key[3] == "codon":
			prot_pos_dict[int(key[1])] = str(vep_dict[key]).strip('[]').replace("'", "").replace(" ", "")
	# 创建蛋白注释字典
	sites_dict = {}
	for key in uniprot_dict:
		if "site" in str(uniprot_dict[key]):
			sites_dict[key] = [str(uniprot_dict[key]).strip('[]').replace("'", "").replace(" ", "")]
	# 合并蛋白注释
	for key in CI_dict:
		if key in sites_dict:
			sites_dict[key].append(CI_dict[key])
		else:
			sites_dict[key] = [CI_dict[key]]
	
	# 初始化用于新列的列表
	prot_start, prot_end, phylop_cons, prot_sites, rna_modified = [], [], [], [], []
	
	# 遍历每一个kmer
	for index in df.index:
		start = df['start'][index].astype(int)
		end = df['end'][index].astype(int)
		lookup_range = list(range(start, end + 1))
		if start > (16570 - kmer_length):  # 处理环状基因组
			lookup_range = list(range(start, 16570)) + list(range(1, end + 1))
		
		# 添加蛋白位置范围——MT-ND6为反向链，需要调换顺序
		if 'MT-ND6' in df['loci'][index]:
			prot_start.append(prot_pos_dict[end])
			prot_end.append(prot_pos_dict[start])
		else:
			prot_start.append(prot_pos_dict[start])
			prot_end.append(prot_pos_dict[end])
			
		# 创建临时列表以收集该范围内的多条注释
		p_cons, p_sites, r_mod = [], [], []
		# 在整个kmer范围内汇总信息
		for i in lookup_range:
			if i != 3107:  # 该位点的参考碱基为N
				# 收集保守性指标
				p_cons.append(float(phylop_dict[i]))
				# 收集功能位点与修饰
				if (i in sites_dict) and (sites_dict[i] not in p_sites):
					p_sites.append(sites_dict[i])
				if ((str(i), "modified") in mod_dict) and (mod_dict[(str(i), "modified")] not in r_mod):
					r_mod.append(mod_dict[(str(i), "modified")])
				if (("\n" + str(i) + "\n") in open('0-required_files/other_annotations/rRNA_bridge_bases.txt').read()) and (
						"rRNA_bridge_base" not in r_mod):
					r_mod.append("rRNA_bridge_base")
		# 计算该kmer保守性指标的平均值
		phylop_cons.append(mean(p_cons))
		# 格式化功能位点与修饰信息
		prot_sites.append(str(p_sites).strip('[').strip(']').strip('\'').replace("', '", ",").replace(", ", "").strip())
		rna_modified.append(str(r_mod).strip('[').strip(']').strip('\'').replace("', '", ",").replace(", ", "").strip())
		
	# 将注释结果添加到数据框
	df['protein_position_start'] = prot_start
	df['protein_position_end'] = prot_end
	df['phylop_mean'] = phylop_cons
	df['protein_sites'] = prot_sites
	df['RNA_modifications'] = rna_modified
	print(df)

	# 写入结果文件
	np.savetxt(
		'output/local_constraint/%skmers_local_constraint.txt' % out_prefix, df, fmt='\t'.join(
			(['%i'] * 3) + (['%.6f'] * 5) + (['%s'] * 2) + (['%.6f'] * 4) + (['%s'] * 2) + (['%.6f'] * 1) + (['%s'] * 2)),
		header="start	end	variant_count	observed	expected	obs:exp	lower_CI	upper_CI	loci	loci_types	pctRank_OEUF	pctRank_OEUF_RNA	pctRank_OEUF_protein	pctRank_OEUF_noncoding	protein_pos_start	protein_pos_end	phylop_mean	protein_sites	RNA_sites")
	
	return df, mitomap_loci, sites_dict


def per_base(
	df: pd.DataFrame,
	kmer_length: int,
	mitomap_loci: Dict[Tuple[int, int], Tuple[str, str]]):
	"""计算mtDNA上每个位点所覆盖kmer的平均百分位排名以及其他注释。

	:param df: 包含全部kmer及其注释信息的pandas数据框
	:param kmer_length: 用户指定的kmer长度
	:param mitomap_loci: 基于MITOMAP注释的位点详细信息字典
	:return: df2，为包含mtDNA每个位点统计信息及注释的pandas数据框
	"""
	
	dict = {}
	for i in list(range(1, 16570)):  # mtDNA上每一个碱基
		# 筛选与位置i重叠的所有kmer
		if i > (16570 - kmer_length):  # 处理环状基因组
			filtered_df = df[
				(
						(df['start'] >= np.int64(i - kmer_length)) &
						(df['end'] <= np.int64(i + kmer_length - 16570))
				) | (
						(df['start'] <= np.int64(i)) &
						(df['end'] >= np.int64(i))
				)]
		elif i <= kmer_length:  # 处理环状基因组
			filtered_df = df[
				(
						(df['start'] >= np.int64(i - kmer_length + 16570)) &
						(df['end'] <= np.int64(i + kmer_length))
				) | (
						(df['start'] <= np.int64(i)) & (df['end'] >= np.int64(i))
				)]
		else:
			filtered_df = df[
				((df['start'] <= np.int64(i)) & (df['end'] >= np.int64(i))) |
				(df['start'] == np.int64(i)) |
				(df['end'] == np.int64(i))
			]
		# 对每个位置计算其重叠kmer的obs:exp上限CI平均值
		mean_OEUF = filtered_df['upper_CI'].mean()
		mean_exp = filtered_df['exp_max_het'].mean()

		# 注释区域，以便分别对蛋白、RNA基因及非编码区域排名
		loci, loci_type = [], []
		for key in mitomap_loci:
			if ((i >= key[0]) and (i <= key[1])) or (i == key[0]) or (i == key[1]):
				loci.append(mitomap_loci[key][0])
				loci_type.append(mitomap_loci[key][1])
		dict[i] = (mean_OEUF, loci, loci_type, mean_exp)

	df2 = pd.DataFrame(dict).T.reset_index()
	df2.columns = ['pos', 'mean_OEUF', 'loci', 'loci_type', 'mean_exp']
	df2['pctRank_mean_OEUF'] = df2['mean_OEUF'].rank(pct=True, ascending=False)
	df2 = df2.sort_values(['pctRank_mean_OEUF'])
	
	# 同时为RNA、蛋白和非编码区域分别计算排名
	df2['loci_type'] = df2['loci_type'].map(str)
	df2['pctRank_mean_OEUF_RNA'] = df2.mean_OEUF[
		df2['loci_type'].str.contains('RNA')].rank(pct=True, ascending=False)
	df2['pctRank_mean_OEUF_protein'] = df2.mean_OEUF[
		df2['loci_type'].str.contains('subunit|Cytochrome')].rank(pct=True, ascending=False)
	df2['pctRank_mean_OEUF_noncoding'] = df2.mean_OEUF[
		~df2['loci_type'].str.contains('RNA|subunit|Cytochrome')].rank(pct=True, ascending=False)
	
	return df2


def annotate(
	input_file: str,
	dataframe: pd.DataFrame,
	sites: Dict[int, List[str]],
	out_prefix: str):
	"""为mtDNA上每个位点提供逐位指标。

	:param input_file: 包含注释信息的文件路径，即annotate_output.py的输出
	:param dataframe: 包含每个位点统计信息和注释的pandas数据框
	:param sites: 键为位置、值为蛋白功能位点注释的字典
	:param out_prefix: 输出文件名前缀
	"""

	f = open('output/local_constraint/%sper_base_local_constraint.txt' % out_prefix, "w")
	header = "POS	REF	ALT	symbol	consequence	protein_position	MLC_pos_score	MLC_var_score	mean_OEUF	pctRank_mean_OEUF	pctRank_OEUF_protein	pctRank_OEUF_RNA	pctRank_OEUF_noncoding	mean_expected	in_phylotree	phyloP_score	protein_sites	RNA_modified	rRNA_bridge_base	apogee_class	mitotip_class	hmtvar_class	mitomap_status	mitomap_plasmy	mitomap_disease	clinvar_interp"
	f.write(header + '\n')
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mean_OEUF = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'mean_OEUF'].values).strip('[]')
		mean_exp = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'mean_exp'].values).strip('[]')
		pct_rank_OEUF = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF'].values).strip('[]')
		pct_rank_OEUF_RNA = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_RNA'].values).strip('[]')
		pct_rank_OEUF_protein = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_protein'].values).strip('[]')
		pct_rank_OEUF_noncoding = str(dataframe.loc[(df2['pos'] == np.int64(row["POS"])), 'pctRank_mean_OEUF_noncoding'].values).strip('[]')
		MLC_pos_score = pct_rank_OEUF
		
		# 这些分数衡量局部功能约束，仅基于RNA/非编码与错义变异
		# 其他功能变异在该范围内按其obs:exp分布进行手动指定
		MLC_var_score = pct_rank_OEUF
		if ("stop_gain" in row["consequence"]) and ('stop_gained&start_lost' not in row["consequence"]):
			MLC_var_score = str(1)
			pct_rank_OEUF_protein = 'nan'
		if (("synonymous" in row["consequence"]) or ("stop_retained" in row["consequence"])) and not any(
				x in row["consequence"] for x in oe_functions.more_severe_than_syn):
			MLC_var_score = str(0)
			pct_rank_OEUF_protein = 'nan'
		# 起始丢失或终止丢失（数量非常少）
		elif any(x in row["consequence"] for x in ["stop_lost", "start_lost", "incomplete_terminal"]):
			MLC_var_score = str(0.70)  # 手动查找的等效值
			pct_rank_OEUF_protein = 'nan'
			
		prot_sites = str(sites[int(row["POS"])]).strip('[]').replace("'", "").replace(" ", "") \
			if int(row["POS"]) in sites else ''
		
		RNA_bridge = "Yes" if ("\n" + row["POS"] + "\n") in open(
			'0-required_files/other_annotations/rRNA_bridge_bases.txt').read() else "No"
		
		f.write(
			row["POS"] + '\t' + row["REF"] + '\t' + row["ALT"] + '\t' +
			row["symbol"] + '\t' + row["consequence"] + '\t' + row["protein_position"] + '\t' +
			str(MLC_pos_score) + '\t' + str(MLC_var_score) + '\t' +
			mean_OEUF + '\t' + pct_rank_OEUF + '\t' +
			pct_rank_OEUF_protein + '\t' + pct_rank_OEUF_RNA + '\t' + pct_rank_OEUF_noncoding + '\t' + mean_exp + '\t' +
			row["in_phylotree"] + '\t' + row["phyloP_score"] + '\t' +
			prot_sites + '\t' + row["RNA_modified"] + '\t' + RNA_bridge + '\t' +
			row["apogee_class"] + '\t' + row["mitotip_class"] + '\t' + row["hmtvar_class"] + '\t' +
			row["mitomap_status"] + '\t' + row["mitomap_plasmy"] + '\t' + row["mitomap_disease"] + '\t' +
			row["clinvar_interp"] + '\n')


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
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calculations")
	parser.add_argument(
		"-kmer_length", type=int, help="Window size to use for analysis in bp")
	args = parser.parse_args()

	# Set default parameters for carrier-count analysis
	if args.input is None:
		args.input = "output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt"
	if args.obs is None:
		args.obs = "carrier_count"
	if args.parameters is None:
		args.parameters = "output/calibration/linear_model_fits.txt"
	if args.prefix is None:
		args.prefix = ""
	if args.exc_sites is None:
		args.exc_sites = [301, 302, 310, 316, 16182]
	print(datetime.datetime.now(), "Analyze local constraint in the mtDNA!")
	print("This may take a little time depending on available computing power...")
	
	(df, mitomap_loci, sites_dict) = kmers(
		kmer_length=args.kmer_length, lookup_dictionary=lookup_dict(input_file=args.input, obs_value=args.obs),
		excluded_sites=args.exc_sites, fit_parameters=args.parameters, out_prefix=args.prefix)
	
	df2 = per_base(df=df, kmer_length=args.kmer_length, mitomap_loci=mitomap_loci)
	
	annotate(input_file=args.input, dataframe=df2, sites=sites_dict, out_prefix=args.prefix)

	print(datetime.datetime.now(), "Script finished!")
