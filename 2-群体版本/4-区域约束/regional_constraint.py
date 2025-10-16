import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import sys
from typing import Dict, List, Tuple, Union

# 添加模块搜索路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '2-计算观测值和预期值'))

# 现在可以正常导入模块
import oe_functions


def res_to_pos():
	"""创建一个字典，用于将蛋白氨基酸残基编号转换为对应的碱基位置列表。

	:return: res_to_pos字典，键为（基因座、残基编号）元组，值为碱基位置列表
	"""
	dict = {}
	for row in csv.DictReader(
			open("0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		if row["Protein_position"]:
			if (row["SYMBOL"], int(row["Protein_position"])) not in dict:
				dict[(row["SYMBOL"], int(row["Protein_position"]))] = [int(row["POS"])]
			else:
				dict[(row["SYMBOL"], int(row["Protein_position"]))].append(int(row["POS"]))
	return dict
	
	
def pos_to_res():
	"""创建一个字典，用于将碱基位置转换为对应的蛋白残基编号。

	:return: pos_to_res字典，键为（基因座、碱基位置）元组，值为残基编号
	"""
	dict = {}
	for row in csv.DictReader(
			open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		if row["Protein_position"] and ((row["SYMBOL"], int(row["POS"])) not in dict):
			dict[(row["SYMBOL"], int(row["POS"]))] = int(row["Protein_position"])
	return dict


def calculate_per_pos(
		locus: str,
		locus_coords: List[int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int]):
	"""汇总每个碱基位置的观测携带者计数、突变似然值及计数。

	:param locus: 基因座名称
	:param locus_coords: 构成该基因座的碱基坐标列表
	:param input_file: 包含突变似然评分和观测携带者计数的注释文件
	:param obs_value: 观测值列标题
	:param excluded_sites: 需要在计算中排除的碱基位置列表
	:return: per_pos字典，键为（碱基位置、突变类别、汇总值类型、变异所在区域）元组，值为对应的汇总结果
	"""
	# 处理同时属于两个基因的变异
	per_gene = {}
	for row in csv.DictReader(open(
			'0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		# 控制区的命名与输入保持一致
		row["SYMBOL"] = "CR" if ((int(row["POS"]) <= 576) or (int(row["POS"]) >= 16024)) else row["SYMBOL"]
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]
	
	# 将基因座坐标转换为字符串列表
	per_pos = oe_functions.initialize_sum_dict(identifier_list=[str(x) for x in locus_coords])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in oe_functions.ori_region) else 'ref_exc_ori'
			if row["POS"] in [str(x) for x in locus_coords]:
				# RNA基因和非编码基因座保留全部碱基；蛋白编码基因仅保留最严重为错义的变异
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
						x in row["consequence"] for x in oe_functions.more_severe_than_missense)):
					per_pos = oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=row["POS"], region=region_to_use,
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=per_pos)
	return per_pos


def calculate_per_res(
		locus: str,
		locus_coords: List[int],
		codon_coords: List[int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int],
		pos_to_res: Dict[Tuple[str, int], int]):
	"""汇总每个蛋白残基的观测携带者计数、突变似然值及计数。

	:param locus: 基因座名称
	:param locus_coords: 构成该基因座的碱基坐标列表
	:param codon_coords: 蛋白残基编号列表（从1开始）
	:param input_file: 包含突变似然评分和观测携带者计数的注释文件
	:param obs_value: 观测值列标题
	:param excluded_sites: 需要排除的碱基位置列表
	:param pos_to_res: 将碱基坐标转换为蛋白残基编号的字典
	:return: per_res字典，键为（残基、突变类别、汇总值类型、变异所在区域）元组，值为对应的汇总结果
	"""
	# 处理同时属于两个基因的变异
	per_gene = {}
	for row in csv.DictReader(open(
			'0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]
	
	# 将坐标转换为字符串列表
	per_res = oe_functions.initialize_sum_dict(identifier_list=[str(x) for x in codon_coords])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in oe_functions.ori_region) else 'ref_exc_ori'
			if row["POS"] in [str(x) for x in locus_coords]:
				# RNA基因和非编码基因座保留全部碱基；蛋白编码基因仅保留最严重为错义的变异
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
						x in row["consequence"] for x in oe_functions.more_severe_than_missense)):
					# 注意：需使用pos_to_res处理属于双基因的位点
					per_res = oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=str(pos_to_res[(locus, int(row["POS"]))]),
						region=region_to_use, observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=per_res)
	return per_res


def shuffle(
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		coords: List[str]):
	"""对真实比对进行随机置换，将观测值、突变似然值和变异计数随机分配至其他位置。

	:param per_unit_dict: 非蛋白基因使用per_pos，蛋白基因使用per_res
	:param coords: 需置换的碱基坐标或蛋白残基编号列表
	:return: 随机后字典，键为（位置、突变类别、汇总值类型、区域）元组，值为对应的汇总结果
	"""
	# 首先重构per_unit_dict便于随机化
	shuffle_dict = {}
	col_names = []  # 记录列名，保持追加顺序一致
	for key in per_unit_dict:
		pos = key[0]
		if pos not in shuffle_dict:
			shuffle_dict[pos] = [per_unit_dict[key]]
			if (str(key[1]) + '-' + str(key[2]) + '-' + str(key[3])) not in col_names:
				col_names.append(str(key[1]) + '-' + str(key[2]) + '-' + str(key[3]))
		else:
			shuffle_dict[pos].append(per_unit_dict[key])
			if (str(key[1]) + '-' + str(key[2]) + '-' + str(key[3])) not in col_names:
				col_names.append(str(key[1]) + '-' + str(key[2]) + '-' + str(key[3]))
	print("after first reformat", shuffle_dict, col_names)

	df = pd.DataFrame.from_dict(shuffle_dict, orient='index', columns=col_names)
	pd.set_option('display.max_columns', None)
	pd.set_option('display.max_rows', None)
	print("before shuffle", df)

	df = df.sample(frac=1).reset_index()  # 在不放回的情况下打乱所有行并重建索引
	print("after shuffle", df)

	df.set_axis(coords, axis='index', inplace=True)
	print("renumber index", df)

	df = df.drop(columns=['index'])  # 删除原始索引列
	shuffle_dict = df.T.to_dict(orient='list')  # 转置以位置为键重新生成字典
	print("new shuffled", shuffle_dict)
	pd.set_option('display.max_rows', 5)

	# 将结果转换回输入所需的字典格式
	per_unit_dict = {}
	# column保存键中各标识，index用于按顺序读取shuffle_dict中的值
	for index, column in enumerate(col_names):
		for key in shuffle_dict:  # key为位置
			if column.split('-')[1] == "count":  # 计数需转换为整数
				per_unit_dict[(
					str(key), column.split('-')[0], column.split('-')[1], column.split('-')[2])] = int(shuffle_dict[key][index])
			else:  # 其余保持浮点类型
				per_unit_dict[(
					str(key), column.split('-')[0], column.split('-')[1], column.split('-')[2])] = float(shuffle_dict[key][index])
	print("after second reformat", per_unit_dict)

	return per_unit_dict

	
def calculate_kmers_step1(
		pos: int,
		kmer_length: int,
		coords: List[int],
		locus_oe: float,
		locus_lower_CI: float,
		fit_parameters: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		outlier_kmers: Dict[Tuple[int, int], Tuple[float, int]],
		sig_threshold: float):
	"""计算所有从基因座位置1开始的kmer的obs:exp值及其p值。

	:param pos: kmer的起始位置（碱基或残基）
	:param kmer_length: 需要评估的kmer长度
	:param coords: 基因座中所有坐标列表，蛋白使用codon_coords，非蛋白使用locus_coords
	:param locus_oe: 该基因座的观测/期望比
	:param locus_lower_CI: 该基因座obs:exp比值的90%置信区间下界
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param per_unit_dict: 字典，键为（位置或残基、突变类别、汇总值类型、区域）元组，值为对应的汇总结果
	:param kmers_dict: 字典，键为（kmer起止位置、突变类别、汇总值类型、区域）元组，值为对应的汇总结果
	:param outlier_kmers: 显著kmer字典，键为（起始、终止）坐标，值为（p值、长度）
	:param sig_threshold: p值阈值，仅保留p小于该值的kmer
	:return: 更新后的kmers_dict与outlier_kmers
	"""
	kmer_start = pos
	kmer_end = pos + kmer_length - 1  # 减1保证覆盖至基因终点
	
	if kmer_end > 16569:  # 处理环状基因组
		kmer_range = list(range(kmer_start, 16570)) + list(range(1, (kmer_end - 16569 + 1)))
		kmer_end = int(kmer_range[-1])
	else:
		kmer_range = list(range(kmer_start, kmer_end + 1))  # kmer覆盖的全部位置

	if kmer_end in coords:  # 限制在基因座范围内
		# 初始化字典条目
		for mut_group in ['G>A_and_T>C', 'other']:  # 两类突变群组
			for value in ['observed_carriers', 'sum_LR_AN', 'sum_callable', 'count']:  # 需要汇总的指标
				for region in ['ref_exc_ori', 'ori']:  # 可能使用的两个突变模型
					kmers_dict[(kmer_start, kmer_end), mut_group, value, region] = 0
		for pos in kmer_range:
			for key in per_unit_dict:
				if int(key[0]) == pos:
					kmers_dict[((kmer_start, kmer_end), key[1], key[2], key[3])] += per_unit_dict[key]
		
		# 计算p值前，先得到观测值与期望值
		expected_carriers = oe_functions.calculate_exp(
			sum_dict=kmers_dict, identifier=(kmer_start, kmer_end), fit_parameters=fit_parameters)
		observed_carriers = oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		observed_carriers = 0.0 if observed_carriers < 0.001 else observed_carriers  # 处理极小数值
		total_all = oe_functions.calculate_total(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		# 期望值基于基因座的obs:exp比进行缩放
		less_than_pvalue = oe_functions.calculate_pvalue(
			observed_carriers=observed_carriers, total=total_all, expected_carriers=(expected_carriers * locus_oe))
		
		# print((kmer_start, kmer_end), observed_carriers, expected_carriers, total_all, less_than_pvalue)
		
		if float(less_than_pvalue) < sig_threshold:  # 仅保留p值低于阈值的kmer
			(lower_CI, upper_CI) = oe_functions.calculate_CI(
				observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers, alpha=0.1)
			if upper_CI < locus_lower_CI:  # 置信区间需与基因座CI不重叠
				outlier_kmers[(kmer_start, kmer_end)] = (less_than_pvalue, kmer_length)

	return kmers_dict, outlier_kmers


def calculate_kmers_step2(
		pos: int,
		kmer_length: int,
		locus: str,
		coords: List[int],
		fit_parameters: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		outlier_kmers: Dict[Tuple[int, int], Tuple[float, int]],
		locus_oe: float,
		locus_lower_CI: float,
		sig_threshold: float):
	"""计算所有不从基因座位置1开始的kmer的obs:exp值及p值。

	:param pos: kmer的起始位置（碱基或残基）
	:param kmer_length: 需要评估的kmer长度
	:param locus: 基因座名称
	:param coords: 基因座中的坐标列表，蛋白使用codon_coords，非蛋白使用locus_coords
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param per_unit_dict: 字典，键为（位置或残基、突变类别、汇总值类型、区域）元组，值为对应的汇总结果
	:param kmers_dict: 字典，键为（kmer起止位置、突变类别、汇总值类型、区域）元组，值为对应的汇总结果
	:param outlier_kmers: 显著kmer字典，键为（起始、终止）坐标，值为（p值、长度）
	:param locus_oe: 基因座的观测/期望比
	:param locus_lower_CI: 基因座obs:exp比的90%置信区间下界
	:param sig_threshold: p值阈值，仅保留p小于该值的kmer
	:return: 更新后的kmers_dict与outlier_kmers
	"""
	kmer_start = pos
	kmer_end = pos + kmer_length - 1
	# 利用step1得到的kmers_dict快速计算
	# 例如：起始于2、终止于10的kmer可由起始于1、终止于10的kmer减去位置1的数值得到
	kmer_start_minus_one = 16569 if (locus == "CR" and pos == 1) else (pos - 1)
	
	if kmer_end > 16569:  # 处理环状基因组
		k_range = list(range(kmer_start, 16570)) + list(range(1, (kmer_end - 16569 + 1)))
		kmer_end = int(k_range[-1])
		
	if kmer_end in coords:  # 限制在基因座范围内
		# 初始化字典条目并记录键列表
		keys = []
		for identifier in [(kmer_start, kmer_end)]:
			for mut_group in ['G>A_and_T>C', 'other']:  # 两类突变群组
				for value in ['observed_carriers', 'sum_LR_AN', 'sum_callable', 'count']:  # 汇总指标
					for region in ['ref_exc_ori', 'ori']:  # 两种突变模型
						kmers_dict[identifier, mut_group, value, region] = 0
						keys.append(mut_group + '-' + value + '-' + region)
		
		# 通过“更长kmer减去单个位点”的方式得到当前kmer的数值
		for key in keys:
			kmers_dict[((kmer_start, kmer_end), key.split('-')[0], key.split('-')[1], key.split('-')[2])] += \
				(kmers_dict[(
					(kmer_start_minus_one, kmer_end), key.split('-')[0], key.split('-')[1], key.split('-')[2])] - per_unit_dict[
					(str(kmer_start_minus_one), key.split('-')[0], key.split('-')[1], key.split('-')[2])])
		
		# 计算p值前，先得到观测值与期望值
		expected_carriers = oe_functions.calculate_exp(
			sum_dict=kmers_dict, identifier=(kmer_start, kmer_end), fit_parameters=fit_parameters)
		observed_carriers = oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		observed_carriers = 0.0 if observed_carriers < 0.001 else observed_carriers  # 处理极小数值
		total_all = oe_functions.calculate_total(sum_dict=kmers_dict, identifier=(kmer_start, kmer_end))
		# 期望值基于基因座的obs:exp比进行缩放
		less_than_pvalue = oe_functions.calculate_pvalue(
			observed_carriers=observed_carriers, total=total_all, expected_carriers=(expected_carriers * locus_oe))
		
		# print(locus, (kmer_start, kmer_end), observed_carriers, expected_carriers, total_all, less_than_pvalue)
		
		if float(less_than_pvalue) < sig_threshold:  # 仅保留p值低于阈值的kmer
			(lower_CI, upper_CI) = oe_functions.calculate_CI(
				observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers, alpha=0.1)
			if upper_CI < locus_lower_CI:  # 要求置信区间与基因座CI不重叠
				outlier_kmers[(kmer_start, kmer_end)] = (less_than_pvalue, kmer_length)

	return kmers_dict, outlier_kmers


def extract_nonoverlapping_kmers(
		key: Tuple[int, int],
		df: pd.DataFrame,
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		locus: str,
		fit_parameters: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		filtered_outlier_kmers: Dict[Tuple[int, int], Tuple[float, float, float, float, float, int, float, int]]):
	"""根据最低p值（必要时再比较长度）提取互不重叠的kmer集合。

	:param key: kmer起止坐标元组；处理控制区时，将1-576的坐标平移到16569之后
	:param df: 包含全部显著kmer的pandas数据框
	:param kmers_dict: 字典，键为（kmer起止位置、突变类别、汇总值类型、区域）元组，值为汇总结果
	:param locus: 基因座名称
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param res_to_pos: 将蛋白残基转换为碱基坐标的字典
	:param filtered_outlier_kmers: 用于存放筛选后kmer的字典
	:return: 更新后的filtered_outlier_kmers
	"""
	kmer_start = key[0]
	kmer_end = key[1]
	# 筛选所有与当前kmer重叠的显著片段，并按p值从小到大（同p值则长度更长优先）排序
	# 过滤条件确保两个区间确实重叠
	pre_filter_df = df[
		(df['kmer_start'] <= np.int64(kmer_end)) &
		(df['kmer_end'] >= np.int64(kmer_start))].sort_values(['less_than_pvalue', 'length'], ascending=(True, False))

	# 提取当前kmer的p值和长度
	kmer_pvalue = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
				pre_filter_df['kmer_end'] == np.int64(kmer_end))]['less_than_pvalue'].iloc[0]
	kmer_length = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
				pre_filter_df['kmer_end'] == np.int64(kmer_end))]['length'].iloc[0]
	# 将当前kmer从数据框中移除，以便比较其他候选
	filter_df = pre_filter_df.drop(pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))].index)
	
	# 依据以下三种情况选出互不重叠的kmer
	# (1) 若不存在任何重叠kmer
	if filter_df.shape[0] == 0:
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 1", key, "to final dict")
		
	# (2) 当前kmer的p值低于所有重叠kmer
	elif kmer_pvalue < filter_df['less_than_pvalue'].iloc[0]:
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 2", key, "to final dict")
	
	# (3) 当前kmer的p值等于最低p值，但长度更长或相同
	elif (kmer_pvalue == filter_df['less_than_pvalue'].iloc[0]) and (kmer_length >= filter_df['length'].iloc[0]):
		filtered_outlier_kmers = annotate_nonoverlapping(
			key=key, kmers_dict=kmers_dict, fit_parameters=fit_parameters, locus=locus, res_to_pos=res_to_pos,
			kmer_pvalue=kmer_pvalue, kmer_length=kmer_length, filtered_outlier_kmers=filtered_outlier_kmers)
		print("added kmer 3", key, "to final dict")

	return filtered_outlier_kmers


def annotate_nonoverlapping(
		key: Tuple[int, int],
		kmers_dict: Dict[Tuple[Tuple[int, int], str, str, str], Union[float, int]],
		fit_parameters: str,
		locus: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		kmer_pvalue: float,
		kmer_length: int,
		filtered_outlier_kmers: Dict[Tuple[int, int], Tuple[float, float, float, float, float, int, float, int]]):
	"""为互不重叠的kmer添加注释，以便写入文件。
	
	:param key: kmer起止坐标元组；控制区需将1-576位置平移至16569之后
	:param kmers_dict: 字典，键为（kmer起止位置、突变类别、汇总值类型、区域）元组，值为汇总结果
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param locus: 基因座名称
	:param res_to_pos: 将蛋白残基转换为碱基坐标的字典
	:param kmer_pvalue: kmer对应的p值
	:param kmer_length: kmer长度；蛋白使用残基数，非蛋白使用碱基数
	:param filtered_outlier_kmers: 存放筛选后kmer的字典
	:return: 更新后的filtered_outlier_kmers
	"""
	# 处理控制区（将1-576位置平移到16569之后）
	key_start = (key[0] - 16569) if (key[0] > 16569) else key[0]
	key_end = (key[1] - 16569) if (key[1] > 16569) else key[1]
	key = (key_start, key_end)
	
	# 计算obs:exp比值及置信区间
	observed_carriers = oe_functions.calculate_obs(sum_dict=kmers_dict, identifier=key)
	observed_carriers = 0 if observed_carriers < 0.001 else observed_carriers  # 处理极小数值
	expected_carriers = oe_functions.calculate_exp(
		sum_dict=kmers_dict, identifier=key, fit_parameters=fit_parameters)
	ratio_oe = observed_carriers / expected_carriers
	total_all = oe_functions.calculate_total(sum_dict=kmers_dict, identifier=key)
	(lower_CI, upper_CI) = oe_functions.calculate_CI(
		observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers, alpha=0.1)
	
	# 将残基范围转换为mtDNA坐标：通常对应第一密码子首位到最后密码子末位
	# 对于反链上的MT-ND6，需取最后密码子末位到第一密码子首位
	if any(x in locus for x in protein_genes):
		key = (min(res_to_pos[(locus, key[1])]), max(res_to_pos[(locus, key[0])])) if (locus == "MT-ND6") \
			else (min(res_to_pos[(locus, key[0])]), max(res_to_pos[(locus, key[1])]))
		kmer_length = kmer_length * 3  # 转换为碱基长度
	
	# 以碱基坐标保存结果
	filtered_outlier_kmers[key] = (
		observed_carriers, expected_carriers, ratio_oe, lower_CI, upper_CI, total_all, kmer_pvalue, kmer_length)
	
	return filtered_outlier_kmers


def check_locus_oe(
		input: str,
		input_noncoding: str,
		locus: str,
		per_unit_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		fit_parameters: str):
	"""计算指定基因座的obs:exp比值，同时用于核对代码结果。

	:param input: 含蛋白基因座obs:exp值的文件路径
	:param input_noncoding: 含非编码基因座obs:exp值的文件路径
	:param locus: 基因座名称
	:param per_unit_dict: 蛋白使用per_res，非蛋白使用per_pos
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:return: 基因座的obs:exp比值及其90%置信区间下界
	"""
	locus_oe = ''
	for file in (input, input_noncoding):
		for row in csv.DictReader(open(file), delimiter='\t'):
			row["locus"] = "CR" if (row["locus"] == "MT-CR") else row["locus"]  # 方便后续解析
			if (row["consequence"] == "SNV" or row["consequence"] == "missense") and (row["locus"] == locus):
				locus_oe = float(row["obs:exp"])

	# 使用per_unit_dict重新计算基因座的obs:exp，用于校验
	loci_dict = oe_functions.initialize_sum_dict(identifier_list=[locus])
	for key in per_unit_dict:
		loci_dict[(locus, key[1], key[2], key[3])] += per_unit_dict[key]

	# 计算观测/期望比
	observed_carriers = oe_functions.calculate_obs(sum_dict=loci_dict, identifier=locus)
	expected_carriers = oe_functions.calculate_exp(
		sum_dict=loci_dict, identifier=locus, fit_parameters=fit_parameters)
	total_all = oe_functions.calculate_total(sum_dict=loci_dict, identifier=locus)
	(lower_CI, upper_CI) = oe_functions.calculate_CI(
		observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers, alpha=0.1)
	
	if round(locus_oe, 2) != round(observed_carriers / expected_carriers, 2):  # 四舍五入以忽略极小差异
		print("Error: gene obs:exp does not equal expected value", locus_oe, "vs", observed_carriers / expected_carriers)
		sys.exit()
	
	return locus_oe, lower_CI


def regional_constraint(
		out_dir: str,
		n_shuffles: int,
		loci_type: str,
		res_to_pos: Dict[Tuple[str, int], List[int]],
		pos_to_res: Dict[Tuple[str, int], int],
		input_file: str,
		obs_value: str,
		excluded_sites: List[int],
		fit_parameters: str,
		loci_input: str,
		noncoding_input: str,
		sig_threshold: float):
	"""识别基因中显著强于整体约束的、互不重叠的区域集合。

	:param out_dir: 输出目录名称
	:param n_shuffles: 随机置换次数，为0时表示使用真实比对结果
	:param loci_type: 'p'表示仅分析蛋白基因，'rn'表示RNA与非编码，'both'表示全部
	:param res_to_pos: 将蛋白残基映射到碱基坐标的字典
	:param pos_to_res: 将碱基坐标映射到蛋白残基的字典
	:param input_file: 包含突变似然评分和观测携带者计数的注释文件路径
	:param obs_value: 观测值列标题
	:param excluded_sites: 需要在计算中排除的碱基位置列表
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param loci_input: 含基因座obs:exp值的文件路径
	:param noncoding_input: 含非编码obs:exp值的文件路径
	:param sig_threshold: p值阈值，仅保留低于该阈值的kmer
	"""
	# 准备输出文件
	out_dir = "real_alignment" if n_shuffles == 0 else out_dir
	out_dir = 'output/regional_constraint/' + out_dir
	kmer_file = open('%s/all_significant_kmers.txt' % out_dir, "w")
	kmer_file.write("kmer_start	kmer_end	locus	pvalue	length	shuffle" + '\n')

	kmer_file2 = open('%s/regional_constraint_intervals.txt' % out_dir, "w")
	header = "kmer_start	kmer_end	locus	protein_pos_start	protein_pos_end	observed_carriers	expected_carriers	ratio_oe	lower_CI	upper_CI	total	kmer_pvalue	length	shuffle	out_dir"
	kmer_file2.write(header + '\n')

	# n_shuffles = 0 表示使用真实比对
	number_shuffles = [0] if (int(n_shuffles) == 0) else list(range(1, n_shuffles + 1))
	
	for shuffle_num in number_shuffles:
		for row in csv.DictReader(open('0-required_files/databases/mitomap_genome_loci.txt'), delimiter='\t'):
			row["Map_Locus"] = row["Map_Locus"] if (row["Map_Locus"] != "MT-CR") else "CR"  # 方便后续解析
			if any(x in row["Map_Locus"] for x in loci_list):
				if (loci_type == "p") and not any(x in row["Map_Locus"] for x in protein_genes):
					continue  # 跳过当前基因座
				elif (loci_type == "rn") and any(x in row["Map_Locus"] for x in protein_genes):
					continue  # 跳过当前基因座
	
				# 设置参数并准备输入
				print("For locus", row["Map_Locus"])
				locus_start = int(row["Starting"])
				locus_end = int(row["Ending"])
				locus = row["Map_Locus"]
				locus_coords = list(range(locus_start, locus_end + 1)) if (locus_start < locus_end) else \
					list(range(locus_start, 16570)) + list(range(1, locus_end + 1))  # 处理跨越m.16569-1的情况
				max_kmer_length = len(locus_coords)
				# 蛋白基因需生成残基编号列表，并更新最大kmer长度
				if any(x in locus for x in protein_genes):
					codon_coords = list(range(1, pos_to_res[(locus, locus_start)] + 1)) if (locus == "MT-ND6") else \
						list(range(1, pos_to_res[(locus, locus_end)] + 1))
					max_kmer_length = len(codon_coords)
				else:
					codon_coords = None
				# 与局部约束一致，蛋白kmer的最小窗口较短
				# 构建所有待评估的kmer长度列表，最大值为基因座长度
				kmer_lengths = list(range(10, max_kmer_length + 1)) if any(x in locus for x in protein_genes) else \
					list(range(20, max_kmer_length + 1))
				
				# 构建per_pos字典
				per_pos = calculate_per_pos(
					locus=locus, locus_coords=locus_coords, input_file=input_file,
					obs_value=obs_value, excluded_sites=excluded_sites)
				# 蛋白基因需额外计算每个残基的汇总量
				per_res = calculate_per_res(
					locus=locus, locus_coords=locus_coords, codon_coords=codon_coords,
					input_file=input_file, obs_value=obs_value, excluded_sites=excluded_sites,
					pos_to_res=pos_to_res) if any(x in locus for x in protein_genes) else None
				# per_unit_dict：蛋白时使用per_res，非蛋白使用per_pos
				per_unit_dict = per_res if any(x in locus for x in protein_genes) else per_pos
				coords = codon_coords if any(x in locus for x in protein_genes) else locus_coords
				
				# 如需随机置换，则构建置换序列
				if shuffle_num > 0:
					per_unit_dict = shuffle(per_unit_dict=per_unit_dict, coords=coords)
					print(datetime.datetime.now(), "Shuffle number", shuffle_num, "for locus", locus, "has completed!")
				
				# 计算基因座obs:exp供区域约束分析使用，同时用于校验
				(locus_oe, locus_lower_CI) = check_locus_oe(
					input=loci_input, input_noncoding=noncoding_input, locus=locus, per_unit_dict=per_unit_dict,
					fit_parameters=fit_parameters)
	
				# 区域约束分析，使用并行处理
				manager = mp.Manager()
				kmers_dict = manager.dict()  # 每次迭代/置换重新初始化
				outlier_kmers = manager.dict()
	
				# 遍历基因座的每个位置并行计算所有可能kmer的obs:exp与p值
				# 函数返回更新后的kmers_dict和outlier_kmers
				for pos in coords:  # 基因座中的每个碱基或残基
					pool = mp.Pool(mp.cpu_count())
					if pos == coords[0]:  # 基因座的首个位置
						for kmer_length in kmer_lengths:  # 当前kmer长度
							pool.apply_async(calculate_kmers_step1, args=(
								pos, kmer_length, coords, locus_oe, locus_lower_CI, fit_parameters, per_unit_dict, kmers_dict, outlier_kmers, sig_threshold))
					else:  # 非首位
						for kmer_length in kmer_lengths:
							pool.apply_async(calculate_kmers_step2, args=(
								pos, kmer_length, locus, coords, fit_parameters, per_unit_dict, kmers_dict, outlier_kmers, locus_oe, locus_lower_CI, sig_threshold))
					pool.close()
					pool.join()  # 阻塞后续操作，直至队列中的进程完成
	
				print(datetime.datetime.now(), "Shuffle number", shuffle_num, "for locus", locus, "has all kmers calculated!")
	
				# 输出所有显著kmer
				for key in sorted(outlier_kmers):  # 并行可能导致顺序错乱，需排序
					kmer_file.write(
							str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' + str(outlier_kmers[key][0])
							+ '\t' + str(outlier_kmers[key][1]) + '\t' + str(shuffle_num) + '\n')
					kmer_file.flush()  # 确保及时写入文件
		
				# 进一步筛选：寻找不存在更低p值重叠区段的显著kmer
				if len(outlier_kmers) > 0:
					# 为方便处理控制区，将1-576位置平移到16569之后
					filtered_outlier_kmers = manager.dict()
					new_dict = {}
					for key in sorted(outlier_kmers):
						key_start = (16569 + key[0]) if (locus == "CR" and key[0] < 577) else key[0]
						key_end = (16569 + key[1]) if (locus == "CR" and key[1] < 577) else key[1]
						new_dict[(key_start, key_end)] = (outlier_kmers[key][0], outlier_kmers[key][1])
					# 基于转换后的坐标构建数据框
					df = pd.DataFrame.from_dict(new_dict, orient='index', columns=['less_than_pvalue', 'length'])
					df.index = pd.MultiIndex.from_tuples(df.index, names=['kmer_start', 'kmer_end'])
					df = df.reset_index()
					pool = mp.Pool(mp.cpu_count())
					for key in list(sorted(new_dict)):  # 遍历每个显著kmer
						pool.apply_async(extract_nonoverlapping_kmers, args=(
							key, df, kmers_dict, locus, fit_parameters, res_to_pos, filtered_outlier_kmers))
					pool.close()
					pool.join()
					
					for key in sorted(filtered_outlier_kmers):  # 并行可能改变顺序，需排序输出
						kmer_file2.write(
							str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
							str(pos_to_res[(locus, key[0])] if (locus, key[0]) in pos_to_res else '') + '\t' +
							str(pos_to_res[(locus, key[1])] if (locus, key[1]) in pos_to_res else '') + '\t' +
							'\t'.join(str(x) for x in filtered_outlier_kmers[key]) + '\t' +
							str(shuffle_num) + '\t' + str(out_dir) + '\n')
						kmer_file2.flush()  # 确保及时写入文件
					
					print("For loci", row["Map_Locus"], "the non-overlapping regions have been extracted!")
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed carrier counts")
	parser.add_argument(
		"-obs", type=str, help="Population dataset providing observed carrier counts")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calculations")
	parser.add_argument(
		"-n_shuffles", type=int, help="Number of shuffled alignments to generate, optional")
	parser.add_argument(
		"-loci_type", type=str, help="Options are 'p' for protein, 'rn' for RNA and non-coding, or 'both' for both")
	parser.add_argument(
		"-loci_input", type=str, help="Path to file with locus obs:exp values")
	parser.add_argument(
		"-noncoding_input", type=str, help="Path to file with noncoding loci obs:exp values")
	parser.add_argument(
		"-out_dir", type=str, help="Name of the output directory")
	parser.add_argument(
		"-sig_threshold", type=float, help="P-value threshold to apply")
	args = parser.parse_args()
	
	# 设置默认参数，适用于gnomAD
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "carrier_count"
	if args.out_dir is None:
		args.out_dir = "real_alignment"
	if args.parameters is None:
		args.parameters = 'output/calibration/linear_model_fits.txt'
	if args.exc_sites is None:
		# 排除gnomAD中易受技术伪影影响的位置：301、302、310、316、3107和16182（3107已在前面排除）
		# 这些位点未在gnomAD中检测到，因此需在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]
	if args.n_shuffles is None:
		args.n_shuffles = 0  # 使用真实比对
	if args.loci_input is None:
		args.loci_input = 'output/oe/genes_obs_exp.txt'
	if args.noncoding_input is None:
		args.noncoding_input = 'output/oe/noncoding_obs_exp.txt'
	
	for path in ['output/regional_constraint/%s' % args.out_dir]:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")

	print(datetime.datetime.now(), "Starting to analyze regional constraint!")
	
	# 待评估的基因座列表
	loci_list = ["MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-RNR1", "MT-RNR2"]
	protein_genes = ["MT-A", "MT-C", "MT-N"]
	
	res_to_pos = res_to_pos()
	pos_to_res = pos_to_res()

	regional_constraint(
		out_dir=args.out_dir, n_shuffles=args.n_shuffles, loci_type=args.loci_type, input_file=args.input,
		obs_value=args.obs, excluded_sites=args.exc_sites,
		fit_parameters=args.parameters, loci_input=args.loci_input, noncoding_input=args.noncoding_input,
		res_to_pos=res_to_pos, pos_to_res=pos_to_res, sig_threshold=args.sig_threshold)

	print(datetime.datetime.now(), "Script finished!")
