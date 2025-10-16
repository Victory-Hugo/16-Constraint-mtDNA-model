import argparse
import csv
import datetime
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import sys
from typing import Dict, IO, List, Tuple

# 添加模块搜索路径
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '2-计算观测值和预期值'))
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

# 现在可以正常导入模块
import oe_functions
import regional_constraint


def calculate_fdr(rc_final: str, shuffle_final: str, shuffle_list: List[str], locus: str):
	"""计算各区域约束片段的FDR。

	:param rc_final: 真实比对中非重叠显著区域的文件路径
	:param shuffle_final: 所有随机置换中非重叠显著区域的文件路径
	:param shuffle_list: 需要遍历的随机置换编号列表
	:param locus: 基因座名称
	:return: FDR估计值字典
	"""
	# 根据真实比对得到的区域，整理待测试的p值阈值及长度（单位：bp）
	p_list = []
	for row in csv.DictReader(open(rc_final), delimiter='\t'):
		if row["locus"] == locus:
			p_list.append((row["locus"], float(row["kmer_pvalue"]), int(row["length"])))
			
	# 使用并行计算估算每个候选区域的FDR
	manager = mp.Manager()
	shuffle_dict = manager.dict()
	pool = mp.Pool(mp.cpu_count())
	for item in p_list:  # 针对每个条件统计假阳性
		pool.apply_async(parallelize_counting, args=(item, shuffle_list, shuffle_final, shuffle_dict))
	pool.close()
	pool.join()
	fdr = {}
	for item in shuffle_dict:
		fdr[item] = shuffle_dict[item] / len(shuffle_list)
		print("FDR estimate for region with following p-value and length", item, "is", fdr[item])

	return fdr


def parallelize_counting(
		item: Tuple[str, float, int], shuffle_list: List[str], shuffle_final: str,
		shuffle_dict: Dict[Tuple[str, float, int], float]):
	"""并行统计满足假阳性条件的置换次数。

	:param item: 需要评估的（基因座、p值阈值、长度）组合
	:param shuffle_list: 需遍历的随机置换编号列表
	:param shuffle_final: 所有置换中非重叠显著区域的文件路径
	:param shuffle_dict: 存储假阳性计数的字典，键为item
	:return: 更新后的shuffle_dict
	"""
	shuffle_counts = []  # 针对每个条件重置计数
	for shuffle_num in shuffle_list:  # 遍历所有置换
		count = 0
		shuffle_file_to_read = open(shuffle_final)
		for row in csv.DictReader(shuffle_file_to_read, delimiter='\t'):
			# 若该区域属于目标基因座，同时满足p值与长度阈值，并且来自当前置换
			if (row["locus"] == item[0]) and (float(row["kmer_pvalue"]) <= item[1]) and (int(row["length"]) == item[2]) and (
					row["shuffle"] == shuffle_num.split('-')[0]) and (row["out_dir"].split('/')[-1] == shuffle_num.split('-')[1]):
				count = 1
		shuffle_file_to_read.close()
		shuffle_counts.append(count)
	# 将满足条件的置换数量写入字典
	shuffle_dict[item] = sum(shuffle_counts)
	
	return shuffle_dict


def filter_by_FDR(
	locus: str, significant_kmers: str, shuffle: str, dict: Dict[Tuple[int, int, str], Tuple[float, int, str]],
		pvalue_threshold: float, filters: Dict[int, float]):
	"""剔除未通过FDR阈值的候选区域，并重新执行贪心算法。

	:param locus: 基因座名称
	:param significant_kmers: 含全部显著（可能重叠）区域的文件路径
	:param shuffle: 当前筛选的随机置换编号
	:param dict: 字典，键为区域起止坐标，值为（p值、长度、置换编号）
	:param pvalue_threshold: 本次分析允许的最高p值
	:param filters: 字典，键为长度，值为对应的p值过滤阈值
	:return: 贪心算法输出字典，键为（kmer起点、终点、置换编号），值为（p值、长度、置换编号）
	"""
	# 过滤未通过FDR阈值的显著kmer，并重新构建贪心算法输入
	outlier_kmers = {}
	for row in csv.DictReader(open(significant_kmers), delimiter='\t'):
		if (row["locus"] == locus) and (float(row["pvalue"]) < pvalue_threshold):  # 全局p值阈值
			# 输入文件中的蛋白基因长度以密码子为单位，此处转换为bp
			length = (int(row["length"]) * 3) if any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) else int(row["length"])
			# 满足真实比对或当前置换
			if (shuffle == "0") or (
					(row["shuffle"] == shuffle.split('-')[0]) and (row["out_dir"].split('/')[-1] == shuffle.split('-')[1])):
				include = "yes"
				for filter_key in sorted(filters):  # 按长度排序依次应用阈值
					if (length == filter_key) and (float(row["pvalue"]) >= filters[filter_key]):
						include = "no"  # 标记为需剔除
				if include == "yes":  # 作为贪心算法的输入
					outlier_kmers[(int(row["kmer_start"]), int(row["kmer_end"]))] = (float(row["pvalue"]), length)
	
	# 使用并行方式运行贪心算法
	if len(outlier_kmers) > 0:
		# 为便于处理控制区，将1-576位置平移至16569之后
		new_dict = {}
		for key in sorted(outlier_kmers):
			key_start = (16569 + key[0]) if (locus == "CR" and key[0] < 577) else key[0]
			key_end = (16569 + key[1]) if (locus == "CR" and key[1] < 577) else key[1]
			new_dict[(key_start, key_end)] = (outlier_kmers[key][0], outlier_kmers[key][1])
		df = pd.DataFrame.from_dict(new_dict, orient='index', columns=['pvalue', 'length'])
		df.index = pd.MultiIndex.from_tuples(df.index, names=['kmer_start', 'kmer_end'])
		df = df.reset_index()
		pool = mp.Pool(mp.cpu_count())
		# 该步骤返回贪心算法的结果，dict参数即函数输入
		for key in list(sorted(new_dict)):  # 遍历每个显著kmer
			pool.apply_async(reextract_nonoverlapping, args=(key, df, dict, shuffle))
		pool.close()
		pool.join()
		
	if dict is not None:  # 避免返回None
		return dict


def reextract_nonoverlapping(
		key: Tuple[int, int],
		df: pd.DataFrame,
		filtered_outlier_kmers: Dict[Tuple[int, int, str], Tuple[float, int, str]],
		shuffle: str):
	"""按最低p值（必要时比较长度）筛选互不重叠的kmer。

	:param key: kmer起止坐标及置换编号；控制区已将1-576位置平移至16569之后
	:param df: 包含全部显著kmer的pandas数据框
	:param filtered_outlier_kmers: 字典，键为kmer起止坐标，值为（p值、长度、置换编号）
	:param shuffle: 当前置换编号
	:return: 更新后的filtered_outlier_kmers
	"""
	kmer_start = key[0]
	kmer_end = key[1]
	# 筛选所有与当前kmer重叠的片段，并按p值升序、长度降序排序
	# 过滤条件确保两个区间确实重叠
	pre_filter_df = df[
		(df['kmer_start'] <= np.int64(kmer_end)) &
		(df['kmer_end'] >= np.int64(kmer_start))].sort_values(['pvalue', 'length'], ascending=(True, False))
	
	# 提取该kmer的p值与长度
	kmer_pvalue = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))]['pvalue'].iloc[0]
	kmer_length = pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))]['length'].iloc[0]
	# 将当前kmer从数据框中移除，便于比较其余候选
	filter_df = pre_filter_df.drop(pre_filter_df[(pre_filter_df['kmer_start'] == np.int64(kmer_start)) & (
			pre_filter_df['kmer_end'] == np.int64(kmer_end))].index)
		
	# 根据以下条件筛选互不重叠的kmer
	# (1) 如果不存在重叠片段
	if filter_df.shape[0] == 0:
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)
		# 保留置换编号以区分不同置换中相同坐标的情况

	# (2) 当前kmer的p值低于全部重叠片段
	elif kmer_pvalue < filter_df['pvalue'].iloc[0]:
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)

	# (3) 当前kmer的p值等于最低p值，但长度更长或相同
	elif (kmer_pvalue == filter_df['pvalue'].iloc[0]) and (kmer_length >= filter_df['length'].iloc[0]):
		filtered_outlier_kmers[(key[0], key[1], shuffle)] = (kmer_pvalue, kmer_length, shuffle)

	if filtered_outlier_kmers is not None:  # 避免返回None
		return filtered_outlier_kmers


def annotate_final(
		key: Tuple[int, int, str], dict_key: Tuple[int, int, str], locus: str, output_file: IO,
		input_file: str, obs_value: str, excluded_sites: List[int], fit_parameters: str,
	pos_to_res: Dict[Tuple[str, int], int],
	real_filtered_outlier_kmers: Dict[Tuple[int, int, str], Tuple[float, int, str]],
	fdr_dict: Dict[Tuple[str, float, int], float]):
	"""为最终确认的区域约束片段添加注释并写入文件。

	:param key: mtDNA坐标下（kmer起点、终点、置换编号）的元组，真实比对置换编号为"0"
	:param dict_key: real_filtered_outlier_kmers中的键，与key一致但用于查找原始结果
	:param locus: 基因座名称
	:param output_file: 输出文件句柄
	:param input_file: 包含突变似然评分和观测最大杂合度的注释文件路径
	:param obs_value: 观测值列标题
	:param excluded_sites: 需要排除的位点列表
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param pos_to_res: 将碱基位置转换为蛋白残基的字典
	:param real_filtered_outlier_kmers: 字典，键为（kmer起点、终点、置换编号），值为（p值、长度、置换编号）
	:param fdr_dict: 各长度和p值组合对应的FDR估计
	"""
	# 计算最终区域指标所需的辅助字典和列表
	per_gene = {}
	for row in csv.DictReader(open(
			'0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		# 控制区命名需与输入一致
		row["SYMBOL"] = "CR" if ((int(row["POS"]) <= 576) or (int(row["POS"]) >= 16024)) else row["SYMBOL"]
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]

	# 记录该区域内的所有位置
	kmer_coords = list(range(key[0], key[1] + 1)) if (key[0] < key[1]) else \
		list(range(key[0], 16570)) + list(range(1, key[1] + 1))  # 处理跨越m.16569-1的情况

	per_kmer = oe_functions.initialize_sum_dict(identifier_list=[(key[0], key[1])])
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in oe_functions.ori_region) else 'ref_exc_ori'
			if int(row["POS"]) in kmer_coords:
				# RNA基因和非编码区域保留全部碱基；蛋白编码基因仅保留最严重为错义的变异
				if not any(x in locus for x in ["MT-A", "MT-C", "MT-N"]) or (
						("missense" in per_gene[(row["POS"], row["REF"], row["ALT"], locus)]) and not any(
					x in row["consequence"] for x in oe_functions.more_severe_than_missense)):
					per_kmer = oe_functions.sum_obs_likelihood(
						mutation=mutation, identifier=(key[0], key[1]), region=region_to_use,
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=per_kmer)
	
	# 计算观测值与期望值
	obs_max_het = oe_functions.calculate_obs(sum_dict=per_kmer, identifier=(key[0], key[1]))
	obs_max_het = 0.0 if obs_max_het < 0.001 else obs_max_het  # 处理极小数值
	exp_max_het = oe_functions.calculate_exp(
		sum_dict=per_kmer, identifier=(key[0], key[1]), fit_parameters=fit_parameters)
	total_all = oe_functions.calculate_total(sum_dict=per_kmer, identifier=(key[0], key[1]))
	ratio_oe = obs_max_het / exp_max_het
	(lower_CI, upper_CI) = oe_functions.calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het)
	# 输出时写入残基坐标等信息
	output_file.write(
		str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
		str(pos_to_res[(locus, key[0])] if (locus, key[0]) in pos_to_res else '') + '\t' +
		str(pos_to_res[(locus, key[1])] if (locus, key[1]) in pos_to_res else '') + '\t' +
		str(obs_max_het) + '\t' + str(exp_max_het) + '\t' + str(ratio_oe) + '\t' +
		str(lower_CI) + '\t' + str(upper_CI) + '\t' + str(total_all) + '\t' +
		str(real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][0]) + '\t' +
		str(real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][1]) + '\t' +
		str(fdr_dict[(
			locus, real_filtered_outlier_kmers[(
				dict_key[0], dict_key[1], "0")][0], real_filtered_outlier_kmers[(dict_key[0], dict_key[1], "0")][1])]) + '\n')
	output_file.flush()  # 确保及时写入文件
	

def apply_FDR_filter(
		rc_final: str, shuffle_final: str, fdr_threshold: float, real_all: str, shuffle_all: str,
		shuffle_list: List[str], input_file: str, obs_value: str, excluded_sites: List[int], fit_parameters: str,
		pvalue_threshold: float, loci_type: str, file: IO):
	"""输出通过FDR筛选的高置信度区域约束区段。

	:param rc_final: 真实比对中非重叠显著区域的文件
	:param shuffle_final: 所有置换中非重叠显著区域的文件
	:param fdr_threshold: FDR阈值
	:param real_all: 真实比对中重叠显著区域的文件
	:param shuffle_all: 所有置换中重叠显著区域的文件
	:param shuffle_list: 需遍历的置换编号列表
	:param input_file: 包含突变似然评分和观测最大杂合度的注释文件路径
	:param obs_value: 观测值列标题
	:param excluded_sites: 需要排除的位点列表
	:param fit_parameters: 包含线性模型系数和截距的文件路径
	:param pvalue_threshold: 全局允许的最高p值
	:param loci_type: 'protein'表示蛋白基因，'rn'表示RNA及非编码
	:param file: 用于写入结果的输出文件
	"""
	# 逐个基因座处理
	loci_list = [
		"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L",
		"MT-ND5", "MT-ND6"] if (loci_type == "protein") else ["MT-RNR1", "MT-RNR2"]
	
	# 初始化所需字典
	pos_to_res = regional_constraint.pos_to_res()
	res_to_pos = regional_constraint.res_to_pos()
	# 记录初始文件路径
	initial_rc_file = rc_final
	initial_shuffle_file = shuffle_final
	manager = mp.Manager()
	
	for locus in loci_list:
		iteration = 0
		fdr_dict = calculate_fdr(
			rc_final=initial_rc_file, shuffle_final=initial_shuffle_file, shuffle_list=shuffle_list, locus=locus)
		
		result = "pass"
		for item in fdr_dict:
			result = "fail" if (fdr_dict[item] >= fdr_threshold) else result
			
		if result == "pass":  # 若该基因座全部候选首次即通过阈值
			print(locus, " passed!")
			real_filtered_outlier_kmers = {}
			for row in csv.DictReader(open(initial_rc_file), delimiter='\t'):
				if row["locus"] == locus:
					real_filtered_outlier_kmers[
						(int(row["kmer_start"]), int(row["kmer_end"]), "0")] = (float(row["kmer_pvalue"]), int(row["length"]), "0")
			# write to file
			for key in real_filtered_outlier_kmers:
				annotate_final(
					key=key, dict_key=key, locus=locus, output_file=file, input_file=input_file,
					obs_value=obs_value, excluded_sites=excluded_sites, fit_parameters=fit_parameters,
					pos_to_res=pos_to_res, real_filtered_outlier_kmers=real_filtered_outlier_kmers, fdr_dict=fdr_dict)
		
		# 若仍有候选未通过FDR阈值，则为该基因座生成过滤条件
		filters = {}
		while result == "fail":  # 只要存在未通过的候选就继续循环
			iteration += 1
			print(datetime.datetime.now(), locus, "needs filtering")
			
			# 构建需应用的过滤条件
			for key in fdr_dict:  # 键格式为（locus, p值, 长度）
				if fdr_dict[key] >= fdr_threshold:
					# 按长度存储p值阈值
					if key[2] not in filters:
						filters[key[2]] = key[1]
					# 如已存在该长度的条件且当前p值更低，则覆盖
					elif (key[2] in filters) and (key[1] < filters[key[2]]):
						filters[key[2]] = key[1]
						
			# 先对真实比对进行过滤并重新执行贪心算法
			real_filtered_outlier_kmers = manager.dict()
			real_filtered_outlier_kmers = filter_by_FDR(
				locus=locus, significant_kmers=real_all, shuffle="0", dict=real_filtered_outlier_kmers,
				pvalue_threshold=pvalue_threshold, filters=filters)
			
			# 再对各置换结果重复过滤与贪心选择
			shuffle_filtered_outlier_kmers = manager.dict()
			for shuffle in shuffle_list:
				shuffle_filtered_outlier_kmers = filter_by_FDR(
					locus=locus, significant_kmers=shuffle_all, shuffle=shuffle, dict=shuffle_filtered_outlier_kmers,
					pvalue_threshold=pvalue_threshold, filters=filters)
			
			# 将贪心算法结果写入文件，为FDR估计准备输入（长度单位为bp）
			rc_final = 'output/regional_constraint/temp_regional_constraint_intervals.txt'
			file2 = open(rc_final, "w")
			file2.write("kmer_start	kmer_end	locus	kmer_pvalue	length" + '\n')
			for key in sorted(real_filtered_outlier_kmers):  # 并行会改变顺序，需排序
				file2.write(
					str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
					str(real_filtered_outlier_kmers[key][0]) + '\t' +
					str(real_filtered_outlier_kmers[key][1]) + '\t' + '\n')
				file2.flush()  # 确保及时写入文件
			
			shuffle_final = 'output/regional_constraint/temp_shuffle_regional_constraint_intervals.txt'
			file3 = open(shuffle_final, "w")
			file3.write("kmer_start	kmer_end	locus	kmer_pvalue	length	shuffle	out_dir" + '\n')
			for key in sorted(shuffle_filtered_outlier_kmers):
				file3.write(
					str(key[0]) + '\t' + str(key[1]) + '\t' + str(locus) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][0]) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][1]) + '\t' +
					str(shuffle_filtered_outlier_kmers[key][2].split('-')[0]) + '\t' +
					str('out_dir/' + shuffle_filtered_outlier_kmers[key][2].split('-')[1]) + '\n')
				file3.flush()  # 确保及时写入文件
		
			# 重新计算FDR
			fdr_dict = calculate_fdr(
				rc_final=rc_final, shuffle_final=shuffle_final, shuffle_list=shuffle_list, locus=locus)
			
			result = "pass"
			for item in fdr_dict:
				result = "fail" if (fdr_dict[item] >= fdr_threshold) else result
				
			if result == "pass":  # 若全部通过，则计算obs/exp等指标并写出
				for key in real_filtered_outlier_kmers:
					# 写出前需执行以下步骤
					# 处理控制区（已将1-576位置平移至16569之后）
					key_start = (key[0] - 16569) if (key[0] > 16569) else key[0]
					key_end = (key[1] - 16569) if (key[1] > 16569) else key[1]
					key = (key_start, key_end, "0")
					# 将残基坐标转换为mtDNA碱基坐标（首个密码子的首位至最后密码子的末位）
					# 因离群字典以残基索引存储起止信息，需进行转换
					residue_key = key
					if any(x in locus for x in ["MT-A", "MT-C", "MT-N"]):
						key = (min(res_to_pos[(locus, key[1])]), max(res_to_pos[(locus, key[0])])) if (locus == "MT-ND6") \
							else (min(res_to_pos[(locus, key[0])]), max(res_to_pos[(locus, key[1])]))
					annotate_final(
						key=(key[0], key[1], "0"), dict_key=residue_key, locus=locus, output_file=file, input_file=input_file,
						obs_value=obs_value, excluded_sites=excluded_sites, fit_parameters=fit_parameters,
						pos_to_res=pos_to_res, real_filtered_outlier_kmers=real_filtered_outlier_kmers, fdr_dict=fdr_dict)
				print("All done for locus", locus)
				continue  # 继续处理下一个基因座
		# 若仍有区域未通过，本循环将继续执行
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	parser.add_argument(
		"-obs", type=str, help="Population dataset from which observed maximum heteroplasmy is obtained")
	parser.add_argument(
		"-parameters", type=str, help="File with parameters from linear model to calculate expected")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="List of base positions to exclude from calculations")
	args = parser.parse_args()
	
	# 设置默认参数，适用于gnomAD
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "carrier_count"
	if args.parameters is None:
		args.parameters = 'output/calibration/linear_model_fits.txt'
	if args.exc_sites is None:
		# 排除gnomAD中易受技术伪影影响的位置：301、302、310、316、3107和16182（3107已在前面排除）
		# 这些位点未在gnomAD中检测到，因此需在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	print(datetime.datetime.now(), "Starting to filter regional constraint by FDR!")
	print("This may take a while...")
	
	# 创建最终结果文件
	file = open('output/regional_constraint/final_regional_constraint_intervals.txt', "w")
	file.write(
		"start	end	locus	protein_pos_start	protein_pos_end	obs_max_het	exp_max_het	ratio_oe	"
		"lower_CI	upper_CI	total	pvalue	length	fdr" + '\n')
	
	# 汇总需要使用的置换编号——先处理蛋白基因
	number_shuffles = 70
	number_out_dir = 14
	shuffle_list = []
	for shuffle_number in range(1, number_shuffles + 1):
		for out_dir in range(1, number_out_dir + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# 手动补齐最后一个置换，使总数达到1000（输出目录15）
	for shuffle_number in range(1, 20 + 1):
		shuffle_list.append(str(shuffle_number) + '-' + str(15))
	
	apply_FDR_filter(
		rc_final='output/regional_constraint/real_alignment/regional_constraint_intervals.txt',
		shuffle_final='output/regional_constraint/protein_shuffle/all_shuffles_nonoverlapping.txt',
		fdr_threshold=0.10,
		real_all='output/regional_constraint/real_alignment/all_significant_kmers.txt',
		shuffle_all='output/regional_constraint/protein_shuffle/all_sig_shuffles.txt',
		shuffle_list=shuffle_list,
		input_file=args.input, obs_value=args.obs, excluded_sites=args.exc_sites, fit_parameters=args.parameters,
		pvalue_threshold=0.005, loci_type="protein", file=file)
	# 根据真实比对候选的最高p值设置pvalue_threshold
	
	# 汇总RNA置换编号（分三批完成）
	number_shuffles = 19
	number_out_dir = 20
	shuffle_list = []
	for shuffle_number in range(1, number_shuffles + 1):
		for out_dir in range(1, number_out_dir + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# 批次2
	for shuffle_number in range(1, 10 + 1):
		for out_dir in range(21, 32 + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	# 批次3
	for shuffle_number in range(1, 10 + 1):
		for out_dir in range(33, 82 + 1):
			shuffle_list.append(str(shuffle_number) + '-' + str(out_dir))
	
	print("Calculating for", len(shuffle_list), "shuffles")
			
	apply_FDR_filter(
		rc_final='output/regional_constraint/real_alignment/regional_constraint_intervals.txt',
		shuffle_final='output/regional_constraint/rn_shuffle/all_shuffles_nonoverlapping.txt',
		fdr_threshold=0.10,
		real_all='output/regional_constraint/real_alignment/all_significant_kmers.txt',
		shuffle_all='output/regional_constraint/rn_shuffle/all_sig_shuffles.txt',
		shuffle_list=shuffle_list,
		input_file=args.input, obs_value=args.obs, excluded_sites=args.exc_sites, fit_parameters=args.parameters,
		pvalue_threshold=0.005, loci_type="rn", file=file)
	
	# 完成后删除临时文件
	os.remove('output/regional_constraint/temp_regional_constraint_intervals.txt')
	os.remove('output/regional_constraint/temp_shuffle_regional_constraint_intervals.txt')
	
	print(datetime.datetime.now(), "Script finished!")
	
