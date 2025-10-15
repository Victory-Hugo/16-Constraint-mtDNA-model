import argparse
import datetime
from oe_functions import *
import os
from typing import List

'''
该脚本用于对线粒体DNA参考序列中的中性变异进行模型校准，分别统计非 OriB-OriH 区域和 OriB-OriH 区域的观测最大杂合度、突变似然得分及计数。主要功能如下：
1. calibrate(input_file, obs_value, output_prefix, excluded_sites):
	- 对参考序列（不含 OriB-OriH 区域）中的中性变异进行累计统计。
	- 中性变异定义为系统发育树中出现的单倍群变异或 phyloP 评分处于最低十分位。
	- 按基因或功能区分组，输出每组的观测最大杂合度、突变似然得分和计数。
2. calibrate_ori(input_file, obs_value, output_prefix, excluded_sites):
	- 对 OriB-OriH 区域中的中性变异进行累计统计。
	- OriB-OriH 区域按约70 bp 划分为等长区块，分别统计每个区块的观测最大杂合度、突变似然得分和计数。
3. 主程序入口：
	- 解析命令行参数，包括输入文件、观测值列名、输出前缀和需排除的位点列表。
	- 设置默认参数并创建输出目录。
	- 依次调用 calibrate 和 calibrate_ori 函数，完成模型校准数据的准备。
依赖：
- oe_functions.py 中的辅助函数（如 initialize_sum_dict, sum_obs_likelihood）
- numpy, csv, argparse, datetime, os 等标准库
输出：
- 校准所用的中性变异列表（neutral_variants_used.txt）
- 各基因/区块的观测值与似然得分统计文件（loci_obs_vs_scores.txt, loci_obs_vs_scores_ori.txt）
适用场景：
- 线粒体DNA变异功能模型的校准与评估
- 变异分布与功能区块的统计分析
'''

def calibrate(input_file: str, obs_value: str, output_prefix: str, excluded_sites: List[int]):
	"""对参考序列（不含 OriB-OriH 区域）中的中性变异，累计观测最大杂合度、突变似然得分以及计数。

	:param input_file: 含突变似然得分与观测最大杂合度的注释文件
	:param obs_value: 观测值列的列名
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	# 首先提取基因/位点列表及其长度
	gene_length = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		for gene in row["symbol"].split(','):  # 拆分是为了处理同属两个基因的位点
			# 将控制区与其他非编码区域单独标记
			gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
			gene = 'other_non-coding' if gene == '' else gene
			# 排除 ori 区域及需剔除的位点
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in excluded_sites):
				if gene not in gene_length:
					gene_length[gene] = 1
				else:
					gene_length[gene] += 1
	
	# 现在初始化用于对各基因/位点求和的字典，并将所有值设为0
	calibration = initialize_sum_dict(identifier_list=list(gene_length.keys()))
	
	# 确定用于识别中性变异的 phyloP 阈值——取底部十分位
	phylop = []
	catch_list = []  # 由于每个位置有三行（3种可能SNV），用列表确保唯一位置
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # 第二个元素 [1] 即底部十分位

	# 写出所有使用到的中性变异
	f = open('output/calibration/neutral_variants_used.txt', "w")
	f.write('variant' + '\t' + 'position' + '\t' + 'source' + '\n')
	
	# 构建字典
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for gene in row["symbol"].split(','):  # 处理落在两个基因内的变异
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in excluded_sites):
				# 中性标准：在系统发育树中出现的单倍群变异或 phyloP 处于最低十分位
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# 将控制区与其他非编码区域单独标记
					gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
					gene = 'other_non-coding' if gene == '' else gene
					# 对每个基因/位点累加对应数值
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=gene, region='ref_exc_ori', dict=calibration,
						observed=row[obs_value], likelihood=row["Likelihood"])
					# 写入文件，位于两个基因的变异会记录两次
					reason = 'haplogroup_variant' if (int(row["in_phylotree"]) == 1) else 'lowest_decile_phyloP'
					f.write('m.' + str(row["POS"]) + row["REF"] + '>' + row["ALT"] + '\t' + str(row["POS"]) + '\t' + reason + '\n')
	
	# 写文件供绘图使用
	f = open('output/calibration/%sloci_obs_vs_scores.txt' % output_prefix, "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # 需要校准的两类突变分组
		for gene in gene_length:
			f.write(
				mut_group + '\t' + gene + '\t' +
				str(calibration[(gene, mut_group, 'obs_max_het', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'sum_LR', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'count', 'ref_exc_ori')]) + '\t' +
				str(gene_length[gene] / 3) + '\n')  # 计数包含每个位点的3种SNV，因此需除以3


def calibrate_ori(input_file: str, obs_value: str, output_prefix: str, excluded_sites: List[int]):
	"""对 OriB-OriH 区域中的中性变异，累计观测最大杂合度、突变似然得分以及计数。

	:param input_file: 含突变似然得分与观测最大杂合度的注释文件
	:param obs_value: 观测值列的列名
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	# 首先按照约70 bp（近似 tRNA 基因长度）将 ori 区域划分为等长区块
	ori_blocks = {}
	block = 1
	for i in range(0, len(ori_region), 70):
		if len(ori_region[i:i + 70]) == 70:
			ori_blocks["block_" + str(block)] = ori_region[i:i + 70]
		else:  # 如有需要，通过扩展最后一个区块确保所有 ori 位点被覆盖
			ori_blocks["block_" + str(block - 1)] = ori_blocks["block_" + str(block - 1)] + \
													ori_region[i:i + ori_region[-1]]
		block += 1
	
	# 初始化用于对各区块求和的字典，并将所有值设为0
	calibration = initialize_sum_dict(identifier_list=list(ori_blocks.keys()))
	
	# 确定用于识别中性变异的 phyloP 阈值——取底部十分位
	phylop = []
	catch_list = []  # 由于每个位置有三行（3种可能SNV），用列表确保唯一位置
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # 第二个元素 [1] 即底部十分位
	
	# 写出所有使用到的中性变异
	f = open('output/calibration/neutral_variants_used.txt', "a")
	
	# 构建字典
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for n_block in ori_blocks:
			# 若变异位于当前循环的 ori 区块
			if (int(row["POS"]) in ori_blocks[n_block]) and (int(row["POS"]) not in excluded_sites):
				# 中性标准：在系统发育树中出现的单倍群变异或 phyloP 处于最低十分位
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# 对每个区块累加对应数值
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=n_block, region='ori', dict=calibration,
						observed=row[obs_value], likelihood=row["Likelihood"])
					# 写入文件
					reason = 'haplogroup_variant' if (int(row["in_phylotree"]) == 1) else 'lowest_decile_phyloP'
					f.write('m.' + str(row["POS"]) + row["REF"] + '>' + row["ALT"] + '\t' + str(row["POS"]) + '\t' + reason + '\n')
	
	# 写文件供绘图使用
	f = open('output/calibration/%sloci_obs_vs_scores_ori.txt' % output_prefix, "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # 需要校准的两类突变分组
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
	
	# 设置 gnomAD 的默认值
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "gnomad_max_hl"
	if args.prefix is None:
		args.prefix = ""
	if args.exc_sites is None:
		# 排除 gnomAD 中的“artifact_prone_sites”：301、302、310、316、3107 和 16182（3107 已排除）
		# 这些位点在 gnomAD 中未调用，因此在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	for path in ['output/calibration']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Preparing files for model calibration!" + '\n')
	
	calibrate(input_file=args.input, obs_value=args.obs, output_prefix=args.prefix, excluded_sites=args.exc_sites)
	calibrate_ori(input_file=args.input, obs_value=args.obs, output_prefix=args.prefix, excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
