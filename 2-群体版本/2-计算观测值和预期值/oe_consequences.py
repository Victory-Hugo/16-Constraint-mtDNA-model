import argparse
import datetime
from oe_functions import *
import os
from typing import List
'''
该脚本用于计算不同功能变异类别（如同义、错义、终止密码子获得/丢失、rRNA、tRNA、间隔区等）的观测值与预期值之比（obs/exp），并给出90%置信区间。主要流程包括：
1. 读取包含突变似然得分和观测最大杂合度的注释文件。
2. 按功能变异类别累加观测值和似然值，排除指定的碱基位点。
3. 基于线性模型参数文件，计算每类变异的预期值。
4. 输出每类变异的obs/exp比值及置信区间到指定文件。
参数说明：
- input_file: 注释文件路径，包含突变似然得分与观测最大杂合度。
- obs_value: 观测值列名，指定用于计算的最大杂合度数据来源。
- fit_parameters: 线性模型参数文件路径，用于计算预期值。
- output_prefix: 输出文件名前缀。
- excluded_sites: 需要排除的碱基位点列表。
命令行参数支持自定义输入文件、观测值来源、模型参数、输出前缀及排除位点，默认使用gnomAD相关数据。输出结果包括每类变异的变异计数、观测值、预期值、obs/exp比值及其置信区间。

'''

def consequences_oe(
		input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算功能变异类别的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测最大杂合度的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file = open('output/oe/%sconsequences_obs_exp.txt' % output_prefix, "w")
	header = "consequence	variant_count	observed	expected	obs:exp	lower_CI	upper_CI"
	file.write(header + '\n')
	
	consequences = [
			"synonymous", "missense", "stop_gain", "start_lost", "stop_lost", "rRNA", "tRNA", "intergenic"]
	
	conseq_sum = initialize_sum_dict(identifier_list=consequences)
	# 针对每个位点，按变异类型累加观测值与似然值
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			name = ''
			# 选取最严重的后果标签
			if ("synonymous" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_syn):
				name = 'synonymous'
			elif ("missense" in row["consequence"]) and not any(x in row["consequence"] for x in more_severe_than_missense):
				name = 'missense'
			elif ("stop_gain" in row["consequence"]) and ('stop_gained&start_lost' not in row["consequence"]):
				name = 'stop_gain'
			elif "start_lost" in row["consequence"]:
				name = 'start_lost'
			# incomplete_terminal_codon_variant 等价于 stop_lost
			elif ("stop_lost" in row["consequence"]) or ("incomplete_terminal" in row["consequence"]):
				name = 'stop_lost'
			elif row["symbol"].startswith('MT-R'):
				name = 'rRNA'
			elif row["symbol"].startswith('MT-T'):
				name = 'tRNA'
			elif "intergenic" in row["consequence"]:
				name = 'intergenic'
			if name != '':
				conseq_sum = sum_obs_likelihood(
					mutation=mutation, identifier=name, region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"],
					callable_samples=row["callable_samples"], dict=conseq_sum)
	
	for consequence in consequences:
		calculate_oe(item=consequence, sum_dict=conseq_sum, fit_parameters=fit_parameters, file=file)


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
		# 排除 gnomAD 中的"artifact_prone_sites"：301、302、310、316、3107 和 16182（3107 已排除）
		# 这些位点在 gnomAD 中未调用，因此在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	for path in ['output/oe']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for each functional class of variation")
	
	consequences_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
	
