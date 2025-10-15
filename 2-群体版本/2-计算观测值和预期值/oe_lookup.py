import argparse
import datetime
from oe_functions import *
import os
'''
该脚本用于计算线粒体基因组指定坐标区间的观测最大杂合度与预期值之比（O/E ratio），并给出其90%置信区间（OEUF）。
主要功能：
- 读取含突变似然得分与观测最大杂合度的注释文件。
- 根据用户指定的坐标区间，累加观测值与似然值，排除指定的碱基位点。
- 利用线性方程参数文件，计算区间内的预期最大杂合度。
- 输出区间的观测值、预期值、O/E比值及OEUF（90%置信区间上界）。
参数说明：
- start (int): 要分析区间的起始坐标（必需）。
- end (int): 要分析区间的终止坐标（必需）。
- input_file (str): 含突变似然得分与观测最大杂合度的注释文件路径。
- obs_value (str): 观测最大杂合度对应的列名。
- fit_parameters (str): 包含线性方程系数与截距的文件路径。
- excluded_sites (List[int]): 需要在计算中排除的碱基位点列表。
命令行参数：
- -start: 起始坐标（必需）。
- -end: 终止坐标（必需）。
- -input: 注释文件路径（可选，默认值见代码）。
- -obs: 观测值列名（可选，默认值见代码）。
- -parameters: 线性方程参数文件路径（可选，默认值见代码）。
- -exc_sites: 排除的碱基位点列表（可选，默认值见代码）。
示例用法：
python3 oe_lookup.py -start 1 -end 16569
注意事项：
- 起始坐标不能大于终止坐标。
- 坐标范围应在1-16569之间（线粒体基因组范围）。
- 默认排除 gnomAD 中的 artifact-prone sites。
输出内容：
- 区间的观测最大杂合度、预期最大杂合度、O/E 比值及 OEUF（90%置信区间上界）。
'''

def oe_lookup(
		start: int, end: int,
		input_file: str, obs_value: str, fit_parameters: str, excluded_sites: List[int]):
	"""计算给定坐标区间的观测值与预期值之比及其90%置信区间。

	:param start: 起始坐标
	:param end: 终止坐标
	:param input_file: 含突变似然得分与观测最大杂合度的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	# 初始化字典，使所有值为0
	loci_sum = initialize_sum_dict(identifier_list=['all SNVs'])
	# 针对该区间内的每个位点，累加观测值与似然值
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in excluded_sites:
			mutation = row["REF"] + '>' + row["ALT"]
			region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
			# 此处包含所有 SNV，且起点必须小于终点
			if (int(row["POS"]) >= start) and (int(row["POS"]) <= end):
				loci_sum = sum_obs_likelihood(
					mutation=mutation, identifier='all SNVs', region=region_to_use,
					observed=row[obs_value], likelihood=row["Likelihood"], dict=loci_sum)
	# 计算比例与置信区间
	for variant_type in ['all SNVs']:
		exp_max_het = calculate_exp(sum_dict=loci_sum, identifier=variant_type, fit_parameters=fit_parameters)
		obs_max_het = calculate_obs(identifier=variant_type, sum_dict=loci_sum)
		ratio_oe = obs_max_het / exp_max_het
		total_all = calculate_total(identifier=variant_type, sum_dict=loci_sum)
		(lower_CI, upper_CI) = calculate_CI(
			obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het)
	print(
		"For interval", start, "to", end, "the observed and expected values are", obs_max_het, "and", exp_max_het,
		". The ratio is ", ratio_oe, "and the OEUF is", upper_CI, ".")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="计算给定坐标区间的观测值与预期值之比及其90%置信区间")
	parser.add_argument(
		"-input", type=str, help="含突变似然得分与观测最大杂合度的注释文件")
	parser.add_argument(
		"-obs", type=str, help="获取观测最大杂合度的人群数据集")
	parser.add_argument(
		"-parameters", type=str, help="包含线性方程系数与截距的文件路径")
	parser.add_argument(
		"-exc_sites", type=int, nargs='+', help="需要在计算中排除的碱基位点列表")
	parser.add_argument(
		"-start", type=int, required=True, help="要计算区间的起始坐标 (必需参数)")
	parser.add_argument(
		"-end", type=int, required=True, help="要计算区间的终止坐标 (必需参数)")
	args = parser.parse_args()
	
	# 设置 gnomAD 的默认值
	if args.input is None:
		args.input = 'output/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'
	if args.obs is None:
		args.obs = "gnomad_max_hl"
	if args.parameters is None:
		args.parameters = 'output/calibration/linear_model_fits.txt'
	if args.exc_sites is None:
		# 排除 gnomAD 中的"artifact_prone_sites"：301、302、310、316、3107 和 16182（3107 已排除）
		# 这些位点在 gnomAD 中未调用，因此在计算中剔除
		args.exc_sites = [301, 302, 310, 316, 16182]
	
	# 验证必需参数
	if args.start is None or args.end is None:
		print("错误：必须提供 -start 和 -end 参数来指定要分析的区间")
		print("用法示例: python3 oe_lookup.py -start 1 -end 16569")
		print("线粒体基因组坐标范围: 1-16569")
		exit(1)
	
	if args.start > args.end:
		print("错误：起始坐标不能大于终止坐标")
		exit(1)
	
	if args.start < 1 or args.end > 16569:
		print("警告：坐标超出线粒体基因组范围 (1-16569)")
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for a given interval")
	
	oe_lookup(
		start=args.start, end=args.end, input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
