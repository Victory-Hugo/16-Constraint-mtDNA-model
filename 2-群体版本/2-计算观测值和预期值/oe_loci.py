import argparse
import datetime
from oe_functions import *
import os
'''该脚本用于计算线粒体基因组各基因或位点的携带者计数与预期值之比（obs/exp），并给出90%置信区间及p值。主要功能如下：
1. 读取注释文件，包含突变似然得分和观测携带者计数。
2. 读取线性模型参数文件，用于计算预期携带者计数。
3. 支持排除指定的碱基位点（如易出错位点）。
4. 针对每个位点和变异类型（同义、错义、终止、SNV）累加观测值和似然值。
5. 计算每个位点的obs/exp比值、置信区间和p值，并输出到指定文件。
6. 支持蛋白编码基因、RNA基因和非编码区域的分类处理。
7. 命令行参数支持输入文件、观测值列名、模型参数文件、输出前缀和排除位点列表。
主要函数：
- loci_oe: 计算每个位点的观测/预期比及置信区间，输出结果文件。
命令行参数说明：
- -input: 注释文件路径（含突变似然得分与观测携带者计数）
- -obs: 观测值列名
- -parameters: 线性模型参数文件路径
- -prefix: 输出文件名前缀
- -exc_sites: 需要排除的碱基位点列表
输出：
- 每个位点/基因的obs/exp比值、置信区间和p值，分别输出到蛋白编码基因和非编码区域结果文件。
'''

def loci_oe(input_file: str, obs_value: str, fit_parameters: str, output_prefix: str, excluded_sites: List[int]):
	"""计算基因和位点的观测值与预期值之比及其90%置信区间。

	:param input_file: 含突变似然得分与观测携带者计数的注释文件
	:param obs_value: 观测值列的列名
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param output_prefix: 输出文件名前缀
	:param excluded_sites: 需要在计算中排除的碱基位点列表
	"""
	file_genes = open('output/oe/%sgenes_obs_exp.txt' % output_prefix, "w")
	file_noncoding = open('output/oe/%snoncoding_obs_exp.txt' % output_prefix, "w")
	header = "locus	description	start	end	consequence	variant_count	observed	expected	obs:exp	lower_CI	upper_CI	pvalue"
	file_genes.write(header + '\n')
	file_noncoding.write(header + '\n')

	# 提取需要计算观测/预期比的位点
	mitomap_loci = {}
	for row in csv.DictReader(open('0-required_files/databases/mitomap_genome_loci.txt'), delimiter='\t'):
		if abs(int(row["Ending"]) - int(row["Starting"])) > 5:  # 跳过过小的位点
			mitomap_loci[(int(row["Starting"]), int(row["Ending"]))] = (row["Map_Locus"], row["Description"])
	
	# 处理位于两个基因中的变异
	per_gene = {}
	for row in csv.DictReader(open(
			'0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf'), delimiter='\t'):
		per_gene[(row["POS"], row["REF"], row["ALT"], row["SYMBOL"])] = row["Consequence"]
	
	for loci in mitomap_loci:
		# 初始化字典，使所有值为0
		loci_sum = initialize_sum_dict(identifier_list=['synonymous', 'missense', 'stop_gain', 'SNV'])
		# 针对每个位点，按变异类别累加观测值与似然值
		for row in csv.DictReader(open(input_file), delimiter='\t'):
			if int(row["POS"]) not in excluded_sites:
				mutation = row["REF"] + '>' + row["ALT"]
				region_to_use = 'ori' if (int(row["POS"]) in ori_region) else 'ref_exc_ori'
				name = ''
				# 如果该位置位于当前循环的位点中，loci[0] 为起点，loci[1] 为终点坐标
				if (loci[0] < loci[1]) and (loci[0] <= int(row["POS"]) <= loci[1]):
					if any(x in row["symbol"] for x in ["MT-A", "MT-C", "MT-N"]):  # 蛋白编码基因
						if ("synonymous" in per_gene[
							(row["POS"], row["REF"], row["ALT"], mitomap_loci[loci][0])]) and not any(
								x in row["consequence"] for x in more_severe_than_syn):
							name = 'synonymous'
						elif ("missense" in per_gene[
							(row["POS"], row["REF"], row["ALT"], mitomap_loci[loci][0])]) and not any(
								x in row["consequence"] for x in more_severe_than_missense):
							name = 'missense'
						elif ("stop_gain" in per_gene[(row["POS"], row["REF"], row["ALT"], mitomap_loci[loci][0])]) and (
								'stop_gained&start_lost' not in row["consequence"]):
							name = 'stop_gain'
					elif not any(x in row["symbol"] for x in ["MT-A", "MT-C", "MT-N"]):  # RNA 基因或非编码区域
						name = 'SNV'
				# 处理跨越人工断点的位点（仅适用于非编码位点）
				elif (loci[0] > loci[1]) and ((int(row["POS"]) >= loci[0]) or (int(row["POS"]) <= loci[1])):
					name = 'SNV'
				if name != '':
					loci_sum = sum_obs_likelihood(
						mutation=mutation, identifier=name, region=region_to_use,
						observed=row[obs_value], likelihood=row["Likelihood"],
						callable_samples=row["callable_samples"], dict=loci_sum)
		# 计算预期值
		for variant_type in ['synonymous', 'missense', 'stop_gain', 'SNV']:
			expected_carriers = calculate_exp(sum_dict=loci_sum, identifier=variant_type, fit_parameters=fit_parameters)
			if expected_carriers > 0:  # 跳过与该位点无关的变异后果
				observed_carriers = calculate_obs(identifier=variant_type, sum_dict=loci_sum)
				ratio_oe = observed_carriers / expected_carriers
				total_all = calculate_total(identifier=variant_type, sum_dict=loci_sum)
				(lower_CI, upper_CI) = calculate_CI(
					observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers)
				# 为保守起见，将观测值设为上置信界（upper_CI * expected_carriers）
				pvalue = calculate_pvalue(
					observed_carriers=(upper_CI * expected_carriers), total=total_all, expected_carriers=expected_carriers)
				
				print("Calculating values for", variant_type, "in", str(mitomap_loci[loci][0]))
				
				output = file_genes if ((variant_type != "SNV") or ("RNA" in mitomap_loci[loci][1])) else file_noncoding
				if (output == file_genes) or ((output == file_noncoding) and (round(expected_carriers) >= 10)):
					output.write(
						mitomap_loci[loci][0] + '\t' + mitomap_loci[loci][1] + '\t' + str(loci[0]) + '\t' + str(loci[1])
						+ '\t' + variant_type + '\t' + str(total_all) + '\t' + str(observed_carriers) + '\t' + str(expected_carriers)
						+ '\t' + str(ratio_oe) + '\t' + str(lower_CI) + '\t' + str(upper_CI) + '\t' + str(pvalue) + '\n')


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
	
	print(datetime.datetime.now(), "Calculate the observed:expected ratio for each gene/locus")
	
	loci_oe(
		input_file=args.input, obs_value=args.obs, fit_parameters=args.parameters, output_prefix=args.prefix,
		excluded_sites=args.exc_sites)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
	
