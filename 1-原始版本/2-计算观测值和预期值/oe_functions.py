import csv
import numpy as np
from scipy.stats import beta
from typing import Dict, List, TextIO, Tuple, Union

'''
该模块用于计算线粒体DNA变异的观测最大杂合度、预期值、观测/预期比值以及置信区间，主要用于评估变异的约束性。核心功能包括：
- initialize_sum_dict: 初始化用于统计不同变异类别的观测值、似然得分和计数的字典。
- sum_obs_likelihood: 累加某类变异的观测最大杂合度、似然得分和计数。
- linear_model: 应用线性模型计算预期观测最大杂合度总和，并进行夹取处理。
- calibrate: 根据线性模型参数，返回某类变异的预期观测最大杂合度总和。
- calculate_exp: 读取线性模型参数文件，计算某类变异的预期观测最大杂合度总和。
- calculate_obs: 计算某类变异的观测最大杂合度总和。
- calculate_total: 计算某类变异的可能变异总数。
- calculate_CI: 基于Beta分布，计算观测/预期比值的90%置信区间。
- calculate_pvalue: 基于Beta分布，计算观测值显著低于预期值的概率。
- calculate_oe: 综合计算观测值、预期值、观测/预期比值及置信区间，并输出结果。
模块还定义了用于区分OriB-OriH区域的变量，以及蛋白编码基因变异后果的严重性分级列表。

'''
def initialize_sum_dict(identifier_list: List[Union[str, Tuple[int, int]]]):
	"""生成一个用于对一组变异的观测最大杂合度、突变似然得分和计数做求和的字典。

	:param identifier_list: 需要求和的每个变异类别的标识符
	:return: 一个以标识符、突变分组、需要求和的值以及包含变异的区域组成的元组为键、以0为初始值的字典
	"""
	# 初始化字典，使所有值为0
	dict = {}
	for identifier in identifier_list:
		for mut_group in ['G>A_and_T>C', 'other']:  # 要拟合的两个突变分组
			for value in ['obs_max_het', 'sum_LR', 'count']:  # 需要求和的三个数值
				for region in ['ref_exc_ori', 'ori']:  # 可能应用的两种突变模型
					dict[identifier, mut_group, value, region] = 0
	return dict


def sum_obs_likelihood(
		mutation: str, identifier: Union[str, Tuple[int, int]], region: str,
		observed: Union[str, float], likelihood: Union[str, float],
		dict: Dict[Tuple[str, str, str, str], Union[float, int]]):
	"""对某一类变异的观测最大杂合度、突变似然得分和计数进行求和。

	:param mutation: Ref>Alt 形式的突变类型
	:param identifier: 需要求和的变异类别标识符
	:param region: 取值应为 'ref_exc_ori' 或 'ori'，用于指示应用哪种突变模型
	:param observed: 用于观测最大杂合度的取值
	:param likelihood: 用于突变似然得分的取值
	:param dict: 需要更新的字典名称
	:return: 一个以标识符、突变分组、需要求和的值和变异所在区域组成的元组为键、以累计值为值的字典
	"""
	# 高度可变的 G>A 和 T>C 单独拟合
	mut_group = 'G>A_and_T>C' if (mutation == 'G>A' or mutation == 'T>C') else 'other'
	# 注意字典应已初始化（所有预期键的初始值均为0）
	dict[(identifier, mut_group, 'obs_max_het', region)] += float(observed)
	dict[(identifier, mut_group, 'sum_LR', region)] += float(likelihood)
	dict[(identifier, mut_group, 'count', region)] += 1
	return dict


def linear_model(count: int, intercept: float, coefficient: float, sum_lr: float):
	"""将针对中性变异拟合的线性模型应用于数据，以计算预期的观测最大杂合度总和。
	会应用一个夹取操作，使预期值不会小于0，也不会大于可能的最大值。

	:param count: 需要计算预期值的变异数量
	:param intercept: 线性方程的截距
	:param coefficient: 线性方程的系数
	:param sum_lr: 需要计算预期值的变异的突变似然得分之和
	:return: 预期的观测最大杂合度总和
	"""
	return max(float(0), min(float(count), intercept + (coefficient * sum_lr)))


def calibrate(
		sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]], mut_group: str, identifier: str, region: str,
		parameters: Dict[Tuple[str, str, str], float]):
	"""应用线性模型，返回预期的观测最大杂合度总和。

	:param sum_dict: 包含计算所需的变异数量和突变得分总和的字典
	:param mut_group: 需要校准的突变分组，取值为 'G>A_and_T>C' 或 'other'
	:param identifier: 需要校准的变异类别标识符
	:param region: 取值应为 'ref_exc_ori' 或 'ori'，用于指示使用哪种校准
	:param parameters: 包含线性方程系数和截距的字典，用于执行校准
	:return: 线性模型的输出，即相应突变分组和区域的预期观测最大杂合度总和
	"""
	return linear_model(
		count=sum_dict[(identifier, mut_group, 'count', region)],
		intercept=parameters[(mut_group, region, 'intercept')],
		coefficient=parameters[(mut_group, region, 'coefficient')],
		sum_lr=sum_dict[(identifier, mut_group, 'sum_LR', region)])


def calculate_exp(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]], fit_parameters: str):
	"""计算某一类变异的预期观测最大杂合度总和。

	:param sum_dict: 包含计算所需的变异数量和突变得分总和的字典
	:param identifier: 需要计算预期值的变异类别标识符
	:param fit_parameters: 包含线性方程系数和截距的文件路径
	:return: 预期的观测最大杂合度总和
	"""
	# 创建包含线性方程系数和截距的字典，用于计算预期值
	parameters = {}
	for row in csv.DictReader(open(fit_parameters), delimiter='\t'):
		parameters[(row["mutation_group"], row["region"], row["item"])] = float(row["value"])
	
	# 对 'G>A_and_T>C' 与 'other' 两类突变分组以及 'ref_exc_ori' 与 'ori' 两种区域分别计算预期值
	# 如果某一分组没有变异，calibrate 函数将因为夹取而返回0
	exp = calibrate(
			sum_dict=sum_dict, mut_group='G>A_and_T>C', identifier=identifier, region='ref_exc_ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='other', identifier=identifier, region='ref_exc_ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='G>A_and_T>C', identifier=identifier, region='ori', parameters=parameters) + \
		calibrate(
			sum_dict=sum_dict, mut_group='other', identifier=identifier, region='ori', parameters=parameters)
	return exp


def calculate_obs(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]]):
	"""计算某一类变异的观测最大杂合度总和。

	:param identifier: 需要计算观测值的变异类别标识符
	:param sum_dict: 包含观测最大杂合度的字典
	:return: 跨越两类突变分组和两类区域类型的观测最大杂合度总和
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'obs_max_het', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'obs_max_het', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'obs_max_het', 'ori')] + sum_dict[
		(identifier, 'other', 'obs_max_het', 'ori')]


def calculate_total(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]]):
	"""计算某一组变异的可能变异总数。

	:param identifier: 需要计算的变异类别标识符
	:param sum_dict: 包含观测最大杂合度的字典
	:return: 跨越两类突变分组和两类区域类型的可能变异总数
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'count', 'ori')] + sum_dict[
		(identifier, 'other', 'count', 'ori')]


def calculate_CI(obs_max_het: float, total: int, exp_max_het: float, max_parameter: float = 2.5):
	"""使用 Beta 分布计算观测值与预期值之比的90%置信区间。

	:param obs_max_het: 观测最大杂合度总和
	:param total: 可能的变异数量
	:param exp_max_het: 预期最大杂合度总和
	:param max_parameter: 变化参数的最大值，默认 2.5
	:return: lower_CI 和 upper_CI，即90%置信区间的上下界
	"""
	# 针对给定的观测值与预期值计算 Beta 分布的密度
	# 使用预期值乘以不同的变化参数，并设定需要评估的参数范围
	varying_parameters = np.arange(0.001, max_parameter, 0.001).tolist()
	
	# 将观测值表示为可能值的比例，使其限定在0-1之间
	prop_poss_maxhet = obs_max_het / total
	# 如果比例为0或1，则调整到合理范围内以便评估；为保守起见使用最小可能观测值
	if prop_poss_maxhet == 0:
		prop_poss_maxhet = 0.1 / total
	elif prop_poss_maxhet == 1:
		prop_poss_maxhet = 1 - (0.1 / total)
	
	beta_distrib = []
	for vp in varying_parameters:
		a = (exp_max_het * vp) + 1  # “成功”次数
		b = total - (exp_max_het * vp) + 1  # “失败”次数
		beta_distrib.append(beta.pdf(prop_poss_maxhet, a, b))

	# 计算该函数的累积分布
	cumulative_beta = []
	for value in np.cumsum(beta_distrib):
		cumulative_beta.append(value / max(np.cumsum(beta_distrib)))

	# 在累计概率等于5%与95%的位置提取对应的变化参数，得到90%置信区间
	lower_list, upper_list = [], []
	for value in cumulative_beta:
		if value < 0.05:
			lower_list.append(cumulative_beta.index(value))
		if value > 0.95:
			upper_list.append(cumulative_beta.index(value))
	
	if (obs_max_het > 0) and (len(lower_list) > 0):
		lower_idx = max(lower_list)  # 概率小于0.05时的最大位置
		lower_CI = varying_parameters[lower_idx]  # 第5百分位
	else:
		lower_CI = 0
	upper_idx = min(upper_list)  # 概率大于0.95时的最小位置
	upper_CI = varying_parameters[upper_idx]  # 第95百分位
	
	return lower_CI, upper_CI


def calculate_pvalue(obs_max_het: float, total: int, exp_max_het: float):
	"""利用 Beta 分布判断观测值是否显著低于预期值。

	:param obs_max_het: 观测最大杂合度总和
	:param total: 可能的变异数量
	:param exp_max_het: 预期最大杂合度总和
	:return: less_than_pvalue，即在给定预期值的情况下，观测到小于或等于该观测值的概率
	"""
	# 将观测值表示为可能值的比例，使其限定在0-1之间
	prop_poss_maxhet = obs_max_het / total
	# 如果比例为0或1，则调整到合理范围内以便评估；为保守起见使用最小可能观测值
	if prop_poss_maxhet == 0:
		prop_poss_maxhet = 0.1 / total
	elif prop_poss_maxhet == 1:
		prop_poss_maxhet = 1 - (0.1 / total)
	a = exp_max_het + 1
	b = total - exp_max_het + 1
	less_than_pvalue = beta.cdf(prop_poss_maxhet, a, b)
	return less_than_pvalue


def calculate_oe(
		item: Union[str, Tuple[int, int]], sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		fit_parameters: str, file: Union[TextIO, None], output: str = "write", max_parameter: float = 2.5):
	"""计算观测值、预期值、观测与预期之比以及90%置信区间。

	:param item: 需要计算预期值的变异类别标识符
	:param sum_dict: 包含计算所需的可能变异数量和突变得分总和的字典
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param file: 当 output 为 "write" 时用于写入结果的文件
	:param output: 指定输出到文件或保存为字典，默认 "write"，可选 "dict"
	:param max_parameter: calculate_CI 使用的变化参数上限，默认 2.5
	"""
	exp_max_het = calculate_exp(sum_dict=sum_dict, identifier=item, fit_parameters=fit_parameters)
	obs_max_het = calculate_obs(sum_dict=sum_dict, identifier=item)
	ratio_oe = obs_max_het / exp_max_het
	total_all = calculate_total(sum_dict=sum_dict, identifier=item)

	(lower_CI, upper_CI) = calculate_CI(
		obs_max_het=obs_max_het, total=total_all, exp_max_het=exp_max_het, max_parameter=max_parameter)
	
	print("Calculating values for", item)
	if output == "write":
		file.write(
			item + '\t' + str(total_all) + '\t' + str(obs_max_het) + '\t' + str(exp_max_het)
			+ '\t' + str(ratio_oe) + '\t' + str(lower_CI) + '\t' + str(upper_CI) + '\n')
	elif output == "dict":
		return total_all, obs_max_het, exp_max_het, ratio_oe, lower_CI, upper_CI


# 待使用的变量

# ori 指代在人工断点两侧具有已知突变特征差异的 OriB-OriH 区域（m.16197-191）
ori_region = list(range(16197, 16570)) + list(range(1, 191 + 1))

# 蛋白编码基因可能的变异后果，按严重程度排序：
# stop_gained（获得终止密码子）
# stop_gained&start_lost、start_lost、stop_lost、incomplete_terminal_codon_variant&coding_sequence_variant（=start lost）
# missense_variant（错义变异）
# synonymous_variant、stop_retained_variant（等同于同义突变，极少出现）
more_severe_than_missense = ["stop_gain", "start_lost", "stop_lost", "incomplete_terminal"]
more_severe_than_syn = ["missense", "stop_gain", "start_lost", "stop_lost", "incomplete_terminal"]
