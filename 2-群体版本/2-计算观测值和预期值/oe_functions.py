import csv
import numpy as np
from scipy.stats import chi2, poisson
from typing import Dict, List, TextIO, Tuple, Union

'''
该模块以携带者计数而非异质性为核心统计量，负责计算线粒体 DNA 变异的观测值、预期值、观测/预期比值及其置信区间，用于评估约束程度。主要功能包括：
- initialize_sum_dict: 初始化用于统计不同变异类别的携带者总数、似然得分和计数的字典。
- sum_obs_likelihood: 累加某类变异的携带者总数、似然加权暴露量和计数。
- linear_model: 使用针对中性变异拟合的线性模型，将暴露量转换为预期携带者计数。
- calibrate: 根据线性模型参数，返回某类变异的预期携带者计数之和。
- calculate_exp: 读取线性模型参数文件，计算预期携带者计数。
- calculate_obs: 计算某类变异的观测携带者计数之和。
- calculate_total: 统计可能的变异总数。
- calculate_CI: 采用 Poisson 近似计算观测/预期比值的 90% 置信区间。
- calculate_pvalue: 基于 Poisson 分布估计观测值显著低于预期值的概率。
- calculate_oe: 汇总观测值、预期值和置信区间，可写入文件或作为结果返回。
模块还定义了 OriB-OriH 区域的坐标以及蛋白编码基因后果的严重性分级列表。
'''
def initialize_sum_dict(identifier_list: List[Union[str, Tuple[int, int]]]):
	"""初始化一个字典，用于累加不同类别变异的携带者计数、似然得分和计数。

	:param identifier_list: 需要求和的每个变异类别的标识符
	:return: 一个以标识符、突变分组、需要求和的值以及包含变异的区域组成的元组为键、以0为初始值的字典
	"""
	# 初始化字典，使所有值为0
	dict = {}
	for identifier in identifier_list:
		for mut_group in ['G>A_and_T>C', 'other']:  # 要拟合的两个突变分组
			for value in ['observed_carriers', 'sum_LR_AN', 'sum_callable', 'count']:  # 需要求和的数值
				for region in ['ref_exc_ori', 'ori']:  # 可能应用的两种突变模型
					dict[identifier, mut_group, value, region] = 0
	return dict


def sum_obs_likelihood(
		mutation: str, identifier: Union[str, Tuple[int, int]], region: str,
		observed: Union[str, float], likelihood: Union[str, float], callable_samples: Union[str, float, int],
		dict: Dict[Tuple[str, str, str, str], Union[float, int]]):
	"""对某一类变异的携带者计数、突变似然得分和计数进行求和。

	:param mutation: Ref>Alt 形式的突变类型
	:param identifier: 需要求和的变异类别标识符
	:param region: 取值应为 'ref_exc_ori' 或 'ori'，用于指示应用哪种突变模型
	:param observed: 观测携带者计数
	:param likelihood: 用于突变似然得分的取值
	:param callable_samples: 当前变异可判读的样本数量
	:param dict: 需要更新的字典名称
	:return: 一个以标识符、突变分组、需要求和的值和变异所在区域组成的元组为键、以累计值为值的字典
	"""
	# 高度可变的 G>A 和 T>C 单独拟合
	mut_group = 'G>A_and_T>C' if (mutation == 'G>A' or mutation == 'T>C') else 'other'
	# 注意字典应已初始化（所有预期键的初始值均为0）
	callable_value = float(callable_samples) if callable_samples is not None else 0.0
	dict[(identifier, mut_group, 'observed_carriers', region)] += float(observed)
	dict[(identifier, mut_group, 'sum_LR_AN', region)] += float(likelihood) * callable_value
	dict[(identifier, mut_group, 'sum_callable', region)] += callable_value
	dict[(identifier, mut_group, 'count', region)] += 1
	return dict


def linear_model(intercept: float, coefficient: float, exposure: float):
	"""将针对中性变异拟合的线性模型应用于数据，以计算预期的携带者计数。

	:param intercept: 线性方程的截距
	:param coefficient: 线性方程的系数
	:param exposure: 变异的加权突变倾向（似然 * 可判读样本数）的总和
	:return: 预期的携带者计数
	"""
	return max(float(0), intercept + (coefficient * exposure))


def calibrate(
		sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]], mut_group: str, identifier: str, region: str,
		parameters: Dict[Tuple[str, str, str], float]):
	"""应用线性模型，返回预期的携带者计数之和。

	:param sum_dict: 包含计算所需的变异数量和突变得分总和的字典
	:param mut_group: 需要校准的突变分组，取值为 'G>A_and_T>C' 或 'other'
	:param identifier: 需要校准的变异类别标识符
	:param region: 取值应为 'ref_exc_ori' 或 'ori'，用于指示使用哪种校准
	:param parameters: 包含线性方程系数和截距的字典，用于执行校准
	:return: 线性模型的输出，即相应突变分组和区域的预期携带者计数之和
	"""
	return linear_model(
		intercept=parameters[(mut_group, region, 'intercept')],
		coefficient=parameters[(mut_group, region, 'coefficient')],
		exposure=sum_dict[(identifier, mut_group, 'sum_LR_AN', region)])


def calculate_exp(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]], fit_parameters: str):
	"""计算某一类变异的预期携带者计数之和。

	:param sum_dict: 包含计算所需的变异数量和突变得分总和的字典
	:param identifier: 需要计算预期值的变异类别标识符
	:param fit_parameters: 包含线性方程系数和截距的文件路径
	:return: 预期的携带者计数之和
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
	"""计算某一类变异的观测携带者计数之和。

	:param identifier: 需要计算观测值的变异类别标识符
	:param sum_dict: 包含观测携带者计数的字典
	:return: 跨越两类突变分组和两类区域类型的观测携带者计数之和
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'observed_carriers', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'observed_carriers', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'observed_carriers', 'ori')] + sum_dict[
		(identifier, 'other', 'observed_carriers', 'ori')]


def calculate_total(
		identifier: Union[str, Tuple[int, int]],
		sum_dict: Dict[Tuple[Union[str, Tuple[int, int]], str, str, str], Union[float, int]]):
	"""计算某一组变异的可能变异总数。

	:param identifier: 需要计算的变异类别标识符
	:param sum_dict: 包含观测携带者计数的字典
	:return: 跨越两类突变分组和两类区域类型的可能变异总数
	"""
	return sum_dict[(identifier, 'G>A_and_T>C', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'other', 'count', 'ref_exc_ori')] + sum_dict[
		(identifier, 'G>A_and_T>C', 'count', 'ori')] + sum_dict[
		(identifier, 'other', 'count', 'ori')]


def calculate_CI(observed_carriers: float, total: int, expected_carriers: float, alpha: float = 0.1):
	"""使用 Poisson 近似计算观测值与预期值之比的90%置信区间。

	:param observed_carriers: 观测携带者数之和
	:param total: 变异数量（保留参数以保持函数签名兼容）
	:param expected_carriers: 预期携带者数之和
	:param alpha: 显著性水平，默认 0.1 对应 90% 置信区间
	:return: lower_CI 和 upper_CI，即90%置信区间的上下界
	"""
	if expected_carriers <= 0:
		return 0.0, 0.0

	obs_int = int(round(observed_carriers))
	lower = 0.0 if obs_int == 0 else chi2.ppf(alpha / 2, 2 * obs_int) / (2 * expected_carriers)
	upper = chi2.ppf(1 - alpha / 2, 2 * (obs_int + 1)) / (2 * expected_carriers)
	return lower, upper


def calculate_pvalue(observed_carriers: float, total: int, expected_carriers: float):
	"""基于 Poisson 分布检验观测值是否显著低于预期值。

	:param observed_carriers: 观测携带者数之和
	:param total: 变异数量（保留参数以保持函数签名兼容）
	:param expected_carriers: 预期携带者数之和
	:return: less_than_pvalue，即在给定预期值的情况下，观测到小于或等于该观测值的概率
	"""
	if expected_carriers <= 0:
		return 1.0
	obs_int = int(round(observed_carriers))
	return poisson.cdf(obs_int, mu=expected_carriers)


def calculate_oe(
		item: Union[str, Tuple[int, int]], sum_dict: Dict[Tuple[str, str, str, str], Union[float, int]],
		fit_parameters: str, file: Union[TextIO, None], output: str = "write", alpha: float = 0.1):
	"""计算观测值、预期值、观测与预期之比以及90%置信区间。

	:param item: 需要计算预期值的变异类别标识符
	:param sum_dict: 包含计算所需的可能变异数量和突变得分总和的字典
	:param fit_parameters: 包含线性方程系数与截距的文件路径
	:param file: 当 output 为 "write" 时用于写入结果的文件
	:param output: 指定输出到文件或保存为字典，默认 "write"，可选 "dict"
	:param alpha: 置信区间的显著性水平，默认 0.1（即 90% 置信区间）
	"""
	expected_carriers = calculate_exp(sum_dict=sum_dict, identifier=item, fit_parameters=fit_parameters)
	observed_carriers = calculate_obs(sum_dict=sum_dict, identifier=item)
	ratio_oe = (observed_carriers / expected_carriers) if expected_carriers > 0 else 0.0
	total_all = calculate_total(sum_dict=sum_dict, identifier=item)

	(lower_CI, upper_CI) = calculate_CI(
		observed_carriers=observed_carriers, total=total_all, expected_carriers=expected_carriers, alpha=alpha)
	
	print("Calculating values for", item)
	if output == "write":
		file.write(
			item + '\t' + str(total_all) + '\t' + str(observed_carriers) + '\t' + str(expected_carriers)
			+ '\t' + str(ratio_oe) + '\t' + str(lower_CI) + '\t' + str(upper_CI) + '\n')
	elif output == "dict":
		return total_all, observed_carriers, expected_carriers, ratio_oe, lower_CI, upper_CI


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
