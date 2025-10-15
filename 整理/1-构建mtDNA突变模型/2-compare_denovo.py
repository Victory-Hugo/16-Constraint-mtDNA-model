import importlib.util
import sys

# 导入带有破折号的模块
spec = importlib.util.spec_from_file_location("composite_likelihood_mito", "4-composite_likelihood_mito.py")
composite_likelihood_mito = importlib.util.module_from_spec(spec)
sys.modules["composite_likelihood_mito"] = composite_likelihood_mito
spec.loader.exec_module(composite_likelihood_mito)

from composite_likelihood_mito import make_type_count_vector, probability_per_ref_nuc, lr_ref_nuc, \
    probability_per_class, lr_class
import csv
import datetime
from typing import Dict, List, TextIO, Union


'''
比较不同类型的线粒体 de novo 突变（生殖系、体细胞组织、体细胞癌症）
在不同样本突变负荷阈值下的突变似然（尤其是 III 类突变，即转换突变）。
程序读取由前一个脚本 compile_denovo.py 生成的整合文件。
通过 filter_denovo() 函数，只保留样本的 de novo 数量小于或等于指定阈值的突变。
对三类数据独立设定最大计数阈值：
    germline_max：生殖系；
    som_tissue_max：体细胞组织；
    som_cancer_max：体细胞癌症。    
计算 III 类突变（转换突变）的似然值
使用外部模块 composite_likelihood_mito 中的函数计算每类突变的概率：
    make_type_count_vector()：统计各突变类型；
    probability_per_ref_nuc()：计算参考碱基比例；
    lr_ref_nuc()：计算各碱基的突变似然；
    lr_class() 和 probability_per_class()：计算突变类别的似然。
程序循环不同阈值（[1, 2, 3, 4, 5, 40]），分别测试：
    仅生殖系；
    生殖系 + 体细胞组织；
    生殖系 + 体细胞组织 + 体细胞癌症；
    每次调用 apply_threshold() 计算对应条件下的似然结果。
'''
# global variables
nucleotides = ["A", "C", "G", "T"]
class_I_mutations = ["C>A", "T>A", "G>T", "A>T"]
class_II_mutations = ["C>G", "T>G", "G>C", "A>C"]
class_III_mutations = ["C>T", "T>C", "G>A", "A>G"]
mut_types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]

# 构建线粒体突变模型时使用的区域定义
# ori 指的是在人工断点 m.16197-191 两端具有已知突变特征差异的 OriB-OriH 区域
start_ori = 16197
end_ori = 191 + 1  # end not included in range so + 1
ori_region = list(range(start_ori, 16570)) + list(range(1, end_ori))
reference_except_ori = list(range(end_ori, start_ori))


def build_denovo_dict(variant: str, dict: Dict[str, int]):
    """生成包含新生变异及其总计数的字典。

    :param variant: 新生变异，采用 RefPosAlt 格式
    :param dict: 字典名称
    :return: 返回一个字典，变异为键，值为该变异的出现总次数
    """
    if variant not in dict:
        dict[variant] = 1
    else:
        dict[variant] += 1
    return dict


def filter_denovo(germline_max: Union[int, None], som_tissue_max: Union[int, None], som_cancer_max: Union[int, None],
                  denovo_counts: Dict[str, int]):
    """过滤由 compile_denovo.py 生成的 all_denovo.txt 文件。
    目的是移除来自样本中超过最大新生变异计数阈值的变异。

    :param germline_max: 用于筛选的生殖系（germline）样本的最大新生变异数
    :param som_tissue_max: 用于筛选的体细胞组织（somatic tissue）样本的最大新生变异数
    :param som_cancer_max: 用于筛选的体细胞癌症（somatic cancer）样本的最大新生变异数
    :param denovo_counts: 一个字典，变异为键，值为该变异的总计数
    :return: 返回更新后的 denovo_counts 字典
    """
    for row in csv.DictReader(open('output_files/denovo/all_denovo.txt'), delimiter='\t'):
        if "germline" in row["sample"]:
            if germline_max is not None:
                if int(row["sample_denovo_count"]) <= germline_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
        elif "somatic_tissue" in row["sample"]:
            if som_tissue_max is not None:
                if int(row["sample_denovo_count"]) <= som_tissue_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
        elif "somatic_cancer" in row["sample"]:
            if som_cancer_max is not None:
                if int(row["sample_denovo_count"]) <= som_cancer_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
    return denovo_counts


def compute_lr_classIII(denovo_counts: Dict[str, int], reference_region: List[int]):
    """计算 III 类（class III）突变的突变似然（likelihood）。

    :param denovo_counts: 一个字典，变异为键，值为该变异的总计数
    :param reference_region: 要分析的线粒体 DNA 区域坐标列表
    :return: 返回 lambda_classIII 字典，键为 III 类突变，值为对应的似然值
    """
    mut_type_counts = make_type_count_vector(denovo_counts=denovo_counts, reference_region=reference_region)
    base_proportions = probability_per_ref_nuc(reference_region=reference_region)
    lambda_ref_nuc = lr_ref_nuc(type_counts=mut_type_counts, base_proportions=base_proportions)
    class_probabilities = probability_per_class(base_proportions=base_proportions)
    lambda_mut_class = lr_class(type_counts=mut_type_counts, class_probabilities=class_probabilities)

    lambda_classIII = {}
    for mut in class_III_mutations:
        lambda_classIII[mut] = lambda_ref_nuc[mut[0]] * lambda_mut_class[mut]
    return lambda_classIII


def apply_threshold(germline_max: int, som_tissue_max: int, som_cancer_max: int, category: str, threshold: int,
                    reference_region: List[int], file: TextIO):
    """对于某一类新生变异（de novo），在一系列阈值上迭代最大样本新生变异数。
    来自样本中新生变异数超过阈值的变异不会被包含在 III 类似然值的计算中。

    :param germline_max: 用于筛选的生殖系（germline）样本的最大新生变异数
    :param som_tissue_max: 用于筛选的体细胞组织（somatic tissue）样本的最大新生变异数
    :param som_cancer_max: 用于筛选的体细胞癌症（somatic cancer）样本的最大新生变异数
    :param category: 当前正在迭代阈值的变异类别
    :param threshold: 本次循环中用于打印的最大样本新生变异数阈值
    :param reference_region: 要分析的线粒体 DNA 区域坐标列表
    :param file: 用于写入 III 类似然值的文本文件对象
    """
    denovo_counts = {}
    denovo_counts = filter_denovo(germline_max=germline_max, som_tissue_max=som_tissue_max,
                                  som_cancer_max=som_cancer_max, denovo_counts=denovo_counts)

    print('\n' + "For ", category, " testing maximum sample de novo count of ", threshold, '\n')
    lambda_classIII = compute_lr_classIII(denovo_counts=denovo_counts, reference_region=reference_region)

    file.write(str(category) + '\t' + str(threshold) + '\t' +
               str(lambda_classIII["C>T"]) + '\t' + str(lambda_classIII["G>A"]) + '\t' +
               str(lambda_classIII["T>C"]) + '\t' + str(lambda_classIII["A>G"]) + '\n')


if __name__ == "__main__":
    print(datetime.datetime.now(), "Starting to compare de novo categories!")
    print('\n' + "This will produce likelihood scores for transitions across sample categories")

    f = open('output_files/denovo/Ts_likelihood_by_category.txt', "w")
    f.write("category	threshold	likelihood_C>T	likelihood_G>A	likelihood_T>C	likelihood_A>G" + '\n')

    # 用作样本 de novo 最大计数阈值的一组取值范围
    threshold_range = [1, 2, 3, 4, 5, 40]  # 已知最大值为 38，因此手动设定需要测试的阈值

    for threshold in threshold_range:
        # first, iterate through germline only
        apply_threshold(germline_max=threshold, som_tissue_max=0, som_cancer_max=0,
                        category="germline", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # then, include all germline plus iterate through somatic tissue
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=threshold, som_cancer_max=0,
                        category="germline + somatic tissue", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # then, include all germline and somatic tissue, plus iterate through somatic cancer
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=max(threshold_range),
                        som_cancer_max=threshold, category="germline + somatic tissue + somatic cancer",
                        threshold=threshold, reference_region=reference_except_ori, file=f)

    print(datetime.datetime.now(), "Script finished!")
