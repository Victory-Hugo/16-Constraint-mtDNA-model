from composite_likelihood_mito import make_type_count_vector, probability_per_ref_nuc, lr_ref_nuc, \
    probability_per_class, lr_class
import csv
import datetime
from typing import Dict, List, TextIO, Union

'''该脚本用于比较不同类别的线粒体 DNA (mtDNA) de novo 突变，并计算 III 类突变（转换突变）的发生可能性。主要功能包括：

1. 定义突变类型、类别和分析区域（包括 OriB-OriH 区域和参考区域）。
2. 统计和过滤 de novo 变异，依据样本最大计数阈值筛选有效变异。
3. 计算 III 类突变的似然分数，结合参考碱基比例和突变类别概率。
4. 按不同阈值遍历各类 de novo 样本，输出 III 类突变的可能性到指定文件。

主要函数说明：
- build_denovo_dict: 生成包含 de novo 变异及其计数的字典。
- filter_denovo: 过滤 de novo 变异，移除样本 de novo 数量超过最大阈值的变异。
- compute_lr_classIII: 计算 III 类突变的发生可能性，返回各类型的似然分数。
- apply_threshold: 按不同阈值遍历 de novo 样本，计算并写入 III 类突变的可能性。

主流程：
- 依次遍历生殖系、体细胞组织、体细胞癌症三类样本的 de novo 最大计数阈值，计算并输出 III 类突变的似然分数。

依赖：
- composite_likelihood_mito.py 中的相关函数
- 输入文件：output/denovo/all_denovo.txt
- 输出文件：output/denovo/Ts_likelihood_by_category.txt

'''
# 全局变量
nucleotides = ["A", "C", "G", "T"]
class_I_mutations = ["C>A", "T>A", "G>T", "A>T"]
class_II_mutations = ["C>G", "T>G", "G>C", "A>C"]
class_III_mutations = ["C>T", "T>C", "G>A", "A>G"]
mut_types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]

# 构建线粒体突变模型时使用的区域定义
# ori 指代 OriB-OriH 区域，该区域在人工断点两侧（m.16197-191）具有已知的突变特征差异
start_ori = 16197
end_ori = 191 + 1  # range 的结束值不包含在范围内，因此需要 +1
ori_region = list(range(start_ori, 16570)) + list(range(1, end_ori))
reference_except_ori = list(range(end_ori, start_ori))


def build_denovo_dict(variant: str, dict: Dict[str, int]):
    """生成一个包含 de novo 变异及其计数的字典。

    :param variant: de novo 变异，采用 RefPosAlt 格式
    :param dict: 字典名称
    :return: 字典，键为变异，值为该变异的总计数
    """
    if variant not in dict:
        dict[variant] = 1
    else:
        dict[variant] += 1
    return dict


def filter_denovo(germline_max: Union[int, None], som_tissue_max: Union[int, None], som_cancer_max: Union[int, None],
                  denovo_counts: Dict[str, int]):
    """过滤 compile_denovo.py 生成的 all_denovo.txt 文件。
    通过移除样本 de novo 数量超过最大阈值的变异，得到最终保留的 de novo 清单。

    :param germline_max: 生殖系 de novo 的样本最大计数阈值
    :param som_tissue_max: 体细胞组织 de novo 的样本最大计数阈值
    :param som_cancer_max: 体细胞癌症 de novo 的样本最大计数阈值
    :param denovo_counts: 字典，键为变异，值为该变异的总计数
    :return: denovo_counts 字典
    """
    for row in csv.DictReader(open('output/denovo/all_denovo.txt'), delimiter='\t'):
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
    """计算 III 类突变的发生可能性。

    :param denovo_counts: 字典，键为变异，值为该变异的总计数
    :param reference_region: 待分析的 mtDNA 区域坐标列表
    :return: lambda_classIII 字典，键为 III 类突变，值为其可能性
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
    """针对某一类 de novo，按不同阈值遍历样本 de novo 最大计数。
    样本 de novo 数量超过阈值的变异将不会用于计算 III 类突变的可能性。

    :param germline_max: 生殖系 de novo 的样本最大计数阈值
    :param som_tissue_max: 体细胞组织 de novo 的样本最大计数阈值
    :param som_cancer_max: 体细胞癌症 de novo 的样本最大计数阈值
    :param category: 当前迭代的 de novo 类别
    :param threshold: 当前遍历的样本 de novo 最大计数，用于写入文件
    :param reference_region: 待分析的 mtDNA 区域坐标列表
    :param file: 用于写入 III 类突变可能性的文本文件对象
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

    f = open('output/denovo/Ts_likelihood_by_category.txt', "w")
    f.write("category	threshold	likelihood_C>T	likelihood_G>A	likelihood_T>C	likelihood_A>G" + '\n')

    # 作为样本 de novo 最大计数的候选阈值集合
    threshold_range = [1, 2, 3, 4, 5, 40]  # 已知最大值为 38，因此手动挑选需要测试的阈值

    for threshold in threshold_range:
        # 首先仅遍历生殖系样本
        apply_threshold(germline_max=threshold, som_tissue_max=0, som_cancer_max=0,
                        category="germline", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # 然后保留所有生殖系，并遍历体细胞组织
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=threshold, som_cancer_max=0,
                        category="germline + somatic tissue", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # 最后在前两类基础上加入体细胞癌症样本
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=max(threshold_range),
                        som_cancer_max=threshold, category="germline + somatic tissue + somatic cancer",
                        threshold=threshold, reference_region=reference_except_ori, file=f)

    print(datetime.datetime.now(), "Script finished!")
