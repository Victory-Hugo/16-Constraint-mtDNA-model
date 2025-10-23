from composite_likelihood_mito import make_type_count_vector, probability_per_ref_nuc, lr_ref_nuc, \
    probability_per_class, lr_class
import csv
import datetime
from typing import Dict, List, TextIO, Union

#* ====================输入文件====================
#* /output/denovo/all_denovo.txt，包含所有 de novo 突变及其样本信息
#* ===============================================

#* ====================输出文件====================
#* /output/denovo/Ts_likelihood_by_category.txt(TSV)，记录阈值与四种转换的似然。
#* ===============================================

#* ====================输出文件格式================
# category	threshold	likelihood_C>T	likelihood_G>A	likelihood_T>C	likelihood_A>G
# germline	1	0.9856147735003647	6.948287620436667	3.9322413969875	2.6381160580030873
# germline + somatic tissue	1	1.1119951210060393	6.993886722124523	3.6526793270258486	2.0089057614263304
#* likelihood_* 是用观测到的 de novo 频率与参考区域碱基比例、突变类别概率等基线做对比得出的似然比。
#* 不同突变类型不直接可比（因为各自的观测数、参考碱基比例不同），但在同一类型内比较不同阈值/样本组合
#* 以 likelihood_C>T 这一列说明：
#* 1. 阈值 1 时：只看生殖系样本，似然 0.986；加入体细胞组织，升到 1.112；再把肿瘤样本也纳入，升到 1.195。说明在最低阈值下，体细胞组织和肿瘤样本比纯生殖系更富集 C>T 转换，尤其肿瘤贡献最大。
#* 2. 阈值 2 时：生殖系降到 0.964，加入体细胞组织变成 1.115，加入肿瘤后 1.191。趋势与阈值 1 类似，再次强调肿瘤数据会推高 C>T 的似然。
#* 3. 阈值 3、4、5、40 时：生殖系保持约 0.949，体细胞组织稳定在 ~1.12–1.20，肿瘤保持 >1.08，说明在放宽阈值后 C>T 估计已经收敛；不再受少数高负荷样本影响。
#* 总体，无论阈值如何选，C>T 在纯生殖系样本中略低于预期（<1），而一旦加入体细胞尤其肿瘤，就明显高于预期（>1），支持“肿瘤/体细胞样本推动 C>T 富集”的结论。
#* ===============================================

#* ====================为什么需要不同的阈值?========
#* 这里的 threshold 指“可接受的单个样本所含 de novo 突变最大数量”。
#* compare_denovo.py 会把同一类别里 de novo 数量超过阈值的样本整行剔除，认为它们可能是杂讯或极端值。
#* 之所以在每个类别上反复尝试 1、2、3、4、5、40 等阈值，是为了摸清模型对这些筛选条件的敏感度：
#* 如果阈值一放宽（比如从 1 到 5）似然值就剧烈波动，说明高负荷样本会显著影响突变率估计，可能需要谨慎处理。
#* 如果阈值变化时似然比较稳定，就能放心选择一个相对宽松的界限，在保留数据量的同时又不被异常样本主导。
#* ===============================================

#* ====================具体逻辑====================
#* 它帮助判断哪些样本类别或高变异样本会扭曲突变频谱，从而决定是否应该剔除。
#* mtDNA 数据里，偶尔会遇到单个样本记录了几十个乃至上百个 de novo 突变，这往往意味着污染、测序批次问题、或者肿瘤样本的强烈突变负荷
#* 如果把这些异常值照单全收，整个突变谱会被它们主导，模型算出的“平均”突变概率就不再代表大部分正常样本。
#* 按样本类别拆开处理 
#* germline 阈值参数 germline_max <- 生殖系 de novo 反映的是母系传递过程中真正传给子代的突变，最贴近“自然状态”。
#* somatic_tissue 阈值参数 som_tissue_max <- 体细胞组织里的 de novo 多数是个体一生中积累的突变，受组织类型、年龄等影响。。
#* somatic_cancer 阈值参数 som_cancer_max <- 体细胞癌症 de novo 反映的是肿瘤细胞中的突变，通常具有更高的突变负荷。
#* ===============================================

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
