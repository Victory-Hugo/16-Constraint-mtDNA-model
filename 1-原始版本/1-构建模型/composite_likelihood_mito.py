import argparse
from compile_denovo import rcrs_pos_to_ref
import csv
import datetime
from typing import Dict, List, Tuple, Union
import os
'''
composite_likelihood_mito.py

该脚本用于构建线粒体DNA（mtDNA）突变的复合似然模型，结合参考碱基、突变类别和序列上下文信息，
计算全基因组范围内每个位点、每种突变类型的突变概率。主要功能包括：

1. 数据读取与预处理
    - 读取de novo突变列表（RefPosAlt格式，附计数）。
    - 读取线粒体参考基因组FASTA文件和所有可能的SNV列表（VCF格式）。

2. 区域定义
    - 支持对OriB-OriH区域（m.16197-191）与其他参考区域分别建模，考虑人工断点和特殊间隔（如m.3107N）。

3. 统计与概率计算
    - 统计各突变类型在指定区域的出现次数。
    - 计算参考碱基（A、C、G、T）在区域内的比例。
    - 计算各突变类别（I、II、III）的概率与似然比。
    - 统计每个位点、每种突变类型的出现次数。

4. 序列上下文建模
    - 统计突变及参考碱基在侧翼窗口（如三核苷酸）内的序列上下文分布。
    - 分别处理转换（Ts）和颠换（Tv）突变，支持跨链合并。

5. 复合似然计算
    - 综合参考碱基、突变类别和序列上下文的似然比，计算每个位点、每种突变类型的复合突变概率。
    - OriB-OriH区域的转换突变单独处理。

6. 文件输出
    - 输出各类统计结果和似然比至指定文件夹，便于后续绘图和分析。

主要函数说明：
- build_dictionary, build_additive_dictionary: 辅助字典构建与计数。
- flip_base, flip_mut: 计算反义链碱基及突变类型。
- handle_mito_genome: 校正环状mtDNA及特殊间隔的侧翼坐标。
- make_denovo_counts, make_type_count_vector: 读取并统计de novo突变数据。
- probability_per_ref_nuc, lr_ref_nuc: 计算参考碱基比例及突变似然比。
- probability_per_class, lr_class: 计算突变类别概率及似然比。
- count_type_per_pos, make_mut_context_vector: 统计每个位点、每种突变类型的上下文分布。
- make_ref_freq_vector: 统计参考序列上下文分布。
- lr_seq_context: 计算序列上下文的似然比。
- composite_likelihood: 综合各类似然比，计算全基因组突变概率。
- write_*: 各类统计结果写入文件。


输出：
- 各类统计文件及最终的mito_mutation_likelihoods.txt，包含每个位点、每种突变类型的复合似然。

适用场景：
- 线粒体突变谱建模、序列上下文分析、突变概率预测、可视化等科研任务。
'''
# 全局变量
nucleotides = ["A", "C", "G", "T"]
class_I_mutations = ["C>A", "T>A", "G>T", "A>T"]
class_II_mutations = ["C>G", "T>G", "G>C", "A>C"]
class_III_mutations = ["C>T", "T>C", "G>A", "A>G"]
mut_types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]

# 构建线粒体突变模型时使用的区域定义
# ori 指代 OriB-OriH 区域，在人工断点两侧（m.16197-191）具有已知的突变特征差异
start_ori = 16197
end_ori = 191 + 1  # 因 range 的右端不包含该值，因此需要加 1
ori_region = list(range(start_ori, 16570)) + list(range(1, end_ori))
reference_except_ori = list(range(end_ori, start_ori))


# 辅助函数
def build_dictionary(key: Union[str, Tuple[str, str, int]], dictionary: Dict[str, int]):
    """生成一个字典，为给定键记录出现次数。

    :param key: 待统计的字符串或元组
    :param dictionary: 字典名称
    :return: 字典，键为待统计对象，值为出现次数
    """
    if key not in dictionary:
        dictionary[key] = 1
    else:
        dictionary[key] += 1
    return dictionary


def build_additive_dictionary(key: Union[str, Tuple[str, str, int]], value: int, dictionary: Dict[str, int]):
    """生成一个字典，为给定键累计计数，确保每个键的值均为非零。

    :param key: 需要累加的字符串或元组
    :param value: 要累加的计数或数值
    :param dictionary: 字典名称
    :return: 字典，键为累加对象，值为累计计数
    """
    if key not in dictionary:
        dictionary[key] = value
    else:
        dictionary[key] += value
    return dictionary


def flip_base(base: str):
    """返回同一位置反义链上的配对碱基。

    :param base: 参考碱基
    :return: flipped_base，配对碱基
    """
    flipped_base = ''
    if base == "A":
        flipped_base = "T"
    elif base == "C":
        flipped_base = "G"
    elif base == "G":
        flipped_base = "C"
    elif base == "T":
        flipped_base = "A"
    return flipped_base


def flip_mut(mut: str):
    """返回反义链同一位置上等效的突变类型。

    :param mut: Ref>Alt 格式的突变类型
    :return: flipped_mut，反义链上的等效突变类型
    """
    flipped_mut = ''
    if mut == "A>C":
        flipped_mut = "T>G"
    elif mut == "A>G":
        flipped_mut = "T>C"
    elif mut == "A>T":
        flipped_mut = "T>A"
    elif mut == "G>A":
        flipped_mut = "C>T"
    elif mut == "G>C":
        flipped_mut = "C>G"
    elif mut == "G>T":
        flipped_mut = "C>A"
    elif mut == "C>A":
        flipped_mut = "G>T"
    elif mut == "C>G":
        flipped_mut = "G>C"
    elif mut == "C>T":
        flipped_mut = "G>A"
    elif mut == "T>A":
        flipped_mut = "A>T"
    elif mut == "T>C":
        flipped_mut = "A>G"
    elif mut == "T>G":
        flipped_mut = "A>C"
    return flipped_mut


def handle_mito_genome(flanking_coord: int, coord: int):
    """针对环状 mtDNA 上人工断点附近或 m.3107N 间隔附近的侧翼碱基，返回正确的侧翼坐标。
    可处理 [-10:10]/{0} 范围内侧翼碱基的校正。

    :param flanking_coord: 侧翼碱基的原始坐标
    :param coord: 参考碱基的坐标
    :return: 纠正后的侧翼碱基坐标字符串，如无需校正则保持不变
    """
    if flanking_coord < 1:  # 处理环状基因组
        flanking_coord = 16569 + flanking_coord  # 因为 flanking_coord 此时为 0 或负数
    elif flanking_coord > 16569:  # 处理环状基因组
        flanking_coord = flanking_coord - 16569
    elif flanking_coord < 3107 and 3118 > coord > 3107:  # 处理 m.3107N，跳过该位置
        flanking_coord = flanking_coord - 1
    elif flanking_coord > 3107 and 3096 < coord < 3107:  # 处理 m.3107N，跳过该位置
        flanking_coord = flanking_coord + 1
    elif flanking_coord == 3107 and flanking_coord < coord:  # 处理 m.3107N
        flanking_coord = 3106
    elif flanking_coord == 3107 and flanking_coord > coord:  # 处理 m.3107N
        flanking_coord = 3108
    flanking_coord = str(flanking_coord)
    return flanking_coord


def write_file_header(file: str, flanking_range: List[int], variable: str):
    """创建文件及其表头，用于写入序列上下文数据以便绘图。

    :param file: 输出文件路径
    :param flanking_range: 需要写入的侧翼位置列表，范围为 [-10:10]/{0}
    :param variable: 用于与侧翼碱基比较的数据变量名
    :return: f，已打开的输出文件对象
    """
    f = open(file, "w")
    f.write("flanking_nucleotide" + '\t')
    for flanking_pos in flanking_range:
        if flanking_pos == flanking_range[-1]:
            f.write(str(flanking_pos) + '\t' + variable + '\n')
        else:
            f.write(str(flanking_pos) + '\t')
    return f


# 复合似然模型相关函数
def make_denovo_counts(denovo_list: str):
    """读取最终的 de novo 变异列表并转换为字典。

    :param denovo_list: 最终使用的 de novo 变异列表路径，RefPosAlt 格式，附带计数
    :return: denovo_counts 字典，键为变异，值为该变异的总计数
    """
    denovo_counts = {}
    for row in csv.DictReader(open(denovo_list), delimiter='\t'):
        denovo_counts[row["denovo"]] = int(row["count"])
    return denovo_counts


def make_type_count_vector(denovo_counts: Dict[str, int], reference_region: List[int]):
    """统计指定区域内各突变类型的数量。

    :param denovo_counts: 字典，键为变异，值为该变异的总计数
    :param reference_region: mtDNA 区域的坐标列表，用于统计各突变类型的出现次数
    :return: mut_type_counts 字典，键为突变类型，值为该区域内的总计数
    """
    mut_type_counts = {}
    # 使用所有可能的 SNV 列表以加速解析
    for row in csv.DictReader(open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        variant = row["REF"] + row["POS"] + row["ALT"]
        mut_type = row["REF"] + ">" + row["ALT"]
        coord = row["POS"]
        # 统计 12 种突变的出现次数
        if variant in denovo_counts:
            if int(coord) in reference_region:
                allele_count = denovo_counts[variant]
                mut_type_counts = build_additive_dictionary(key=mut_type, value=allele_count,
                                                            dictionary=mut_type_counts)
    print('\n' + "Number of observed de novo mutations: ", sum(mut_type_counts.values()))
    print('\n' + "Mutation type counts for region: ", mut_type_counts)
    return mut_type_counts


def probability_per_ref_nuc(reference_region: List[int]):
    """统计 mtDNA 区域内各碱基（A、C、G、T）的出现次数，并换算成比例。
    该比例用于表示每个参考碱基发生突变的概率。

    :param reference_region: mtDNA 区域的坐标列表，用于统计各碱基出现次数
    :return: base_proportions 字典，键为碱基，值为其在该区域内的占比
    """
    # 使用线粒体参考基因组的 FASTA 文件，去除换行以便逐位解析
    fasta = open('0-required_files/synthetic_vcf/NC_012920.1_noheader.fasta').read().replace('\n', '')

    dictionary = {}
    # 构建字典以统计区域内 A、C、G、T 的数量
    for coord in reference_region:
        # FASTA 下标从 0 开始，因此需要减 1；若为 N 则跳过 m.3107
        if fasta[coord - 1] != "N":
            dictionary = build_dictionary(key=fasta[coord - 1], dictionary=dictionary)
    total_number_bases = sum(dictionary.values())

    # 计算比例，键为 A、C、G、T
    base_proportions = {}
    for key in dictionary:
        base_proportions[key] = dictionary[key] / total_number_bases
    print('\n' + "The proportion of each base is: ", base_proportions)
    return base_proportions


def lr_ref_nuc(type_counts: Dict[str, int], base_proportions: Dict[str, float]):
    """计算每个参考碱基 n(t) = A、C、G、T 的突变似然比。

    :param type_counts: 字典，键为突变类型，值为对应计数
    :param base_proportions: 字典，键为碱基，值为其占比
    :return: lambda_ref_nuc 字典，键为参考碱基，值为突变似然比
    """
    sum_vtype = sum(type_counts.values())
    lambda_ref_nuc = {}
    # 对每个碱基，统计以其为参考碱基的突变总数，再计算似然比
    for base in nucleotides:
        total_mut_at_base = 0
        for key in type_counts:
            if key.startswith(base + ">"):
                total_mut_at_base += type_counts[key]
        lambda_ref_nuc[base] = total_mut_at_base / (sum_vtype * base_proportions[base])
    print('\n' + "The likelihood of mutation at each reference nucleotide is: ", lambda_ref_nuc)
    return lambda_ref_nuc


def probability_per_class(base_proportions: Dict[str, float]):
    """利用替换碱基的频率信息，计算每个突变类别的概率。

    :param base_proportions: 字典，键为碱基，值为其占比
    :return: class_probabilities 字典，键为突变类别，值为对应的发生概率
    """
    class_probabilities = {}
    probC_to_A = base_proportions["C"] * (base_proportions["A"] / (1 - base_proportions["C"]))  # 类型 1
    probG_to_T = base_proportions["G"] * (base_proportions["T"] / (1 - base_proportions["G"]))  # 类型 7
    probT_to_A = base_proportions["T"] * (base_proportions["A"] / (1 - base_proportions["T"]))  # 类型 4
    probA_to_T = base_proportions["A"] * (base_proportions["T"] / (1 - base_proportions["A"]))  # 类型 10
    class_probabilities["I"] = probC_to_A + probG_to_T + probT_to_A + probA_to_T

    probC_to_G = base_proportions["C"] * (base_proportions["G"] / (1 - base_proportions["C"]))  # 类型 2
    probG_to_C = base_proportions["G"] * (base_proportions["C"] / (1 - base_proportions["G"]))  # 类型 8
    probT_to_G = base_proportions["T"] * (base_proportions["G"] / (1 - base_proportions["T"]))  # 类型 6
    probA_to_C = base_proportions["A"] * (base_proportions["C"] / (1 - base_proportions["A"]))  # 类型 12
    class_probabilities["II"] = probC_to_G + probG_to_C + probT_to_G + probA_to_C

    probC_to_T = base_proportions["C"] * (base_proportions["T"] / (1 - base_proportions["C"]))  # 类型 3
    probG_to_A = base_proportions["G"] * (base_proportions["A"] / (1 - base_proportions["G"]))  # 类型 9
    probT_to_C = base_proportions["T"] * (base_proportions["C"] / (1 - base_proportions["T"]))  # 类型 5
    probA_to_G = base_proportions["A"] * (base_proportions["G"] / (1 - base_proportions["A"]))  # 类型 11
    class_probabilities["III"] = probC_to_T + probG_to_A + probT_to_C + probA_to_G

    print('\n' + "The probability of mutation classes I, II and III is: ", class_probabilities)
    return class_probabilities


def lr_class(type_counts: Dict[str, int], class_probabilities: Dict[str, float]):
    """计算各突变类别（I、II、III）的似然比。

    :param type_counts: 字典，键为突变类型，值为计数
    :param class_probabilities: 字典，键为突变类别，值为对应概率
    :return: lambda_mut_class 字典，键为突变类型，值为突变类别的似然比
    """
    # 统计每个突变类别的突变数量
    total_I = total_II = total_III = 0
    for key in type_counts:
        if key in class_I_mutations:
            total_I += type_counts[key]
        elif key in class_II_mutations:
            total_II += type_counts[key]
        elif key in class_III_mutations:
            total_III += type_counts[key]

    # 计算各类别的似然比
    sum_vtype = sum(type_counts.values())
    lambda_mut_class = {}
    for mut in mut_types:
        if mut in class_I_mutations:
            lambda_mut_class[mut] = total_I / (sum_vtype * class_probabilities["I"])
        elif mut in class_II_mutations:
            lambda_mut_class[mut] = total_II / (sum_vtype * class_probabilities["II"])
        elif mut in class_III_mutations:
            lambda_mut_class[mut] = total_III / (sum_vtype * class_probabilities["III"])
    print('\n' + "The likelihood of each mutation class is: ", lambda_mut_class)
    return lambda_mut_class


def count_type_per_pos(denovo_counts: Dict[str, int]):
    """统计每个位点上各突变类型的出现次数，为 make_mut_context_vector 提供输入。
    由于键中包含位点信息，返回的字典可用于任意区域。

    :param denovo_counts: 字典，键为变异，值为该变异的总计数
    :return: pos_by_type_counts 字典，记录每个位点 12 种突变类型的计数
    """
    pos_by_type_counts = {}
    # 使用所有可能的 SNV 列表以加速解析
    for row in csv.DictReader(open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        variant = row["REF"] + row["POS"] + row["ALT"]
        mut_type = row["REF"] + ">" + row["ALT"]
        coord = row["POS"]
        # 统计每个位点 12 种突变的次数
        if variant in denovo_counts:
            allele_count = denovo_counts[variant]
            pos_by_type_counts[(coord, mut_type)] = allele_count
    return pos_by_type_counts


def make_mut_context_vector(pos_by_type_counts: Dict[Tuple[str, str], int],
                            flanking_range: List[int], reference_region: List[int], mut_group: str = None):
    """生成突变序列上下文向量 Vseq t,p,n，统计每种突变在各侧翼位置的侧翼碱基数量，可同时处理转换和颠换。

    对四种转换分别计算序列上下文，并考虑参考链的方向。对于颠换，由于计数较少，按嘧啶突变（C 或 T）汇总处理。

    :param pos_by_type_counts: 字典，记录每个位点 12 种突变类型的计数
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param reference_region: 待计算的 mtDNA 区域坐标列表
    :param mut_group: 可选值 Ts 或 Tv，分别代表转换或颠换，默认同时包含二者
    :return: vseq_vector 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为对应计数；对 Tv 来说，突变键为嘧啶突变（C>A、C>G、T>A、T>G）
    """
    # 首先生成坐标到参考碱基的字典
    rcrs_pos2ref = rcrs_pos_to_ref()
    vseq_vector = {}
    for mut in mut_types:
        for coord in reference_region:
            if (str(coord), mut) in pos_by_type_counts:  # 若该坐标存在突变
                for flanking_pos in flanking_range:
                    if mut in class_III_mutations:  # 转换（Ts）
                        if mut_group == "Tv":
                            continue  # 如果只分析颠换，则跳过
                        # 统计每种突变周围侧翼位置的 A、C、G、T 数量
                        # 先将 [-10:10]/{0} 范围内的侧翼位置转换为基因组坐标
                        # 再将侧翼坐标映射为参考基因组中的侧翼碱基
                        # 将每个侧翼位置的侧翼碱基计数保存下来
                        flanking_coord = coord + flanking_pos
                        # 该函数可处理人工断点和 m.3107N 间隔附近的坐标
                        flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                        flanking_nuc = rcrs_pos2ref[flanking_coord]  # 侧翼位置的参考碱基
                        vseq_vector = build_additive_dictionary(key=(mut, flanking_nuc, flanking_pos),
                                                                value=pos_by_type_counts[(str(coord), mut)],
                                                                dictionary=vseq_vector)
                    else:  # 颠换（Tv）
                        if mut_group == "Ts":
                            continue  # 如果只分析转换，则跳过
                        # 由于计数较少，对 Tv 进行跨链合并
                        # 与上述逻辑类似，但以嘧啶（C 或 T）为参考记录突变
                        if rcrs_pos2ref[str(coord)] == "C" or rcrs_pos2ref[str(coord)] == "T":
                            pyr_mut = mut
                            flanking_coord = coord + flanking_pos
                            flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                            flanking_nuc = rcrs_pos2ref[flanking_coord]
                            vseq_vector = build_additive_dictionary(key=(pyr_mut, flanking_nuc, flanking_pos),
                                                                    value=pos_by_type_counts[(str(coord), mut)],
                                                                    dictionary=vseq_vector)
                        elif rcrs_pos2ref[str(coord)] == "A" or rcrs_pos2ref[str(coord)] == "G":
                            pyr_mut = flip_mut(mut)  # 求反向互补
                            flanking_coord = coord - flanking_pos  # 求反向互补
                            flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                            flanking_nuc = flip_base(rcrs_pos2ref[flanking_coord])  # 求反向互补
                            vseq_vector = build_additive_dictionary(key=(pyr_mut, flanking_nuc, flanking_pos),
                                                                    value=pos_by_type_counts[(str(coord), mut)],
                                                                    dictionary=vseq_vector)
    print('\n' + "The count of sequence context around mutations is: ", vseq_vector)
    return vseq_vector


def write_vseq_for_plotting(vseq_vector: Dict[Tuple[str, str, int], int], mut_type_counts: Dict[str, int],
                            flanking_range: List[int], region_name: str, mut_group: str = None):
    """将突变序列上下文向量 Vseq t,p,n 写入文件以便绘图，并先将计数转换为频率。

    :param vseq_vector: 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为计数；对 Tv 而言，突变类型指嘧啶突变（C>A、C>G、T>A、T>G）
    :param mut_type_counts: 字典，键为突变类型，值为在该区域内的总计数
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param region_name: 区域名称，用于生成文件名
    :param mut_group: 可选值 Ts 或 Tv，分别表示转换或颠换，默认同时输出二者
    """
    # 先写出 Ts，再写 Tv
    if mut_group is None or mut_group == "Ts":
        f = write_file_header(
            file='output/sequence_context_vectors/Vseq_mut_%s_%s.txt' % ("Ts", region_name),
            flanking_range=flanking_range, variable="mutation_type")
        for flanking_nuc in nucleotides:
            for mut in class_III_mutations:  # 转换（Ts）
                f.write(flanking_nuc + '\t')
                for flanking_pos in flanking_range:
                    # 转换为频率
                    vseq_freq = vseq_vector[(mut, flanking_nuc, flanking_pos)] / mut_type_counts[mut]
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(vseq_freq) + '\t' + mut + '\n')  # 行末尾写入类型
                    else:
                        f.write(str(vseq_freq) + '\t')

    if mut_group is None or mut_group == "Tv":
        f = write_file_header(
            file='output/sequence_context_vectors/Vseq_mut_%s_%s.txt' % ("Tv", region_name),
            flanking_range=flanking_range, variable="mutation_type")
        for flanking_nuc in nucleotides:
            for mut in mut_types:
                if mut not in class_III_mutations:  # 颠换（Tv）
                    if mut[0] == "C" or mut[0] == "T":  # 仅遍历嘧啶颠换
                        f.write(flanking_nuc + '\t')
                        for flanking_pos in flanking_range:
                            # vseq_vector 仅包含嘧啶颠换（C>A、C>G、T>A、T>G）
                            # 但 mut_type_counts 覆盖 8 种颠换，因此需要累加互补突变（如 C>A 与 G>T）
                            vseq_freq = vseq_vector[(mut, flanking_nuc, flanking_pos)] / (
                                        mut_type_counts[mut] + mut_type_counts[flip_mut(mut)])
                            if flanking_pos == flanking_range[-1]:
                                f.write(str(vseq_freq) + '\t' + mut + '\n')  # 行末尾写入类型
                            else:
                                f.write(str(vseq_freq) + '\t')


def make_ref_freq_vector(mut_group: str, flanking_range: List[int], reference_region: List[int]):
    """生成参考序列上下文 f(ref) n(t),p,n'，统计每个参考碱基在各侧翼位置的侧翼碱基数量。

    转换（Ts）与颠换（Tv）分别计算。对于 Ts，分别统计四种参考碱基 A、C、G、T 的序列上下文；
    对于 Tv，由于需要跨链合并，仅对嘧啶参考碱基 C、T 进行计算。

    :param mut_group: 取值 Ts 或 Tv，分别表示转换或颠换
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param reference_region: 待统计的 mtDNA 区域坐标列表
    :return: mut_group_fref_vector 字典，键为 (参考碱基, 侧翼碱基, 侧翼位置)，值为对应计数；对 Tv 来说，参考碱基为嘧啶（C 或 T）
    """
    # 首先生成坐标到参考碱基的字典
    rcrs_pos2ref = rcrs_pos_to_ref()
    mut_group_fref_vector = {}
    for coord in reference_region:
        if coord != 3107:  # 该位置的参考碱基为 N
            ref_nuc = rcrs_pos2ref[str(coord)]
            for flanking_pos in flanking_range:
                if mut_group == "Ts":
                    # 统计每个参考碱基（A、C、G、T）周围的 A、C、G、T 侧翼碱基数量
                    # 先将 [-10:10]/{0} 范围内的侧翼位置转换为基因组坐标
                    # 再将侧翼坐标映射为参考基因组中的侧翼碱基
                    # 保存每个参考碱基在各侧翼位置的侧翼碱基计数
                    flanking_coord = coord + flanking_pos
                    # 处理人工断点与 m.3107N 间隔附近的坐标
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    flanking_nuc = rcrs_pos2ref[flanking_coord]
                    mut_group_fref_vector = build_dictionary(key=(ref_nuc, flanking_nuc, flanking_pos),
                                                             dictionary=mut_group_fref_vector)
                elif mut_group == "Tv":
                    # 由于计数较少，对 Tv 进行跨链合并
                    # 与上述逻辑类似，但使用嘧啶（C 或 T）作为参考碱基进行记录
                    flanking_coord = pyr_ref_nuc = ''
                    if ref_nuc == "C" or ref_nuc == "T":
                        pyr_ref_nuc = ref_nuc
                        flanking_coord = coord + flanking_pos
                    elif ref_nuc == "A" or ref_nuc == "G":
                        pyr_ref_nuc = flip_base(ref_nuc)
                        flanking_coord = coord - flanking_pos  # 求反向互补
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    flanking_nuc = rcrs_pos2ref[flanking_coord]
                    if ref_nuc == "A" or ref_nuc == "G":
                        flanking_nuc = flip_base(flanking_nuc)  # 求反向互补
                    mut_group_fref_vector = build_dictionary(key=(pyr_ref_nuc, flanking_nuc, flanking_pos),
                                                             dictionary=mut_group_fref_vector)
    print('\n' + "The count of sequence context around reference nucleotides for mut group ", mut_group,
          " is: ", mut_group_fref_vector)
    return mut_group_fref_vector


def write_fref_for_plotting(mut_group: str, mut_group_fref_vector: Dict[Tuple[str, str, int], int],
                            base_proportions: Dict[str, float], flanking_range: List[int],
                            reference_region: List[int], region_name: str):
    """将参考序列上下文 f(ref) n(t),p,n' 写入文件以便绘图，并先将计数转换为频率。

    :param mut_group: 取值 Ts 或 Tv，分别表示转换或颠换
    :param mut_group_fref_vector: 字典，键为 (参考碱基, 侧翼碱基, 侧翼位置)，值为计数；对 Tv 而言，参考碱基为嘧啶（C 或 T）
    :param base_proportions: 字典，键为碱基，值为其在该区域内的占比
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param reference_region: mtDNA 区域的坐标列表，用于统计碱基出现次数
    :param region_name: 区域名称，用于生成文件名
    """
    f = write_file_header(file='output/sequence_context_vectors/f_ref_%s_%s.txt' % (mut_group, region_name),
                          flanking_range=flanking_range, variable="reference_nucleotide")
    for flanking_nuc in nucleotides:
        if mut_group == "Ts":
            for ref_nuc in nucleotides:
                f.write(flanking_nuc + '\t')
                prob_ref_nuc = base_proportions[ref_nuc]
                for flanking_pos in flanking_range:
                    # 转换为频率，再除以参考碱基（A、C、G、T）的出现概率
                    f_ref = mut_group_fref_vector[(ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(f_ref / prob_ref_nuc) + '\t' + ref_nuc + '\n')  # 行末尾写入参考碱基
                    else:
                        f.write(str(f_ref / prob_ref_nuc) + '\t')
        elif mut_group == "Tv":
            for pyr_ref_nuc in ["C", "T"]:
                f.write(flanking_nuc + '\t')
                prob_ref_nuc = ''
                if pyr_ref_nuc == "C":
                    prob_ref_nuc = base_proportions["C"] + base_proportions["G"]
                elif pyr_ref_nuc == "T":
                    prob_ref_nuc = base_proportions["A"] + base_proportions["T"]
                for flanking_pos in flanking_range:
                    # 转换为频率，再除以嘧啶参考碱基（C、T）的出现概率
                    f_ref = mut_group_fref_vector[(pyr_ref_nuc, flanking_nuc, flanking_pos)] / \
                            (len(reference_region) - 1)
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(f_ref / prob_ref_nuc) + '\t' + pyr_ref_nuc + '\n')  # 行末尾写入参考碱基
                    else:
                        f.write(str(f_ref / prob_ref_nuc) + '\t')


def lr_seq_context(mut_type_counts: Dict[str, int], vseq_vector: Dict[Tuple[str, str, int], int],
                   Ts_fref_vector: Union[Dict[Tuple[str, str, int], int], None],
                   Tv_fref_vector: Union[Dict[Tuple[str, str, int], int], None],
                   base_proportions: Dict[str, float], flanking_range: List[int],
                   reference_region: List[int], mut_group: str = None):
    """计算 12 种突变类型的序列上下文似然比。

    :param mut_type_counts: 字典，键为突变类型，值为该区域内的总计数
    :param vseq_vector: 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为计数；对 Tv 而言，突变类型为嘧啶突变（C>A、C>G、T>A、T>G）
    :param Ts_fref_vector: 字典，键为 (参考碱基, 侧翼碱基, 侧翼位置)，值为计数
    :param Tv_fref_vector: 字典，键为 (参考碱基, 侧翼碱基, 侧翼位置)，值为计数；对 Tv 而言，参考碱基为嘧啶（C 或 T）
    :param base_proportions: 字典，键为碱基，值为其在该区域内的占比
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param reference_region: mtDNA 区域的坐标列表，用于统计碱基出现次数
    :param mut_group: 可选值 Ts 或 Tv，分别表示转换或颠换，默认同时处理二者
    :return: lambda_seq_context 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为序列上下文似然比，涵盖全部 12 种突变类型
    """
    lambda_seq_context = {}
    for mut in mut_types:
        if mut in class_III_mutations:  # 转换（Ts）
            if mut_group == "Tv":
                continue  # 如果只分析颠换，则跳过
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # 分子：突变序列上下文 Vseq t,p,n
                    numerator = vseq_vector[(mut, flanking_nuc, flanking_pos)]
                    # 分母：Ts_fref_vector，即参考序列上下文 f(ref) n(t),p,n'
                    # 先除以位点数获得频率，再除以参考碱基的概率作为缩放因子
                    ref_nuc = mut[0]
                    f_ref = Ts_fref_vector[(ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    prob_ref = base_proportions[ref_nuc]
                    f_ref_scaled = f_ref / prob_ref
                    # 最后乘以对应突变类型的计数
                    denominator = mut_type_counts[mut] * f_ref_scaled
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = numerator / denominator
        elif mut in ["C>A", "C>G", "T>A", "T>G"]:  # 嘧啶颠换
            if mut_group == "Ts":
                continue  # 如果只分析转换，则跳过
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # 分子：突变序列上下文 Vseq t,p,n
                    # 注意 vseq_vector 包含全部 4 种转换，但仅包含嘧啶颠换（C>A、C>G、T>A、T>G）
                    numerator = vseq_vector[(mut, flanking_nuc, flanking_pos)]
                    # 分母：Tv_fref_vector，即参考序列上下文 f(ref) n(t),p,n'
                    # 先除以位点数获得频率，再除以嘧啶参考碱基的概率作为缩放因子
                    pyr_ref_nuc = mut[0]
                    f_ref = Tv_fref_vector[(pyr_ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    prob_ref = ''
                    # 先折叠至嘧啶类别
                    if pyr_ref_nuc == "C":
                        prob_ref = base_proportions["C"] + base_proportions["G"]
                    elif pyr_ref_nuc == "T":
                        prob_ref = base_proportions["A"] + base_proportions["T"]
                    f_ref_scaled = f_ref / prob_ref
                    # 然后乘以对应突变类型的计数
                    # 但 mut_type_counts 包含 8 种颠换，需要将互补突变（如 C>A 与 G>T）合并
                    denominator = (mut_type_counts[mut] + mut_type_counts[flip_mut(mut)]) * f_ref_scaled
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = numerator / denominator
    for mut in mut_types:  # 需在嘧啶颠换计算完成后处理
        if mut_group == "Ts":
            continue  # 如果只分析转换，则跳过
        if mut in ["G>T", "G>C", "A>T", "A>C"]:  # 嘌呤颠换
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # 与互补链上的取值相同
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = \
                        lambda_seq_context[(flip_mut(mut), flip_base(flanking_nuc), (flanking_pos * -1))]
    print('\n' + "The likelihood of sequence context around mutation types is: ", lambda_seq_context)
    return lambda_seq_context


def write_lambda_seq_context_for_plotting(lambda_seq_context: Dict[Tuple[str, str, int], float],
                                          lambda_ref_nuc: Dict[str, float], lambda_mut_class: Dict[str, float],
                                          flanking_range: List[int], region_name: str,  mut_group: str = None):
    """将序列上下文的似然比写入文件以便绘图，同时包含参考碱基和突变类别的似然比。

    :param lambda_seq_context: 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为序列上下文似然比
    :param lambda_ref_nuc: 字典，键为参考碱基，值为该碱基的突变似然比
    :param lambda_mut_class: 字典，键为突变类型，值为突变类别的似然比
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param region_name: 区域名称，用于生成文件名
    :param mut_group: 取值 Ts 或 Tv，分别表示转换或颠换，默认同时输出二者
    """
    f = write_file_header(
        file='output/sequence_context_vectors/lambda_seq_context_%s.txt' % region_name,
        flanking_range=flanking_range,
        variable="mutation_type" + '\t' + "lambda_ref_nucleotide" + '\t' + "lambda_mutation_class")
    for mut in mut_types:
        if mut_group == "Ts" and mut not in class_III_mutations:
            continue  # 若仅关注转换，则跳过其他类型
        elif mut_group == "Tv" and mut in class_III_mutations:
            continue  # 若仅关注颠换，则跳过转换
        for flanking_nuc in nucleotides:
            f.write(str(flanking_nuc) + '\t')
            for flanking_pos in flanking_range:
                if flanking_pos == flanking_range[-1]:
                    f.write(str(lambda_seq_context[(mut, flanking_nuc, flanking_pos)]) + '\t' + str(mut) + '\t' +
                            str(lambda_ref_nuc[mut[0]]) + '\t' + str(lambda_mut_class[mut]) + '\n')
                else:
                    f.write(str(lambda_seq_context[(mut, flanking_nuc, flanking_pos)]) + '\t')


def inputs_for_composite_likelihood(denovo_list: str, reference_region: List[int], region_name: str,
                                    flanking_range: List[int], mut_group: str = None):
    """生成复合似然计算所需的各类字典。
    支持仅计算转换、仅计算颠换，或同时计算二者（默认）。用于将 Ori 区域单独处理。

    :param denovo_list: 最终 de novo 变异列表路径，RefPosAlt 格式并附计数
    :param reference_region: 待分析的 mtDNA 区域坐标列表
    :param region_name: 区域名称，用于生成文件名
    :param flanking_range: 侧翼位置列表，范围为 [-10:10]/{0}
    :param mut_group: 取值 Ts 或 Tv，分别表示转换或颠换，默认同时包含二者
    :return: lambda_ref_nuc、lambda_mut_class、lambda_seq_context 三个字典，供后续函数使用
    """
    print('Start with calculating likelihood of mutation at reference nucleotide and likelihood of mutation class')

    denovo_counts = make_denovo_counts(denovo_list=denovo_list)
    mut_type_counts = make_type_count_vector(denovo_counts=denovo_counts, reference_region=reference_region)
    base_proportions = probability_per_ref_nuc(reference_region=reference_region)
    lambda_ref_nuc = lr_ref_nuc(type_counts=mut_type_counts, base_proportions=base_proportions)
    class_probabilities = probability_per_class(base_proportions=base_proportions)
    lambda_mut_class = lr_class(type_counts=mut_type_counts, class_probabilities=class_probabilities)

    print('\n' + 'Now calculate likelihood of mutation sequence context')

    pos_by_type_counts = count_type_per_pos(denovo_counts=denovo_counts)
    vseq_vector = make_mut_context_vector(pos_by_type_counts=pos_by_type_counts, flanking_range=flanking_range,
                                          reference_region=reference_region, mut_group=mut_group)
    write_vseq_for_plotting(vseq_vector=vseq_vector, mut_type_counts=mut_type_counts, flanking_range=flanking_range,
                            region_name=region_name, mut_group=mut_group)

    # 以下函数需要指定突变类别，以便分别处理
    Ts_fref_vector = Tv_fref_vector = {}
    if mut_group is None or mut_group == "Ts":
        Ts_fref_vector = make_ref_freq_vector(mut_group="Ts", flanking_range=flanking_range,
                                              reference_region=reference_region)
        write_fref_for_plotting(mut_group="Ts", mut_group_fref_vector=Ts_fref_vector,
                                base_proportions=base_proportions, flanking_range=flanking_range,
                                reference_region=reference_region, region_name=region_name)
    if mut_group is None or mut_group == "Tv":
        Tv_fref_vector = make_ref_freq_vector(mut_group="Tv", flanking_range=flanking_range,
                                              reference_region=reference_region)
        write_fref_for_plotting(mut_group="Tv", mut_group_fref_vector=Tv_fref_vector,
                                base_proportions=base_proportions, flanking_range=flanking_range,
                                reference_region=reference_region, region_name=region_name)

    lambda_seq_context = lr_seq_context(mut_type_counts=mut_type_counts, vseq_vector=vseq_vector,
                                        Ts_fref_vector=Ts_fref_vector, Tv_fref_vector=Tv_fref_vector,
                                        base_proportions=base_proportions, flanking_range=flanking_range,
                                        reference_region=reference_region, mut_group=mut_group)
    write_lambda_seq_context_for_plotting(lambda_seq_context=lambda_seq_context, lambda_ref_nuc=lambda_ref_nuc,
                                          lambda_mut_class=lambda_mut_class, flanking_range=flanking_range,
                                          region_name=region_name, mut_group=mut_group)

    return lambda_ref_nuc, lambda_mut_class, lambda_seq_context


def composite_likelihood(lambda_ref_nuc: Dict[str, float], lambda_mut_class: Dict[str, float],
                         lambda_seq_context: Dict[Tuple[str, str, int], float],
                         ori_lambda_ref_nuc: Dict[str, float], ori_lambda_mut_class: Dict[str, float],
                         ori_lambda_seq_context: Dict[Tuple[str, str, int], float],
                         context_size: int):
    """计算参考 mtDNA 序列中每个位点、每种突变类型的复合似然。

    :param lambda_ref_nuc: 字典，键为参考碱基，值为该碱基的突变似然比
    :param lambda_mut_class: 字典，键为突变类型，值为突变类别的似然比
    :param lambda_seq_context: 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为序列上下文似然比，涵盖全部 12 种突变类型
    :param ori_lambda_ref_nuc: 字典，键为参考碱基，值为 Ori 区域内的突变似然比
    :param ori_lambda_mut_class: 字典，键为突变类型，值为 Ori 区域内的突变类别似然比
    :param ori_lambda_seq_context: 字典，键为 (突变类型, 侧翼碱基, 侧翼位置)，值为 Ori 区域内的序列上下文似然比
    :param context_size: 序列上下文窗口大小，默认 3（即三核苷酸）
    """
    f = open('output/mutation_likelihoods/mito_mutation_likelihoods.txt', "w")
    f.write("POS	REF	ALT	Likelihood" + '\n')

    # 窗口在参考位点两侧各包含多少个碱基；由于 range 右端不包含，需要额外加 1
    pos_range = list(range(1, ((context_size - 1) // 2) + 1))

    # 先生成坐标到参考碱基的映射字典
    rcrs_pos2ref = rcrs_pos_to_ref()

    # 遍历线粒体基因组中所有可能的突变
    for row in csv.DictReader(open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        coord = int(row["POS"])
        ref = row["REF"]
        alt = row["ALT"]
        mut = ref + ">" + alt

        if coord != 3107:  # 该位置的参考碱基为 N
            # product_mut_type 表示参考碱基的突变似然与突变类别似然的乘积；Ori 区域单独处理
            if coord in ori_region:
                product_mut_type = ori_lambda_ref_nuc[ref] * ori_lambda_mut_class[mut]
            else:
                product_mut_type = lambda_ref_nuc[ref] * lambda_mut_class[mut]

            # 对每个参考位点的突变类型初始化序列上下文乘积
            product_seq_context = 0

            for number in pos_range:  # 取决于上下文窗口大小，这里的数值范围为 1-10
                # 将侧翼位置转换为具体坐标
                flanking_coord_list = [coord - number, coord + number]
                for flanking_coord in flanking_coord_list:
                    # 标记该侧翼位置位于参考位点的左侧还是右侧
                    flanking_pos = ''
                    if flanking_coord < coord:
                        flanking_pos = -1 * number
                    elif flanking_coord > coord:
                        flanking_pos = number
                    # 该函数可处理人工断点和 m.3107N 间隔附近的坐标
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    # 该侧翼坐标对应的参考碱基
                    flanking_nuc = rcrs_pos2ref[flanking_coord]

                    # 计算各侧翼位置的序列上下文似然乘积；Ori 区域内的转换单独处理
                    if (coord in ori_region) and (mut in class_III_mutations):  # OriB-OriH 区域的转换
                        if product_seq_context == 0:
                            product_seq_context = ori_lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                        else:
                            product_seq_context = product_seq_context * \
                                                  ori_lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                    else:
                        if product_seq_context == 0:
                            product_seq_context = lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                        else:
                            product_seq_context = product_seq_context * \
                                                  lambda_seq_context[(mut, flanking_nuc, flanking_pos)]

            # 最终复合似然：突变类型似然与序列上下文似然的乘积
            product_pos = product_mut_type * product_seq_context
            # 写入文件
            f.write(str(coord) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(product_pos) + '\n')


if __name__ == "__main__":
    print(datetime.datetime.now(), '\n' + "Starting to build composite likelihood model for mtDNA!" + '\n')

    parser = argparse.ArgumentParser()
    parser.add_argument("-context_size", type=int,
                        help="how many nucleotides to include when calculating sequence context, ie 3 is trinucleotide")
    parser.add_argument("-denovo_list", type=str,
                        help="path to the de novo mutation list, in RefPosAlt format with their counts")
    args = parser.parse_args()

    # 设置默认参数
    if args.context_size is None:
        args.context_size = 3
    if args.denovo_list is None:
        args.denovo_list = 'output/denovo/final_denovo.txt'

    for path in ['output/sequence_context_vectors/', 'output/mutation_likelihoods/']:
        if not os.path.exists(path):
            os.makedirs(path)
    print("Creating required directories")

    print('\n' + 'This script will calculate the likelihood of mutation in the mitochondrial genome')
    print('\n' + 'This will be done using window size of ' + str(args.context_size) + ' for sequence context')
    # 使用双重整除以确保结果为整数
    flanking_range = [i for i in
                      list(range(-((args.context_size - 1) // 2), (((args.context_size - 1) // 2) + 1)))
                      if i != 0]
    print('\n' + 'Mutation likelihoods in the OriB-OriH region from ' + str(start_ori) + '-' + str(end_ori) +
          ' will be calculated separately')

    print('\n' + 'First, start with the reference region excluding the OriB-OriH' + '\n')

    (lambda_ref_nuc, lambda_mut_class, lambda_seq_context) = \
        inputs_for_composite_likelihood(denovo_list=args.denovo_list, reference_region=reference_except_ori,
                                        region_name="", flanking_range=flanking_range)

    print('\n' + 'Next, calculate for the OriB-OriH region' + '\n')

    # 仅对转换（Ts）执行该计算
    (ori_lambda_ref_nuc, ori_lambda_mut_class, ori_lambda_seq_context) = \
        inputs_for_composite_likelihood(denovo_list=args.denovo_list, reference_region=ori_region,
                                        region_name="OriB-OriH", flanking_range=flanking_range, mut_group="Ts")

    print('\n' + 'Now, calculate position mutability across mitochondrial genome and write to file' + '\n')

    composite_likelihood(lambda_ref_nuc=lambda_ref_nuc, lambda_mut_class=lambda_mut_class,
                         lambda_seq_context=lambda_seq_context, ori_lambda_ref_nuc=ori_lambda_ref_nuc,
                         ori_lambda_mut_class=ori_lambda_mut_class, ori_lambda_seq_context=ori_lambda_seq_context,
                         context_size=args.context_size)

    print(datetime.datetime.now(), '\n' + "Finished building composite likelihood model for mtDNA!")
