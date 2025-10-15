import csv
import datetime
from typing import Dict, TextIO, Union
import os


# 辅助函数
def build_sample_dict(sample: str, dict: Dict[str, int]):
    """生成包含每个样本 de novo 数量的字典。

    :param sample: 样本名称
    :param dict: 字典名称
    :return: 一个字典，键为样本名，值为该样本的 de novo 计数
    """
    if sample not in dict:
        dict[sample] = 1
    else:
        dict[sample] += 1
    return dict


def write_denovo(file: TextIO, sample: str, variant: str, dict: Union[Dict[str, int], int]):
    """将 de novo 突变及其所在样本和样本的 de novo 数量写入文件。

    :param file: 输出文件对象
    :param sample: 样本名称
    :param variant: de novo 突变（RefPosAlt 格式）
    :param dict: 一个字典，键为样本名，值为样本的 de novo 数量，也可以直接提供整数计数
    """
    if type(dict) is int:
        sample_count = dict
    else:
        sample_count = dict[sample]
    file.write(variant + '\t' + sample + '\t' + str(sample_count) + '\n')


def rcrs_pos_to_ref():
    """生成一个字典，将每个位点与 rCRS 中的参考碱基对应起来。

    :return: 字典，键为 rCRS 位点，值为该位点的参考碱基
    """
    dictionary = {}
    for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        dictionary[row["POS"]] = row["REF"]
    return dictionary


# 提取 de novo 的函数
def extract_denovo(file: TextIO):
    """解析来自文献和自有数据集的 de novo 列表。
    每个文件的格式各不相同，因此需单独处理。
    共三类：生殖系（germline）、体细胞组织（somatic tissue）和体细胞癌症（somatic cancer）de novo 突变。

    :param file: 要写入 de novo 突变、样本及样本 de novo 数量的 txt 文件
    """
    # 首先构建 sample_counts 字典，统计每个样本的 de novo 数量
    # 然后将每个 de novo 及对应样本信息写入文件
    sample_counts = {}

    # 生殖系（GERMLINE）

    # 数据集1 - 来自 Wei 等人 2019 Science

    for row in csv.reader(open('required_files/input_denovo/germline/PMID31123110_Data_S1.txt'), delimiter='\t'):
        if row[6] == "de novo":
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID31123110" + "-germline"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID31123110_Data_S1.txt'), delimiter='\t'):
        if row[6] == "de novo":
            write_denovo(file=file, sample=(row[0] + "-PMID31123110" + "-germline"), variant=row[2], dict=sample_counts)

    # 数据集2 - 来自 Rebolledo-Jaramillo 等人 2014 PNAS

    for row in csv.reader(open('required_files/input_denovo/germline/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "child":  # “child” 类别表示生殖系 de novo 突变，见原文
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID25313049" + "-germline"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "child":
            write_denovo(file=file, sample=(row[0] + "-PMID25313049" + "-germline"), variant=(row[2] + row[1] + row[3]),
                         dict=sample_counts)

    # 数据集3 - 来自 Zaidi 等人 2019 PNAS

    catch_list = []  # 用于捕获同一个个体被列出两次的 de novo（每个组织类型一行）
    for row in csv.reader(open('required_files/input_denovo/germline/PMID31757848_TableS3.txt'), delimiter='\t'):
        if row[6] == "1":  # 启发式标记法识别，详见原文
            sample = row[2] + "-PMID31757848" + "-germline"
            variant = row[9] + row[4] + row[10]
            if (sample, variant) not in catch_list:
                catch_list.append((sample, variant))
                sample_counts = build_sample_dict(sample=sample, dict=sample_counts)

    catch_list = []
    for row in csv.reader(open('required_files/input_denovo/germline/PMID31757848_TableS3.txt'), delimiter='\t'):
        if row[6] == "1":
            sample = row[2] + "-PMID31757848" + "-germline"
            variant = row[9] + row[4] + row[10]
            if (sample, variant) not in catch_list:
                catch_list.append((sample, variant))
                write_denovo(file=file, sample=sample, variant=variant, dict=sample_counts)

    # 数据集4 - 来自 Li 等人 2016 Gen Res

    for row in csv.reader(open('required_files/input_denovo/germline/PMID26916109_Dataset_S1_Heteroplasmy_list.txt'),
                          delimiter='\t'):
        if row[8] == "Yes":  # de novo
            sample_counts = build_sample_dict(sample=(row[0] + row[1] + "-PMID26916109" + "-germline"),
                                              dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID26916109_Dataset_S1_Heteroplasmy_list.txt'),
                          delimiter='\t'):
        if row[8] == "Yes":  # de novo
            write_denovo(file=file, sample=(row[0] + row[1] + "-PMID26916109" + "-germline"),
                         variant=(row[5] + row[3] + row[6]), dict=sample_counts)

    # 数据集5 - 自有 SPARK 分析结果
    # 由于数据限制无法公开
    if os.path.exists('required_files/input_denovo/germline/SPARK_mtDNA_de_novo.txt'):
        for row in csv.DictReader(open('required_files/input_denovo/germline/SPARK_mtDNA_de_novo.txt'), delimiter='\t'):
            s_counts = int(row["number"])  # 直接读取样本的 de novo 数量
            write_denovo(
                file=file, sample=("sample_from" + "-SPARK" + "-germline"), variant=row["de_novo"], dict=s_counts)
        
        
    # 体细胞组织（SOMATIC TISSUE）

    # 来自 germline 数据集2（Rebolledo-Jaramillo 等人 2014 PNAS）

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "somatic-gain":  # “somatic-gain” 表示体细胞 de novo 突变
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID25313049" + "-somatic_tissue"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "somatic-gain":
            write_denovo(file=file, sample=(row[0] + "-PMID25313049" + "-somatic_tissue"),
                         variant=(row[2] + row[1] + row[3]), dict=sample_counts)

    # 来自 germline 数据集3（Zaidi 等人 2019 PNAS）

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID31757848_TableS4.txt'), delimiter='\t'):
        sample_counts = build_sample_dict(sample=(row[0] + row[4] + "-PMID31757848" + "-somatic_tissue"),
                                          dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID31757848_TableS4.txt'), delimiter='\t'):
        if not (row[0].startswith('Table S4')) and not (row[0].startswith('SampleID')):  # 跳过表头行
            write_denovo(file=file, sample=(row[0] + row[4] + "-PMID31757848" + "-somatic_tissue"),
                         variant=(row[7] + row[5] + row[8]), dict=sample_counts)

    # 数据集6 - gtex，来自 Ludwig 等人 2019 Cell

    for row in csv.reader(
            open('required_files/input_denovo/somatic_tissue/PMID30827679_NIHMS1518665-supplement-10.txt'),
            delimiter='\t'):
        if not (row[0].startswith('Table')) and not (row[0].startswith('Mutation')):
            sample_counts = build_sample_dict(sample=(row[1] + "-PMID30827679" + "-somatic_tissue"), dict=sample_counts)

    # 从 GTEx 数据集中将 Yoruba 参考序列转换为 rCRS 需要 rCRS 查找表
    rcrs_pos2ref = rcrs_pos_to_ref()
    # Yoruba 序列（L3e2b1a1）与 rCRS（H2a2a1）的差异位点
    # 参考 https://haplogrep.i-med.ac.at/2014/09/08/rcrs-vs-rsrs-vs-hg19/
    yoruba_rCRS_diff = [73, 150, 195, 263, 309, 315, 408, 750, 1438, 2352, 2483, 2706, 3107, 4769, 5580, 7028, 8701,
                        8860, 9377, 9540, 10398, 10819, 10873, 11017, 11719, 11722, 12705, 12850, 14212, 14580, 14766,
                        14905, 15301, 15326, 15932, 16172, 16183, 16189, 16193, 16223, 16320, 16519]

    for row in csv.reader(
            open('required_files/input_denovo/somatic_tissue/PMID30827679_NIHMS1518665-supplement-10.txt'),
            delimiter='\t'):
        if not (row[0].startswith('Table')) and not (row[0].startswith('Mutation')):
            pos = int(row[0].split('_')[0])
            alt = row[0].split('_')[1]

            # 该数据似乎基于 Yoruba 参考序列，因此需转换为 rCRS
            # 参考 https://www.mitomap.org/foswiki/bin/view/MITOMAP/YorubanConversion
            if 311 <= pos <= 316:  # 在 309 之前相同
                pos = pos - 1
            elif 318 <= pos <= 3108:
                pos = pos - 2
            elif 3109 <= pos <= 16190:
                pos = pos - 1
            elif 16192 <= pos <= 16571:
                pos = pos - 2
            elif pos == 310 or pos == 317 or pos == 16191:
                continue
            if pos in yoruba_rCRS_diff:
                continue  # 跳过这些差异位点

            ref = rcrs_pos2ref[str(pos)]
            write_denovo(file=file, sample=(row[1] + "-PMID30827679" + "-somatic_tissue"),
                         variant=(ref + str(pos) + alt), dict=sample_counts)

    # 体细胞癌症（SOMATIC CANCER）

    # 数据集7 - 来自 Yuan 等人 2020 Nature Genetics

    for row in csv.DictReader(open('required_files/input_denovo/somatic_cancer/TCMA-MutationSNV.tsv'), delimiter='\t'):
        sample_counts = build_sample_dict(sample=(row["sample_id"] + '-PMID32024997' + '-somatic_cancer'),
                                          dict=sample_counts)

    for row in csv.DictReader(open('required_files/input_denovo/somatic_cancer/TCMA-MutationSNV.tsv'), delimiter='\t'):
        sample = row["sample_id"] + '-PMID32024997' + '-somatic_cancer'
        # 删除作者指出的7个高突变样本（定义为 >13 个体细胞突变）
        if sample_counts[sample] < 14:
            write_denovo(file=file, sample=sample, variant=(row["ref"] + row["position"] + row["var"]),
                         dict=sample_counts)


if __name__ == "__main__":
    print(datetime.datetime.now(), "开始整合所有 de novo!")

    if not os.path.exists('output_files/denovo'):
        os.makedirs('output_files/denovo')
    print(datetime.datetime.now(), "正在创建所需目录")

    # 创建输出文件，用于记录所有 de novo 突变及其样本信息
    f = open('output_files/denovo/all_denovo.txt', "w")
    f.write("denovo	sample	sample_denovo_count" + '\n')

    extract_denovo(file=f)

    print(datetime.datetime.now(), "脚本运行完成!")
