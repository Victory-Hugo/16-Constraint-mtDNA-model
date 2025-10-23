import argparse
from compare_denovo import filter_denovo
import datetime

#* ====================输入文件====================
#* output/denovo/all_denovo.txt(TSV)
#* ===============================================

#* ====================输出文件====================
#* output/denovo/final_denovo.txt(TSV，含denovo / count列)。
#* 它确保用来建模的 de novo 集合来自“正常”样本，减少异常突变负荷对 mtDNA 变异率估计的干扰。
#* ===============================================

#* ====================输出文件格式================
# denovo	count
# T12586C	1
# T13819C	1
# A1409G	2
# A214G	8
# T6351C	2
# G1747A	3
#* ===============================================

#* ====================脚本功能====================
#* 真正的筛选发生在 compare_denovo.py 里的 filter_denovo 函数中。
# 举个简化例子：假设 all_denovo.txt 有两行
# denovo    sample             sample_denovo_count
# A1234G    sample1-germline   2
# A1234G    sample2-germline   10
# 默认阈值对生殖系是 None（全部保留），返回 {'A1234G': 2+10}。如果调用时把 germline_max 设为 5，那么第二行会被丢弃（因为 10>5），只剩第一行，返回 {'A1234G': 1}。
# filter_denovo.py 得到这个字典后，就写成：
# denovo    count
# A1234G    1
# 这就是“筛选 + 汇总”的过程。
#* ===============================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-germline_max", type=int,
                        help="Maximum sample de novo count for inclusion, for germline de novo")
    parser.add_argument("-som_tissue_max", type=int,
                        help="Maximum sample de novo count for inclusion, for somatic tissue de novo")
    parser.add_argument("-som_cancer_max", type=int,
                        help="Maximum sample de novo count for inclusion, for somatic cancer de novo")
    args = parser.parse_args()

    # 设置默认值，注意 germline_max 和 som_tissue_max 默认为 None，此处不希望对它们应用阈值
    if args.som_cancer_max is None:
        args.som_cancer_max = 1

    print(datetime.datetime.now(), "Starting to filter de novo!")
    print('\n' + "This will produce a final list of de novo that are used to calculate mutational likelihood scores")

    f = open('output/denovo/final_denovo.txt', "w")
    f.write("denovo	count" + '\n')

    denovo_counts = {}
    denovo_counts = filter_denovo(germline_max=args.germline_max, som_tissue_max=args.som_tissue_max,
                                  som_cancer_max=args.som_cancer_max, denovo_counts=denovo_counts)

    for key in denovo_counts:
        f.write(key + '\t' + str(denovo_counts[key]) + '\n')

    print(datetime.datetime.now(), "Script finished!")
