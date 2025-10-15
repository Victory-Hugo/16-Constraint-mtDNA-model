import argparse
"""
该脚本用于过滤和统计不同类型的de novo变异（新生突变），并将结果写入输出文件。

功能说明:
- 通过命令行参数指定三类de novo变异的最大计数阈值:
    - germline_max: 生殖系de novo的最大计数阈值
    - som_tissue_max: 体细胞组织de novo的最大计数阈值
    - som_cancer_max: 体细胞癌症de novo的最大计数阈值（默认为1）
- 调用`filter_denovo`函数进行实际过滤和计数
- 输出最终用于计算突变可能性分数的de novo列表及其计数到指定文件

参数:
    -germline_max (int): 生殖系de novo最大计数阈值
    -som_tissue_max (int): 体细胞组织de novo最大计数阈值
    -som_cancer_max (int): 体细胞癌症de novo最大计数阈值

输出:
    output/denovo/final_denovo.txt: 包含每个de novo及其计数的文件

使用方法:
    python filter_denovo.py -germline_max <值> -som_tissue_max <值> -som_cancer_max <值>
"""
from compare_denovo import filter_denovo
import datetime


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
