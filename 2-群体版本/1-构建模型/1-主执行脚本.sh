#!/bin/bash
#!不建议一次性运行，因为代码中有绝对路径
#!千万不要修改python代码文件名
#!如下路径千万不要修改，也不要重命名文件名、文件路径等:


# 构建突变模型以计算期望值
cd "/mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-群体版本"

python3 1-构建模型/compile_denovo.py \
    --vcf example.vcf.gz \
    --max-af 0.01 \
    --min-carriers 1

python3 1-构建模型/composite_likelihood_mito.py

python3 1-构建模型/annotate_mutations.py \
    --variant-stats \
    output/denovo/variant_carrier_stats.txt
