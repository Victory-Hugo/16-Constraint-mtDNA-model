#!/bin/bash
#!不建议一次性运行，因为代码中有绝对路径
#!千万不要修改python代码文件名
#!如下路径千万不要修改，也不要重命名文件名、文件路径等:
# ├── 0-required_files
# │   ├── README.md
# │   ├── databases
# │   ├── input_denovo
# │   ├── insilicos
# │   ├── other_annotations
# │   ├── pdb
# │   ├── svg_input
# │   └── synthetic_vcf
# ├── 1-构建模型
# │   ├── 1-主执行脚本.sh
# │   ├── annotate_mutations.py
# │   ├── compare_denovo.py
# │   ├── compile_denovo.py
# │   ├── composite_likelihood_mito.py
# │   ├── filter_denovo.py
# │   └── 蒙特卡洛模拟异质性.R
# ├── 2-计算观测值和预期值
# │   ├── 1-主执行脚本.sh
# │   ├── README.md
# │   ├── calibrate_model_step1.py
# │   ├── calibrate_model_step2.R
# │   ├── oe_consequences.py
# │   ├── oe_functions.py
# │   ├── oe_loci.py
# │   ├── oe_lookup.py
# │   └── oe_other.py


#*========================第一部分========================
#*这部分使用的是gnomAD数据集，使用默认参数，没有排除任何位点。
#*========================第一部分========================
# 在 gnomAD 中跨整个线粒体基因组及各基因计算不同功能类别变异的观测/预期比
# 这些脚本已设置为使用 gnomAD 分析所需的默认参数，无需额外提供
cd /mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-整理后/

python3 2-计算观测值和预期值/calibrate_model_step1.py
Rscript 2-计算观测值和预期值/calibrate_model_step2.R
python3 2-计算观测值和预期值/oe_consequences.py
python3 2-计算观测值和预期值/oe_loci.py
python3 2-计算观测值和预期值/oe_other.py

#* 如果需要计算,请给出起始和终止位置：
python3 \
    2-计算观测值和预期值/oe_lookup.py \
    -start 1000 \
    -end 2000

#*========================第二部分========================
#*这部分使用的是gnomAD数据集，使用默认参数，没有排除任何位点。
#*========================第二部分========================
# 在含异质性的复现数据集 HelixMTdb 中计算不同功能类别变异的观测/预期比
# 这部分使用的是HelixMTdb数据集作为复现验证数据集，有以下关键差异：

# 数据源不同：-obs 'helix_max_hl' 参数指定使用HelixMTdb数据
# 输出目录不同：-prefix 'replication_dataset/helix_' 输出到复现数据集目录
# 排除特定位点：-exc_sites 参数排除了300-316、513-525、16182-16194这些位点的变异

# 这是一个典型的主分析+复现验证的研究设计：
    # 主分析：使用gnomAD数据集进行线粒体DNA约束模型的观测/预期比计算
    # 复现验证：使用独立的HelixMTdb数据集验证模型的稳健性和可重复性
    # 通过在两个不同数据集上运行相同的分析流程，可以确保研究结果的可靠性和泛化能力。



mkdir output/oe/replication_dataset/
mkdir output/calibration/replication_dataset/


python3 2-计算观测值和预期值/calibrate_model_step1.py \
    -obs 'helix_max_hl' \
    -prefix 'replication_dataset/helix_' \
    -exc_sites 300 301 302 303 304 305 306 307 308 309 310 \
    311 312 313 314 315 316 513 514 515 516 517 518 519 520 \
    521 522 523 524 525 16182 16183 16184 16185 16186 16187 \
    16188 16189 16190 16191 16192 16193 16194

Rscript 2-计算观测值和预期值/calibrate_model_step2.R \
    'replication_dataset/helix_' \
    'no'

python3 2-计算观测值和预期值/oe_consequences.py \
    -obs 'helix_max_hl' \
    -parameters 'output/calibration/replication_dataset/helix_linear_model_fits.txt' \
    -prefix 'replication_dataset/helix_' \
    -exc_sites 300 301 302 303 304 305 306 307 308 309 310 \
    311 312 313 314 315 316 513 514 515 516 517 518 519 520 \
    521 522 523 524 525 16182 16183 16184 16185 16186 16187 \
    16188 16189 16190 16191 16192 16193 16194
