#!/bin/bash
#!/bin/bash
#!不建议一次性运行，因为代码中有绝对路径
#!千万不要修改python代码文件名
#!如下路径千万不要修改，也不要重命名文件名、文件路径等:
# ├── 0-required_files
# │   ├── databases
# │   ├── input_denovo
# │   ├── insilicos
# │   ├── other_annotations
# │   ├── pdb
# │   ├── svg_input
# │   ├── synthetic_vcf
# │   └── 阅读我.md
# ├── 1-构建模型
# │   ├── 1-主执行脚本.sh
# │   ├── __pycache__
# │   ├── annotate_mutations.py
# │   ├── compare_denovo.py
# │   ├── compile_denovo.py
# │   ├── composite_likelihood_mito.py
# │   ├── filter_denovo.py
# │   └── 蒙特卡洛模拟异质性.R
# ├── 1-警告.md
# ├── 2-计算观测值和预期值
# │   ├── 1-主执行脚本.sh
# │   ├── __pycache__
# │   ├── calibrate_model_step1.py
# │   ├── calibrate_model_step2.R
# │   ├── oe_consequences.py
# │   ├── oe_functions.py
# │   ├── oe_loci.py
# │   ├── oe_lookup.py
# │   ├── oe_other.py
# │   └── 阅读我.md
# ├── 3-局部约束
# │   ├── 1-主执行脚本.sh
# │   ├── __pycache__
# │   ├── local_constraint.py
# │   └── 阅读我.md
# ├── 4-区域约束
# │   ├── 1-主执行脚本.sh
# │   ├── __pycache__
# │   ├── annotate_regional_constraint.py
# │   ├── apply_FDR_filter.py
# │   ├── regional_constraint.py
# │   └── 阅读我.md

cd /mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-群体版本

python3 3-局部约束/local_constraint.py