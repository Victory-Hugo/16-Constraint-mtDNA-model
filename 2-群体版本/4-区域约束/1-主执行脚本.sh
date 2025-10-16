#!/bin/bash
#!/bin/bash
#!不建议一次性运行，因为代码中有绝对路径
#!千万不要修改python代码文件名
#!如下路径千万不要修改，也不要重命名文件名、文件路径等:



cd /mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-群体版本
#* ================第一步================
python3 4-区域约束/regional_constraint.py
#* ================第二步================
python3 4-区域约束/annotate_regional_constraint.py
#* ================第三步================
python3 4-区域约束/apply_FDR_filter.py
