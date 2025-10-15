## 概述

**`calibrate_model_step1.py`**：
针对每个基因或位点，计算中性变异的观测最大杂合度总和及其突变似然得分总和，用于拟合线性模型。用户可通过参数 *-input* 指定输入文件（含似然得分与观测值）；默认输入为 `annotate_mutations.py` 的输出。

**`calibrate_model_step2.R`**：
使用 `calibrate_model_step1.py` 的输出结果，对中性变异数据进行线性模型拟合。

**`oe_functions.py`**：
提供计算观测值与预期值（observed:expected）比值的函数。

**`oe_consequences.py`**：
计算不同功能类别变异的观测/预期比值及其置信区间。用户可通过参数 *-input* 指定输入文件，默认使用 `annotate_mutations.py` 的输出；其他参数（*-obs*, *-parameters*, *-prefix*, *-exc_sites*）默认使用 gnomAD 设置。

**`oe_loci.py`**：
计算不同基因或位点中，不同功能类别变异的观测/预期比值及置信区间。用户可通过参数 *-input* 指定输入文件，默认使用 `annotate_mutations.py` 的输出；其他参数（*-obs*, *-parameters*, *-prefix*, *-exc_sites*）默认使用 gnomAD 设置。

**`oe_other.py`**：
计算其他类别变异的观测/预期比值及置信区间。用户可通过参数 *-input* 指定输入文件，默认使用 `annotate_mutations.py` 的输出；其他参数（*-obs*, *-parameters*, *-prefix*, *-exc_sites*）默认使用 gnomAD 设置。

**`calculate_oe_gnomad.sh`**：
执行上述分析流程，用于评估 gnomAD 数据中的变异约束性。

**`replicate_oe_helix.sh`**：
在验证数据集 HelixMTdb 上执行相同分析，用于评估复制性。用户可为 HelixMTdb 提供特定参数，以覆盖 gnomAD 的默认设置。

**`oe_lookup.py`**：
使用 gnomAD 数据，在指定的坐标区间内计算观测/预期比值及其置信区间。
