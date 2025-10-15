## 概述

本系列脚本用于构建新的线粒体突变模型。

`compile_denovo.py`：解析来自文献及内部数据集的 **de novo** 突变列表，输出包含所有 **de novo** 突变及其来源的汇总文件。文献中使用的数据集存储于 `required_files` 目录中。**需注意，因数据使用限制，未公开的内部数据集未包含于此代码库中，因此仅使用公开数据集时的结果会略有差异。** 所用 **de novo** 数据集的详细信息见论文补充材料。

`compare_denovo.py`：比较不同来源（生殖系、体细胞组织、体细胞癌症）及不同样本 **de novo** 数量下的转换突变似然。

`filter_denovo.py`：移除来源于异常样本的 **de novo** 突变，输出用于计算突变率的 **de novo** 突变列表。用户可通过参数 *-germline_max*、*-som_tissue_max*、*-som_cancer_max* 指定各来源中样本 **de novo** 数量的最大阈值；默认值依据分析结果设定。

`composite_likelihood_mito.py`：实现线粒体复合似然模型计算，为每个位点生成单核苷酸变异的突变似然评分。用户可通过参数 *-context_size* 指定序列上下文范围（默认 3，即三核苷酸上下文），并通过 *-denovo_list* 指定所使用的 **de novo** 列表路径（默认为 `filter_denovo.py` 的输出）。

`annotate_mutations.py`：对 `composite_likelihood_mito.py` 的输出进行功能注释，包括变异后果、基因定位及计算预测结果。用户可通过参数 *-input* 指定输入文件（默认为 `composite_likelihood_mito.py` 的输出）。

`build_model.sh`：整合上述分析流程，生成带有功能注释和突变似然评分的线粒体变异列表。运行命令如下：

```
bash -e build_model/build_model.sh
```

`simulate_heteroplasmy.R`：构建生殖系线粒体 DNA 突变与异质性漂变的计算模型，用于验证突变率与最大异质性之间的相关性。

#### 模拟运行示例

运行生殖系线粒体 DNA 突变及异质性漂变的计算模型。在 12 核处理器上运行时间约为 5 分钟。

```
Rscript simulate_heteroplasmy.R
```

生成的文本输出文件用于脚本 `figure_scripts/FigureS4.R` 的后续绘图分析。
