## 概述

### `databases`（数据库）

包含的文件与数据来源：

* gnomAD 数据集，来自 [gnomAD 浏览器](https://gnomad.broadinstitute.org/downloads)。
* HelixMTdb 数据集，来自 [Helix 官网](https://www.helix.com/pages/mitochondrial-variant-database)。
* MITOMAP 数据集，来自 [MITOMAP](https://www.mitomap.org/MITOMAP/resources)。
* 疾病相关变异，来自 [MITOMAP](https://www.mitomap.org/MITOMAP/resources)。
* 疾病相关变异，来自 [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)。
* 单倍群（Haplogroup）相关变异，提取自 [PhyloTree](https://www.phylotree.org/)。
* 线粒体基因组位点列表，来自 [MITOMAP](https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci)。

---

### `input_denovo`（denovo突变输入）

用于构建突变模型的线粒体denovo突变。
所有公开发表的数据集均包含在此；由于数据限制，未发表的数据集未包含。
关于denovo突变数据集及其整理方式的细节，请参见论文补充材料。

---

### `insilicos`（计算注释）

* 100 种脊椎动物比对生成的 PhyloP 保守性评分，来自 [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way/)。
* APOGEE 功能预测评分，来自 [MitImpact](https://mitimpact.css-mendel.it/)。
* MitoTip 预测评分，来自 [MITOMAP](https://www.mitomap.org/MITOMAP/MitoTipScores)。
* HmtVar 疾病风险评分，来自 [HmtVar](https://www.hmtvar.uniba.it/)。

---

### `other_annotations`（其他注释）

* UniProt 蛋白知识库注释，来自 [UniProt FTP](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/)。
  对于线粒体编码蛋白，使用 *UP000005640_9606_beds* 目录中所有人类 bed 文件提取。
* 复合物 I 质子转移相关残基，整理自 [PMID:32972993](https://pubmed.ncbi.nlm.nih.gov/32972993/)。
* tRNA 位置编号，源文件来自 [MitoTIP](https://github.com/sonneysa/MitoTIP/)，并手动转换为各 mtDNA 坐标的 tRNA 位置列表。
* RNA 修饰与 tRNA 结构域，来自 [Lake 等人研究](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356)。
* RNA 碱基类型，依据 [Lake 等人](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356) 手动整理的数据文件 `all_RNA_bases.tsv`。
* 参与 rRNA:rRNA 桥接的碱基，整理自 [PMID:25838379](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4501431/)。
* 人类与黑猩猩线粒体参考序列比对，使用 [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) 的参考序列生成。
* 英国生物银行中线粒体基因组变异数据，来自 [Hong, Battle 等人](https://www.nature.com/articles/s41467-023-41785-7)。

#### *以下为核基因相关数据：*

* 核基因错义突变约束指标，来自 [Karczewski 等人](https://www.nature.com/articles/s41586-020-2308-7)。
* 必需基因数据，来自 [IMPC](https://www.ebi.ac.uk/mi/impc/essential-genes-search/)。
* 发育障碍基因数据，来自 [DECIPHER](https://www.deciphergenomics.org/ddd/ddgenes)。

---

### `synthetic_vcf`（合成VCF）

该目录下文件的生成方式参考 [原始仓库说明](https://github.com/broadinstitute/gnomad-mitochondria/tree/main/gnomad_mitochondria/manuscript_analyses)。

---

### `svg_input`（SVG输入）

用于生成 tRNA 和 rRNA 二级结构图的输入文件，来源见 [Lake 等人研究](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356)。

---

### `pdb`（蛋白结构文件）

用于生成三维图像、动画以及计算区域约束距离的结构文件，来自 [Protein Data Bank](https://www.rcsb.org/) 或 [UniProt](https://www.uniprot.org/)。
