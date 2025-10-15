library(dplyr)
library(stringr)

# comparing nuclear missense constraint to mtDNA missense constraint

# inititate output file
sink(sprintf("../output_files/other/assess_missense_constraint.txt"), append = FALSE)
print("Mitochondrial vs nuclear genome missense constraint")

# mean mitochondrial missense constraint
mito <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep = "\t")
print("Mitochondrial missense constraint - mean observed:expected ratio")
print(round(mean(mito[mito$consequence == "missense", c("obs.exp")]), 2))
print("Mitochondrial missense constraint - mean OEUF")
print(round(mean(mito[mito$consequence == "missense", c("upper_CI")]), 2))

# mean nuclear missense constraint
file <- read.delim(file = '../required_files/other_annotations/PMID32461654_supplementary_dataset_11_full_constraint_metrics.tsv', header = TRUE, sep = "\t")

# first bin by missense OEUF decile & filter to canonical, remove NA and those with flags
file <- file[file$canonical == "true" & !is.na(file$oe_mis_upper) & file$constraint_flag == "",] %>% 
  mutate(oe_mis_upper_decile = ceiling((rank(oe_mis_upper)/length(oe_mis_upper)) * 10))

print("Nuclear missense constraint - mean observed:expected ratio")
print(round(mean(file[, c("oe_mis")]), 2))
print("Nuclear missense constraint - mean OEUF")
print(round(mean(file[, c("oe_mis_upper")]), 2))
print("Nuclear missense constraint - lowest decile OEUF summary")
print(summary(file[file$oe_mis_upper_decile == "1", c("oe_mis_upper")]))

# which nuclear genes are essential in human cell lines, or required for development or viability in mice?
impc <- read.delim(file = '../required_files/other_annotations/impc_essential_genes_full_dataset.csv', header = TRUE, sep = ",")

# collapse since multiple rows per gene & filter those with low orthology
impc_unique <- as.data.frame(impc[impc$orthologue_category != "LOW", c("human_ensembl_gene_acc_id", "fusil_bin_code") ] %>% 
                                group_by(human_ensembl_gene_acc_id) %>% 
                                summarise(across(everything(), str_c, collapse=","))) 
file$human_ensembl_gene_acc_id <- file$gene_id
merged <- merge(file, impc_unique, by = "human_ensembl_gene_acc_id", all.x = TRUE)

# aggregate into a group and remove outliers
merged$essential <- ifelse(grepl("CL|DL|SV", merged$fusil_bin_code) & !grepl("outlier", merged$fusil_bin_code) & !is.na(merged$fusil_bin_code), "yes", "no")

print("Proportion of lowest decile that are essential/required for development or viability")
print(round(nrow(merged[merged$oe_mis_upper_decile == "1" & merged$essential == "yes",]) / nrow(merged[merged$oe_mis_upper_decile == "1",]), 3))
print("Proportion of all genes that are essential/required for development or viability")
print(round(nrow(merged[merged$essential == "yes",]) / nrow(merged), 3))

# how many nuclear genes have dominant missense that can underlie developmental disorders in humans?
ddd <- read.delim(file = '../required_files/other_annotations/DDG2P_13_4_2023.csv', header = TRUE, sep = ",")

# filter to genes with dominant missense underlying DD, and removed low confidence 
ddd <- ddd[ddd$confidence.category != "limited" & grepl("monoallelic", ddd$allelic.requirement) & grepl("missense", ddd$variant.consequence),]
ddd$gene <- ddd$gene.symbol
# collapse since multiple rows for some genes
ddd <- as.data.frame(ddd[,c("gene", "organ.specificity.list")] %>% 
                       group_by(gene) %>% 
                       summarise(across(everything(), str_c, collapse=","))) 
ddd$ddd <- "yes"
merged_two <- merge(merged, ddd, by = "gene", all.x = TRUE)

print("Proportion of genes associated with dominant missense developmental disorders that are in lowest missense decile")
print(round(nrow(merged_two[merged_two$oe_mis_upper_decile == "1" & merged_two$ddd == "yes" & !is.na(merged_two$ddd),]) / nrow(merged_two[merged_two$ddd == "yes" & !is.na(merged_two$ddd),]), 3))
print("Proportion of all genes associated with dominant missense developmental disorders")
print(round(nrow(merged_two[merged_two$ddd == "yes" & !is.na(merged_two$ddd),]) / nrow(merged_two), 3))

sink()
