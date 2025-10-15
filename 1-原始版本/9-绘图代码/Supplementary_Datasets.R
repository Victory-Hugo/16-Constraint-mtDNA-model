# Supplementary Dataset 1 - constraint metrics for each gene

file <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
file$consequence <- ifelse(file$consequence == "SNV", "RNA_variant", as.character(file$consequence))
write.table(file[, c("locus", "start", "end", "consequence", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_1.tsv', col.names = c("Symbol", "Start_position", "End_position", "Consequence", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 2 - regional constraint intervals

file <- read.delim(file = '../output_files/regional_constraint/final_regional_constraint_intervals.txt', header = TRUE, sep= "\t")
file$obs_max_het <- round(file$obs_max_het, 3)
file$exp_max_het <- round(file$exp_max_het, 3)
file$ratio_oe <- round(file$ratio_oe, 3)
write.table(file[, c("locus", "start", "end", "protein_pos_start", "protein_pos_end", "obs_max_het", "exp_max_het", "ratio_oe", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_2.tsv', col.names = c("Symbol", "Start_position", "End_position", "Protein_residue_start", "Protein_residue_end", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 3 - curated missense from VCGS

mcri <- read.delim(file='../required_files/other_annotations/curated_mtDNA_variants_combined_final_copy.txt', header = TRUE, sep = "\t")
mcri$var <- paste(mcri$REF, mcri$POS, mcri$ALT)
mcri$Classification <- factor(ifelse(grepl("Class 5|Class 4", mcri$Classification), "Pathogenic & Likely pathogenic",
                            ifelse(grepl("Class 1|Class 2", mcri$Classification), "Benign & Likely Benign", paste("VUS ", as.character(mcri$Classification), sep = ""))),
                     levels = c("Benign & Likely Benign", "VUS Class 3c", "VUS Class 3b", "VUS Class 3a", "Pathogenic & Likely pathogenic"),
                     labels = c("Benign & Likely Benign", "VUS of low clinical significance", "VUS ", 
                                "VUS of high clinical significance", "Pathogenic & Likely pathogenic"))
# annotate with consequences
file <- read.delim(file='../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file$var <- paste(file$REF, file$POS, file$ALT)

# restrict to missense, as most severe consequence 
mcri <- merge(mcri, file[, c("var", "consequence")], by = "var", all.x = T) 
mcri <- mcri[grepl("missense", mcri$consequence) & !grepl("stop_gain|start_lost|stop_lost|incomplete_terminal", mcri$consequence),]

write.table(mcri[order(mcri$POS), c("POS", "REF", "ALT", "Classification")], 
            file = 'supplementary_datasets/supplementary_dataset_3.tsv', col.names = c("Position", "Reference", "Alternate", "Classification_group"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 4 - Constraint metrics for each position within the tRNA secondary structure

file <- read.delim(file = '../output_files/oe/tRNA_position_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
# manually reorder
write.table(file[c(2:23,75,24:48,76,49:52,71:72,53:59,73:74,60:70), c("tRNA_position", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_4.tsv', col.names = c("tRNA_position", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 5 - Constraint metrics for non-coding elements

file <- read.delim(file = '../output_files/oe/noncoding_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
write.table(file[, c("locus", "description", "start", "end", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_5.tsv', col.names = c("Locus", "Description", "Start_position", "End_position", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 6 - MLC scores for each base position

file <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
file$MLC_pos_score <- round(file$MLC_pos_score, 5)
write.table(unique(file[, c("POS", "MLC_pos_score")]), 
            file = 'supplementary_datasets/supplementary_dataset_6.tsv', col.names = c("Position", "MLC_pos_score"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 7 - MLC scores for every SNV

file <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
file$MLC_var_score <- round(file$MLC_var_score, 5)
write.table(file[, c("POS", "REF", "ALT", "consequence", "MLC_var_score")], 
            file = 'supplementary_datasets/supplementary_dataset_7.tsv', col.names = c("Position", "Reference", "Alternate", "Consequence", "MLC_score"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 8 - catalog of neutral variants used to validate the model
file <- read.delim(file = '../output_files/calibration/neutral_variants_used.txt', header = TRUE, sep = "\t")
write.table(unique(file[order(file$position), c("variant", "source")]), 
            file = 'supplementary_datasets/supplementary_dataset_8.tsv', col.names = c("Variant", "Criteria"), row.names = FALSE, sep = '\t', quote = FALSE)


# Supplementary Dataset 9 - curated confirmed MITOMAP variants by plasmy

file <- read.delim(file = '../required_files/other_annotations/cnfm_by_status_curated.txt', header = TRUE, sep = "\t")
file$curated_homoplasmy_report <- gsub(" ", "", file$curated_homoplasmy_report)
write.table(file[grepl("missense|transcript", file$consequence) & !grepl("gain|lost|terminal", file$consequence), c("var", "symbol", "mitomap_plasmy", "curated_status", "curated_homoplasmy_report")], 
            file = 'supplementary_datasets/supplementary_dataset_9.tsv', col.names = c("Variant", "Symbol", "MITOMAP_plasmy", "Curated_status", "Curated_homoplasmy_report"), row.names = FALSE, sep = '\t', quote = FALSE)
