library(dplyr)
library(ggExtra)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(scales)
library(tidyr)

# SUPPLEMENTARY FIGURE 8

# Figure S8a - show that sites with the same AF can have different mutability in gnomAD

# read in annotations
all <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# for the kmers used to compute MLC, compute the sum mutational likelihood and sum allele counts
mlc <- read.delim(file = "../output_files/local_constraint/kmers_local_constraint.txt", header = TRUE, sep = '\t')

# specify variant types used for MLC - only these will be included
ok_list <- c('non_coding_transcript_exon_variant','missense_variant', 'intergenic_variant', 
             'missense_variant,synonymous_variant', 'missense_variant,missense_variant', 'stop_retained_variant',
             'non_coding_transcript_exon_variant,non_coding_transcript_exon_variant', 'synonymous_variant,missense_variant')

# set as numeric
all$Likelihood <- as.numeric(all$Likelihood)
all$gnomad_ac_het <- as.numeric(all$gnomad_ac_het)
all$gnomad_max_hl <- as.numeric(all$gnomad_max_hl)

# function to calculate the sum mutational likelihood and sum allele counts for each kmer
df <- data.frame()
for(i in 1:nrow(mlc)){
  # handle across the artificial break
  if(as.numeric(mlc$X..start[i]) < as.numeric(mlc$end[i])){
    kmer <- seq(as.numeric(mlc$X..start[i]), as.numeric(mlc$end[i]))
  }
  else{
    kmer <- c(seq(1, as.numeric(mlc$end[i])), seq(as.numeric(mlc$X..start[i]), 16569))
  }
  print(paste(mlc$X..start[i], mlc$end[i]))
  sum_lr = 0
  sum_gnomad_het = 0
  sum_obs = 0
  for(pos in kmer){
    sum_lr = sum_lr + sum(all[all$POS == pos & (all$consequence %in% ok_list), c("Likelihood")])
    sum_gnomad_het = sum_gnomad_het + sum(all[all$POS == pos & (all$consequence %in% ok_list), c("gnomad_ac_het")])
    sum_obs = sum_obs + sum(all[all$POS == pos & (all$consequence %in% ok_list), c("gnomad_max_hl")])
  }
  df <- rbind(df, c(mlc$X..start[i], mlc$end[i], sum_lr, sum_gnomad_het,  sum_obs))
}

colnames(df) <- c("start", "end", "sum_lr", "sum_gnomad_het", "sum_obs")

# merge with existing data
df$range <- paste(df$start, "-", df$end)
mlc$range <- paste(mlc$X..start, "-", mlc$end)
for_plot <- merge(mlc[,c("range", "upper_CI", "variant_count", "loci", "protein_sites", "observed", "RNA_sites")], df, by = "range", all.x = TRUE)
for_plot$check <- for_plot$sum_obs/for_plot$observed # note handles variants in two genes slightly differently but not needed here

# plot across the RNA genes which have the same number of variants in each kmer and highlight a highly constrained RNA modification and RNA sites in the functional study
plotA <- ggplot(for_plot[grepl("MT-RNR|MT-T", for_plot$loci) & !grepl(",", for_plot$loci) & for_plot$sum_gnomad_het < 50,], aes(x = sum_lr, y = sum_gnomad_het)) +
  geom_point(aes(color = upper_CI), size = 0.65, position = "jitter") + 
  paper_theme +
  ylim(0, 50) +
  scale_color_continuous(limits = c(0, 1)) +
  labs(x = "Mutation likelihood", y = "Population AC het gnomAD", color = "OEUF") +
  geom_point(data = for_plot[grepl("MT-RNR2", for_plot$loci) & grepl("2_-O-methylation", for_plot$RNA_sites) &
                               for_plot$start <= 3039 & for_plot$end >= 3040,],
             color = "red", size = 0.6)  +
  geom_point(data = for_plot[grepl("MT-RNR2", for_plot$loci) & (for_plot$start <= 3047 & for_plot$end >= 3075),], 
             color = "green", size = 0.7)  +
  theme(legend.key.width = unit(0.5, 'cm'),
        legend.text = element_text(size = 7), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.position = "left")


# Figure S8b - show that sites with the same AF can have different mutability in a different dataset (UKB)

# collate annotations
file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$var <- paste(file$POS, file$REF, ">", file$ALT, sep = "")

# read in source data from Nature Communications manuscript
natcom <- read.delim(file = '../required_files/other_annotations/PMID37777527_Fig1_source_data.txt', header = TRUE, sep = "\t")
natcom$var <- paste(natcom$POS, natcom$REF, ">", natcom$ALT, sep = "")

# merge (specifically counts)
all <- merge(file, natcom[,c("var", "count_het")], by = "var", all.x = T)

# set missing values to 0 and compute af
all[is.na(all$count_het), c("count_het")] <- "0"
all$het_af_ukb <- as.numeric(all$count_het)/194871  # number of individuals included in Nat Comms

# set as numeric
all$Likelihood <- as.numeric(all$Likelihood)
all$count_het <- as.numeric(all$count_het)
all$gnomad_max_hl <- as.numeric(all$gnomad_max_hl)

# function to calculate the sum mutational likelihood and sum allele counts for each kmer
df <- data.frame()
for(i in 1:nrow(mlc)){
  # handle across the artificial break
  if(as.numeric(mlc$X..start[i]) < as.numeric(mlc$end[i])){
    kmer <- seq(as.numeric(mlc$X..start[i]), as.numeric(mlc$end[i]))
  }
  else{
    kmer <- c(seq(1, as.numeric(mlc$end[i])), seq(as.numeric(mlc$X..start[i]), 16569))
  }
  print(paste(mlc$X..start[i], mlc$end[i]))
  sum_lr = 0
  sum_ukb_het = 0
  for(pos in kmer){
    sum_lr = sum_lr + sum(all[all$POS == pos & (all$consequence %in% ok_list), c("Likelihood")])
    sum_ukb_het = sum_ukb_het + sum(all[all$POS == pos & (all$consequence %in% ok_list), c("count_het")])
  }
  df <- rbind(df, c(mlc$X..start[i], mlc$end[i], sum_lr, sum_ukb_het))
}

colnames(df) <- c("start", "end", "sum_lr", "sum_ukb_het")

# merge with existing data
df$range <- paste(df$start, "-", df$end)
mlc$range <- paste(mlc$X..start, "-", mlc$end)
for_plot <- merge(mlc[,c("range", "upper_CI", "variant_count", "loci", "protein_sites", "observed", "RNA_sites")], df, by = "range", all.x = TRUE)

# plot across the RNA genes which have the same number of variants in each kmer and highlight a highly constrained RNA modification and RNA sites in the functional study
plotB <- ggplot(for_plot[grepl("MT-RNR|MT-T", for_plot$loci) & !grepl(",", for_plot$loci),], aes(x = sum_lr, y = sum_ukb_het)) +
  geom_point(aes(color = upper_CI), size = 0.65, position = "jitter") + 
  paper_theme +
  ylim(0, 200) +
  scale_color_continuous(limits = c(0, 1)) +
  labs(x = "Mutation likelihood", y = "Population AC het UKB", color = "OEUF\n") +
  geom_point(data = for_plot[grepl("MT-RNR2", for_plot$loci) & grepl("2_-O-methylation", for_plot$RNA_sites) &
                               for_plot$start <= 3039 & for_plot$end >=3040,],
             color = "red", size = 0.6)  +
  geom_point(data = for_plot[grepl("MT-RNR2", for_plot$loci) & (for_plot$start <= 3047 & for_plot$end >= 3075),], 
             color = "green", size = 0.7)  +
  theme(legend.position="none")


# collate plot
ggarrange(plotA, plotB, ncol = 2, widths = c(0.6, 0.45), labels = c("a", "b"), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS8.jpeg", width = 140, height = 50, dpi = 600, units = c("mm"))


