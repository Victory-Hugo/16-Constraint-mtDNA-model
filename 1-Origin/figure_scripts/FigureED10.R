library(ggplot2)
library(ggExtra)
library(ggpubr)
library(scales)

# EXTENDED DATA FIGURE 10

# Figure ED10a - scatter plot of MLC vs UKB heteroplasmy frequency

# first, create df with all annotations
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$var <- paste(scores$POS, scores$REF, ">", scores$ALT, sep = "")
file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$var <- paste(file$POS, file$REF, ">", file$ALT, sep = "")
all_scores <- merge(file, scores[, c("var", "MLC_pos_score", "MLC_var_score", "pctRank_mean_OEUF")], by = "var", all.x = T)

# read in source data from Nature Communications manuscript
natcom <- read.delim(file = '../required_files/other_annotations/PMID37777527_Fig1_source_data.txt', header = TRUE, sep = "\t")
natcom$var <- paste(natcom$POS, natcom$REF, ">", natcom$ALT, sep = "")

# merge (specifically heteroplasmy count)
all <- merge(all_scores, natcom[,c("var", "count_het")], by = "var", all.x = T)

# set missing values to 0 and compute
all[is.na(all$count_het), c("count_het")] <- "0"
all$het_af_ukb <- as.numeric(all$count_het)/194871  # number of individuals included in Nat Comms

# set this up for plotting
mysqrt_trans <- function() { #this is a hack to enable y=0 tick mark on plot per https://gist.github.com/DarwinAwardWinner/21652acf017880c271e95cc2e35574f4
  trans_new("mysqrt", 
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}


# now plot
p_track1.1 <- ggplot(all[all$het_af_ukb > 0 & (grepl("intergenic|non_coding_transcript", all$consequence) | all$consequence == "synonymous_variant" | all$consequence == "synonymous_variant,synonymous_variant"),], 
                     aes(x = as.numeric(MLC_var_score), y = as.numeric(het_af_ukb))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.75, label.y = 0.1, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'UKB allele frequency', title = "\nAll synonymous, RNA & non-coding SNVs") +
  geom_smooth(method = 'lm', size = 0.7) +
  paper_theme +
  scale_y_continuous(trans = "sqrt", breaks = c(0.0001, 0.001, 0.005, 0.01), 
                     labels = label_number(accuracy = 0.0001), limits = c(0, 0.01)) +
  theme(plot.title = element_text(hjust = 0.5, size = 6.5)) +
  xlim(0, 1)

p_track1.2 <- ggplot(all[all$het_af_ukb > 0 & grepl("intergenic|non_coding_transcript", all$consequence),], 
                     aes(x = as.numeric(MLC_var_score), y = as.numeric(het_af_ukb))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.75, label.y = 0.1, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'UKB allele frequency', title = "\nAll RNA & non-coding SNVs") +
  geom_smooth(method = 'lm', size = 0.7) +
  paper_theme +
  scale_y_continuous(trans = "sqrt", breaks = c(0.0001, 0.001, 0.005, 0.01), 
                     labels = label_number(accuracy = 0.0001), limits = c(0, 0.01)) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 6.5)) +
  xlim(0, 1)

p_track2 <- ggplot(all[all$het_af_ukb > 0 & grepl("missense|stop_gained", all$consequence),], 
                   aes(x = as.numeric(MLC_var_score), y = as.numeric(het_af_ukb))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.75, label.y = 0.1, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'UKB allele frequency', title = "\nAll missense & stop gain SNVs") +
  geom_smooth(method = 'lm', size = 0.7) +
  paper_theme +
  scale_y_continuous(trans = "sqrt", breaks = c(0.0001, 0.001, 0.005, 0.01), 
                     labels = label_number(accuracy = 0.0001), limits = c(0, 0.01)) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 6.5)) +
  xlim(0, 1)


# Figure ED10b - scatter plot of MLC vs UKB heteroplasmy phyloP

# plot
r_track1.1 <- ggplot(all[all$het_af_ukb > 0 & (grepl("intergenic|non_coding_transcript", all$consequence) | all$consequence == "synonymous_variant" | all$consequence == "synonymous_variant,synonymous_variant"),], 
                     aes(x = as.numeric(MLC_var_score), y = as.numeric(phyloP_score))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.8, label.y = -19, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'phyloP score') +
  paper_theme +
  geom_smooth(method = 'lm', size = 0.7) +
  scale_y_continuous(limits = c(-20, 10), labels = c("     -20", "      -10", "      0", "     10"))

r_track1.2 <- ggplot(all[all$het_af_ukb > 0 & grepl("intergenic|non_coding_transcript", all$consequence),], 
                     aes(x = as.numeric(MLC_var_score), y = as.numeric(phyloP_score))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.8, label.y = -19, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'phyloP score') +
  paper_theme +
  geom_smooth(method = 'lm', size = 0.7) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001), limits = c(-20, 10))

r_track2 <- ggplot(all[all$het_af_ukb > 0 & grepl("missense|stop_gained", all$consequence),], 
                   aes(x = as.numeric(MLC_var_score), y = as.numeric(phyloP_score))) +
  geom_point(size = 0.7) +
  stat_cor(method = "pearson", label.x = 0.8, label.y = -19, digits = 1, show.legend = FALSE, size = 2, aes(label = ..rr.label..)) +  # the label only shows the R^2 value
  labs(x = 'MLC score', y = 'phyloP score') +
  paper_theme +
  geom_smooth(method = 'lm', size = 0.7) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001), limits = c(-20, 10))

write.table(all[order(all$POS),c("POS", "REF", "ALT", "consequence", "het_af_ukb", "MLC_var_score", "phyloP_score")],
            col.names = c("pos", "ref", "alt", "consequence", "het_af_ukb", "MLC_score", "phyloP"),
            file = 'final_figures_source_data/FigureED10a-b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# collate plot
ggarrange(
  ggarrange(ggExtra::ggMarginal(p_track1.1, type = "histogram", margins = "y", yparams = list(binwidth = 0.0025)),
            ggExtra::ggMarginal(p_track1.2, type = "histogram", margins = "y", yparams = list(binwidth = 0.0025)),
            ggExtra::ggMarginal(p_track2, type = "histogram", margins = "y", yparams = list(binwidth = 0.0025)), nrow = 1, widths = c(1.05, 0.85, 0.85)),
  ggarrange(ggExtra::ggMarginal(r_track1.1, type = "histogram", margins = "y", yparams = list(binwidth = 0.5)),
            ggExtra::ggMarginal(r_track1.2, type = "histogram", margins = "y", yparams = list(binwidth = 0.5)),
            ggExtra::ggMarginal(r_track2, type = "histogram", margins = "y", yparams = list(binwidth = 0.5)), nrow = 1, widths = c(1.05, 0.85, 0.85)),
  nrow = 2, heights = c(0.475, 0.4), labels = c("a", "b"), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED10.jpeg", width = 180, height = 105, dpi = 600, units = c("mm"))


