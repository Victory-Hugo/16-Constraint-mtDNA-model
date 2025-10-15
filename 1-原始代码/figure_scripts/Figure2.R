library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(png)

# FIGURE 2

# Figure 2a - plot the missense observed:expected ratio and 90% confidence interval for each protein, order by complex

file <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$complex <- ifelse(grepl("MT-A", file$locus), "V", 
                       ifelse(grepl("MT-CO", file$locus), "IV", 
                              ifelse(grepl("MT-CY", file$locus), "III", 
                                     ifelse(grepl("MT-N", file$locus), "I", "error"))))

plotA <- ggplot(file[file$consequence == "missense", ], aes(y = fct_rev(locus), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.85, colour = "#777575") +
  geom_point(shape = 18, size = 5.5, color = '#377eb8') +
  labs(x = "Missense observed:expected ratio") + 
  xlim(0, 1) +
  paper_theme +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(rows = vars(complex), scales = "free", space = 'free')

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[file$consequence == "missense", c("locus", "consequence", "variant_count", "obs.exp", "lower_CI", "upper_CI", "complex")],
            file = 'final_figures_source_data/Figure2a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 2b - linear and 3D figure of regional constraint using MT-ND1 as an example

load("extended_data_figures/rdata/MT-ND1.Rdata")  # linear
linear <- plot +
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"),
        axis.title.x = element_blank())

chimera <- readPNG("figures/Figure2b_ND1.png")  # 3D
chimera_fig <- ggplot() + 
  background_image(chimera) +
  theme(plot.margin = unit(c(0.0, 1.2, 0.1, 1.2), "cm"))

plotB <- ggarrange(linear, chimera_fig, nrow = 2, ncol = 1, heights = c(0.2, 0.8))


# Figure 2c - odds ratio enrichment of pathogenic vs benign missense variants

# read in file and subset to missense (most severe)
file <- read.delim(file = '../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file <- file[grepl("missense", file$consequence) & !grepl("stop_gain|start_lost|stop_lost|incomplete_terminal", file$consequence),]
file$in_rc <- factor(ifelse(file$in_rc == "no" & !is.na(file$min_distance_to_rc) & (as.numeric(file$min_distance_to_rc) <= 6), "proximal", as.character(file$in_rc)), 
                     levels = c("yes", "proximal", "no"), labels = c("Within", "Proximal", "Outside"))

# use ClinVar pathogenic and likely pathogenic variants and MITOMAP confirmed disease variants as pathogenic, and ClinVar benign variants as benign
file$group <- factor(ifelse(grepl("Cfrm", file$mitomap_status) | grepl("athogenic", file$clinvar_interp), "Pathogenic",
                     ifelse(grepl("Benign", file$clinvar_interp), "Benign", "neither")),
                     levels = c("Pathogenic", "Benign", "neither"))

odds_ratio <- function(df) {
  output <- data.frame(matrix(vector(), ncol = 7))
  colnames(output) <-c("bin", "value", "se", "lower_CI", "upper_CI", "number_pathogenic", "number_benign")
  output[1,] <- c(NA, NA, NA, NA, NA, NA, NA)
  
  n_pathogenic <- nrow(df[df$group == "Pathogenic", ]) 
  n_benign <- nrow(df[df$group == "Benign", ])
  
  for(group in list("Within", "Proximal", "Outside") ){
    a = nrow(df[df$in_rc == group & df$group == "Pathogenic", ])
    b = nrow(df[df$in_rc == group & df$group == "Benign", ])
    c = n_pathogenic - a
    d = n_benign - b
    OR <- round((a / b) / (c / d), digits = 4)
    OR_se <- round(sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d)), digits = 4)
    OR_lowerCI <- round(exp(log(OR) - (1.96 * OR_se)), digits = 4)
    OR_upperCI <- round(exp(log(OR) + (1.96 * OR_se)), digits = 4)
    
    row <- as.data.frame(t(c(group, OR, OR_se, OR_lowerCI, OR_upperCI, a, b)))
    colnames(row) <- c("bin", "value", "se", "lower_CI", "upper_CI", "number_pathogenic", "number_benign")
    output <- rbind(output, row)
  }
  output <- output[!is.na(output$bin), ]
  return(output)
}

or <- odds_ratio(file)

fisher.test(data.frame(
  "pathogenic" = as.numeric(or[, c("number_pathogenic")]),
  "benign" = as.numeric(or[, c("number_benign")]),
  row.names = c("Within", "Proximal", "Outside"),
  stringsAsFactors = FALSE
))

plotC <- ggplot(or, aes(x = bin, y = as.numeric(value), fill = bin)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = as.numeric(lower_CI), ymax = as.numeric(upper_CI)), width = 0.5, position = position_dodge(.9), size = 0.5) +
  scale_y_sqrt(expand = c(0, 0), breaks = c(0, 1, 5, 10, 20, 30, 40, 50)) + 
  paper_theme +
  scale_x_discrete(labels = c("Outside", "Proximal", "Within")) +
  scale_fill_manual(values = c("#377eb8", "#926281", "#ff4040"), guide = FALSE) +
  labs(y = "OR (pathogenic vs benign)", x = "\nRegional missense constraint") +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))

as.data.frame(file[file$group != "neither",] %>% group_by(group) %>% summarize(n = n()))
as.data.frame(file[file$group != "neither",] %>% group_by(group, in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))

write.table(or[, c("bin", "value", "lower_CI", "upper_CI")],
            col.names = c("group", "OR_value", "lower_CI", "upper_CI"),
            file = 'final_figures_source_data/Figure2c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 2d - plot regional missense constraint data as stacked bar to show proportions for clinical classifications

mcri <- read.delim(file = '../required_files/other_annotations/curated_mtDNA_variants_combined_final_copy.txt', header = TRUE, sep = "\t")
mcri$var <- paste(mcri$REF, mcri$POS, mcri$ALT)
file <- read.delim(file='../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file$var <- paste(file$REF, file$POS, file$ALT)

# annotate, and restrict to missense as most severe consequence 
file <- merge(file[grepl("missense", file$consequence) & !grepl("stop_gain|start_lost|stop_lost|incomplete_terminal", file$consequence),], mcri[, c("var", "Classification")], by = "var", all.x = T) 
file$in_rc <- factor(ifelse(file$in_rc == "no" & !is.na(file$min_distance_to_rc) & (as.numeric(file$min_distance_to_rc) <= 6), "proximal", as.character(file$in_rc)), 
                     levels = c("yes", "proximal", "no"), labels = c("Within", "Proximal", "Outside"))
# relabel classification
file$group <- factor(ifelse(grepl("Class 5|Class 4", file$Classification), "Pathogenic & Likely pathogenic",
                            ifelse(grepl("Class 1|Class 2", file$Classification), "Benign & Likely Benign", 
                                   paste("VUS ", as.character(file$Classification), sep = ""))),
                     levels = c("Benign & Likely Benign", "VUS Class 3c", "VUS Class 3b", "VUS Class 3a", "Pathogenic & Likely pathogenic"),
                     labels = c("Benign\n& Likely\nBenign", "VUS of low\nclinical\nsignificance", "VUS\n", 
                                "VUS of high\nclinical\nsignificance", "Pathogenic\n& Likely\npathogenic"))

file[!is.na(file$group),] %>% group_by(group) %>% summarize(n())

plotD <- ggplot(file[!is.na(file$Classification), ], aes(group, fill = in_rc)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(y = "Proportion of curated missense", x = "Missense variant classification") + 
  paper_theme +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.25, 0.2, 0.25, 0.4), "cm")) +
  scale_fill_manual(values = c("#ff4040","#926281", "#377eb8"), 
                    labels = c("Within", "Proximal", "Outside"), name = "Regional\nmissense\nconstraint", guide = FALSE) 

# remove new lines for writing
file$group <- gsub("\n", " ", file$group)

write.table(file[!is.na(file$Classification), c("POS", "REF", "ALT", "group", "in_rc")],
            col.names = c("position", "ref", "alt", "classification_group", "in_rc"),
            file = 'final_figures_source_data/Figure2d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# compile panel
ggarrange(
  ggarrange(plotA, plotB, nrow = 1, ncol = 2, labels = c("a", "b"), font.label = list(size = 10)),
  ggarrange(plotC, plotD, nrow = 1, ncol = 2, widths = c(0.4, 0.6), labels = c("c", "d"), font.label = list(size = 10)), 
  nrow = 2, ncol = 1, heights = c(1, 0.65))

ggsave("figures/Figure2.jpeg", width = 180, height = 140, dpi = 600, units = c("mm"))


