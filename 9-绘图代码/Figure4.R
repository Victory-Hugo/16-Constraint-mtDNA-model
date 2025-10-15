library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(png)
library(stringr)
library(readr)

# FIGURE 4

# Figure 4a - plot of local constraint across mtDNA

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")

svg <- readPNG("figures/Figure4a_top.png") 
plotA_top <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.05, 1.25, 0.05, 0.9), "cm"),
        panel.background = element_blank())

plotA1 <- ggplot(scores, aes(x = as.numeric(POS), y = as.numeric(MLC_pos_score))) +
  geom_line(aes(color = pctRank_mean_OEUF)) + 
  scale_color_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 16569)) + 
  ylim(0, 1) + 
  labs(color = "MLC score", x = "Position", y = "MLC score") + 
  paper_theme +
  theme(plot.margin = unit(c(0, 0.25, 0, 0.25), "cm"),
        axis.title.x = element_text(vjust = 5),
        legend.key.height = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "dark grey", size = 0.2)

svg <- readPNG("figures/Figure4a_mtDNA.png") 
plotA2 <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0, 1.73, 0.25, 0.8), "cm"))


# Figure 4b - stacked bar plot to show the proportion of each locus type across local constraint score bins

# divide into four bins
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))
# assign biotype
scores$biotype <- ifelse(grepl("intergenic",scores$consequence), "Non-coding",
                         ifelse(grepl("MT-R",scores$symbol), "rRNA",
                                ifelse(grepl("MT-T",scores$symbol), "tRNA",
                                       ifelse(grepl("MT-A|MT-C|MT-N",scores$symbol), "Protein", "Error"))))
# biotype colors
mycolors = c('#ffcc00', '#377eb8', '#984ea3', '#ff7f00')

plotB <- ggplot(scores[!duplicated(scores$POS),], aes(rank_bin, fill = biotype, color = biotype)) + 
  geom_bar(position = "fill", color = "black", size = 0.25, width = 0.9) +
  labs(x = 'MLC score quartile', y = 'Proportion', fill = 'Locus type') +
  paper_theme + 
  scale_fill_manual(values = mycolors) +
  scale_colour_manual(values = mycolors, guide = FALSE) +
  theme(plot.margin = unit(c(0.3, 0.25, 0.25, 0.5), "cm"))

write.table(unique(scores[, c("POS", "MLC_pos_score", "rank_bin", "biotype")]),
            file = 'final_figures_source_data/Figure4a-b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 4c - bar plot to show the odds ratio enrichment of pathogenic vs benign variants across local constraint score bins

# divide into four bins
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))

# assign pathogenic or benign status, including pathogenic and likely pathogenic in clinvar
scores$group <- ifelse(grepl("Cfrm", scores$mitomap_status) | grepl("athogenic", scores$clinvar_interp), "Pathogenic",
                        ifelse(grepl("Benign", scores$clinvar_interp), "Benign", "neither")) 

odds_ratio <- function(df) {
  output <- data.frame(matrix(vector(), ncol = 7))
  colnames(output) <-c("bin", "value", "se", "lower_CI", "upper_CI", "number_pathogenic", "number_benign")
  output[1,] <- c(NA, NA, NA, NA, NA, NA, NA)
  
  n_pathogenic <- nrow(df[df$group == "Pathogenic", ]) 
  n_benign <- nrow(df[df$group == "Benign", ])
  
  for(group in list("0.0-0.25", "0.25-0.50", "0.50-0.75", "0.75-1.0") ){
    a = nrow(df[df$rank_bin == group & df$group == "Pathogenic", ])
    b = nrow(df[df$rank_bin == group & df$group == "Benign", ])
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

# restrict to missense (most severe) and RNA base changes (ie in genes) - position score same as variant score for these
or <- odds_ratio(scores[grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),])

fisher.test(data.frame(
  "pathogenic" = as.numeric(or[, c("number_pathogenic")]),
  "benign" = as.numeric(or[, c("number_benign")]),
  row.names = c("0.0-0.25", "0.25-0.50", "0.50-0.75", "0.75-1.0"),
  stringsAsFactors = FALSE
))

plotC <- ggplot(or, aes(x = bin, y = as.numeric(value), fill = bin)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = as.numeric(lower_CI), ymax = as.numeric(upper_CI)), width = 0.5, position = position_dodge(.9), size = 0.5) +
  scale_y_sqrt(expand = c(0, 0), breaks = c(0, 1, 5, 10)) + 
  paper_theme +
  labs(y = "OR (pathogenic vs benign)", x = "MLC score quartile") +
  scale_fill_manual(values = c("#542eff", "#cfb1ff", "#ffbfaa", "#ff4124"), guide = FALSE) +
  theme(plot.margin = unit(c(0.3, 0.25, 0.25, 0.5), "cm"))

write.table(or[, c("bin", "value", "lower_CI", "upper_CI")],
            col.names = c("group", "OR_value", "lower_CI", "upper_CI"),
            file = 'final_figures_source_data/Figure4c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# relevant statistics for manuscript
sum(as.numeric(or$number_pathogenic))
sum(as.numeric(or$number_benign))


# compile figure panel
ggarrange(
  ggarrange(plotA_top, plotA1, NULL, plotA2, nrow = 4, labels = c("a", "", "", "", ""), heights = c(1.3, 1, -0.05, 0.6), font.label = list(size = 10)),
  ggarrange(plotB, plotC, nrow = 1, ncol = 2, labels = c("b", "c"), widths = c(0.575, 0.425), font.label = list(size = 10)), 
  nrow = 2, ncol = 1, heights = c(1.4, 0.65)) 

ggsave("figures/Figure4.jpeg", width = 180, height = 140, dpi = 600, units = c("mm"))
