library(ggplot2)
library(ggpubr)

# SUPPLEMENTARY DATA 2

# Figure S2a - plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD

file <- read.delim(file = '../output_files/calibration/loci_obs_vs_scores.txt', header = TRUE, sep = "\t")

plotA <- ggplot(data = file, aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 1) + 
  labs(x = "Mutation likelihood", y = "Observed neutral variation") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.8, .95),
        legend.key.size = unit(0.8, "mm"), 
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0.55, 0.25, 0.25, 0.25), "cm")) + 
  scale_color_discrete(labels = c("G>A and T>C", "All other variants")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(700, 650), digits = 2, show.legend = FALSE, size = 2.25) +
  geom_smooth(method = 'lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))


# FigureS2b - plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD, in tRNA genes only

plotB <- ggplot(data = file[grepl("MT-T", file$symbol), ], aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 1) +
  labs(x = "Mutation likelihood", y = "Observed neutral variation in tRNAs") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.8, .95),
        legend.key.size = unit(0.8, "mm"), 
        legend.text = element_text(size = 7),
        axis.title.y = element_text(size = 7), 
        plot.margin = unit(c(0.55, 0.25, 0.25, 0.25), "cm")) + 
  scale_color_discrete(labels = c("G>A and T>C", "All other variants")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(28, 26), digits = 2, show.legend = FALSE, size = 2.25) +
  geom_smooth(method='lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))


# Figure S2c - plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD, in the OriB-OriH region

ori_file <- read.delim(file = '../output_files/calibration/loci_obs_vs_scores_ori.txt', header = TRUE, sep = "\t")

plotC <- ggplot(data = ori_file, aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 1) +
  labs(x = "Mutation likelihood", y = "Observed neutral variation in OriB-OriH") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.8, .95),
        legend.key.size = unit(0.8, "mm"), 
        legend.text = element_text(size = 7),
        axis.title.y = element_text(size = 7), 
        plot.margin = unit(c(0.55, 0.25, 0.25, 0.25), "cm")) + 
  ylim(c(1, 110)) + 
  scale_color_discrete(labels = c("G>A and T>C", "All other variants")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(110, 102.5), digits = 2, show.legend = FALSE, size = 2.25) +
  geom_smooth(method = 'lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))


# collate panel
ggarrange(plotA, plotB, plotC, ncol = 3, nrow = 1, labels = c("a", "b", "c"), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS2.jpeg", width = 180, height = 60, dpi = 600, units = c("mm"))  

