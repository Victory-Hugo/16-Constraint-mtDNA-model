library(ggplot2)
library(ggpubr)
library(png)

# EXTENDED DATA FIGURE 8

# Figure ED8a - plot MLC position scores in the control region

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")

# schematic with annotated loci in control region 
svg <- readPNG("extended_data_figures/FigureED8a-left.png") 
plot_top_left <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.3, 0.25, 0.2, 0.9), "cm"))

svg <- readPNG("extended_data_figures/FigureED8a-right.png") 
plot_top_right <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.3, 1.75, 0.2, 0), "cm"))

plot_top <- ggarrange(plot_top_left, plot_top_right, ncol = 2, widths = c(1, 1))

plot1 <- ggplot(scores[scores$POS > 16023,], aes(x = as.numeric(POS), y = as.numeric(MLC_pos_score))) +
  geom_line(aes(color = pctRank_mean_OEUF)) + 
  scale_color_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(16024, 16100, 16200, 16300, 16400, 16500)) + 
  ylim(0, 1) + 
  labs(fill = "MLC score", x = "Position", y = "MLC score") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0, 0.25, 0.5), "cm")) +
  guides(color = FALSE) 

plot2 <- ggplot(scores[scores$POS < 578,], aes(x = as.numeric(POS), y = as.numeric(MLC_pos_score))) +
  geom_line(aes(color = pctRank_mean_OEUF)) + 
  scale_color_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300, 400, 500, 577)) + 
  ylim(0, 1) + 
  labs(color = "MLC score", x = "Position", y = "MLC score") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(0.5, 0.5, 0.25, 0), "cm"),
        legend.key.height = unit(0.4, 'cm')) 

plotA <- ggarrange(plot_top, NULL,
                   ggarrange(plot1, plot2, ncol = 2, nrow = 1, widths = c(5.5, 6)),
                   ncol = 1, nrow = 3, heights = c(1, -0.05, 1.75))

write.table(unique(scores[(scores$POS > 16023) | (scores$POS < 578), c("POS", "MLC_pos_score")]), 
            col.names = c("pos", "MLC_score"),
            file = 'final_figures_source_data/FigureED8a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED8b-d - plot allele frequency data across control region

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# gnomAD
plot1 <- ggplot(file[file$POS > 16023,], aes(x = as.numeric(POS), y = as.numeric(gnomad_af_hom))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(16024, 16100, 16200, 16300, 16400, 16500)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "gnomAD AF") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0, 0.25, 0.5), "cm")) +
  guides(color = FALSE) 

plot2 <- ggplot(file[file$POS < 578,], aes(x = as.numeric(POS), y = as.numeric(gnomad_af_hom))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300, 400, 500, 577)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "gnomAD AF") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(0.5, 2.6, 0.25, 0), "cm")) 

plotB <- ggarrange(plot1, plot2, ncol = 2, nrow = 1, widths = c(5.5, 6))


# HelixMTdb
plot1 <- ggplot(file[file$POS > 16023,], aes(x = as.numeric(POS), y = as.numeric(helix_af_hom))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(16024, 16100, 16200, 16300, 16400, 16500)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "HelixMTdb AF") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0, 0.25, 0.5), "cm")) +
  guides(color = FALSE) 

plot2 <- ggplot(file[file$POS < 578,], aes(x = as.numeric(POS), y = as.numeric(helix_af_hom))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300, 400, 500, 577)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "HelixMTdb AF") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(0.5, 2.6, 0.25, 0), "cm")) 

plotC <- ggarrange(plot1, plot2, ncol = 2, nrow = 1, widths = c(5.5, 6))


# MITOMAP
plot1 <- ggplot(file[file$POS > 16023,], aes(x = as.numeric(POS), y = as.numeric(mitomap_af))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(16024, 16100, 16200, 16300, 16400, 16500)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "MITOMAP AF") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0, 0.25, 0.5), "cm")) +
  guides(color = FALSE) 

plot2 <- ggplot(file[file$POS < 578,], aes(x = as.numeric(POS), y = as.numeric(mitomap_af))) +
  geom_line() + 
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 100, 200, 300, 400, 500, 577)) + 
  scale_y_sqrt(limits = c(0,1), breaks = c(0.01, 0.1, 0.25, 0.5, 1)) + 
  labs(x = "Position", y = "MITOMAP AF") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(0.5, 2.6, 0.25, 0), "cm")) 

plotD <- ggarrange(plot1, plot2, ncol = 2, nrow = 1, widths = c(5.5, 6))


write.table(file[(file$POS > 16023) | (file$POS < 578), c("POS", "REF", "ALT", "gnomad_af_hom", "helix_af_hom", "mitomap_af")], 
            col.names = c("pos", "ref", "alt", "gnomad_af_hom", "helix_af_hom", "mitomap_af"),
            file = 'final_figures_source_data/FigureED8b-d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# collate figure 
ggarrange(plotA, plotB, plotC, plotD, ncol = 1, nrow = 4, labels = c("a", "b", "c", "d"), heights = c(1.5, 1, 1, 1), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED8.jpeg", width = 180, height = 165, dpi = 600, units = c("mm"))

