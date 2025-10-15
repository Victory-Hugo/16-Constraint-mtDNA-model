library(ggplot2)
library(ggpubr)
library(png)

# SUPPLEMENTARY FIGURE 7 

# Figure S7a - schematic overview of method used for MLC

fig <- readPNG("supplementary_figures/FigureS7a.png")  
plotA <- ggplot() + 
  background_image(fig) +
  theme(panel.background = element_rect(fill = 'white', color = 'white'))


# Figure S7b - plot MLC position scores vs mean OEUF values

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")

plotB <- ggplot(scores[!duplicated(scores$POS),], aes(x = MLC_pos_score, y = mean_OEUF)) +
  geom_point(size = 0.25, stroke = 0) +
  labs(x = 'MLC score', y = 'mean OEUF') + 
  paper_theme #+
  #theme(plot.margin = unit(c(0, 5, 0.25, 5), "cm"))


# compile panel
ggarrange(plotA, plotB, nrow = 1, ncol = 2, widths = c(0.55, 0.45), labels = c("a", "b"), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS7.jpeg", width = 180, height = 50, dpi = 600, units = c("mm"))