library(ggplot2)
library(ggpubr)

# SUPPLEMENTARY FIGURE 3 

# confidence interval characteristics

file <- read.delim(file = '../output_files/other/CI_analyses.txt', header = TRUE, sep = "\t")
file$proportion <- round(file$exp/file$total, 2)

# define lines for plotting
data_vline <- data.frame(ratio = c(0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 1.5, 1.5, 1.5),
                         proportion = c(0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 0.1, 0.5, 1.0),
                         vline = c(0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 1.5, 1.5, NA))

# plot across ratios and proportions 
plotA <- ggplot(data = file[file$vp == 2 & file$proportion %in% c(1.0, 0.50, 0.10) & file$ratio %in% c(1.5, 1.0, 0.50, 0.10) & file$exp < 500,], aes(x = ratio, y = exp)) +
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), color = "grey", size = 0.5) +
  geom_point(aes(lower_CI), color = "dark red", size = 0.5) +
  geom_point(aes(upper_CI), color = "dark green", size = 0.5) +
  facet_grid(proportion ~ ratio, labeller = label_both) + 
  labs(x = "observed:expected ratio", y = "Expected value") + 
  paper_theme + 
  theme(strip.text.y = element_text(size = 6)) +
  geom_vline(data = data_vline, aes(xintercept = vline), linetype = "dashed", colour = "black", size = 0.5) +
  geom_hline(aes(yintercept = 10), linetype = "solid", colour = "black", size = 0.1)

# showing effect of vp
plotB <- ggplot(data = file[file$proportion == 0.10 & file$ratio %in% c(1.5, 1.0, 0.5, 0.1) & file$exp < 500,], aes(x = ratio, y = exp, color = as.character(vp))) +
  geom_point(aes(lower_CI), size = 0.5) +
  geom_point(aes(upper_CI), size = 0.5) +
  scale_color_manual(values = rev(c('#3690c0', '#034e7b', '#a6bddb')), name = "Parameter range") +
  facet_grid(~ratio, labeller = label_both) + 
  labs(x = "observed:expected ratio", y = "Expected") + 
  paper_theme + 
  theme(strip.text.y = element_text(size = 6)) +
  theme(legend.position = 'bottom') +
  #scale_y_continuous(breaks = c(0, 50, 100)) +
  geom_hline(aes(yintercept = 10), linetype = "solid", colour = "black", size = 0.1)

# collate figure 
ggarrange(plotA, plotB, ncol = 1, nrow = 2, labels = c("a", "b"), heights = c(0.7, 0.4), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS3.jpeg", width = 180, height = 150, dpi = 600, units = c("mm"))

