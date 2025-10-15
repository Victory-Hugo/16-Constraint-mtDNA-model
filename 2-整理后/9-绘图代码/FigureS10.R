library(ggplot2)
library(ggpubr)
library(png)


# SUPPLEMENTARY FIGURE 10 - analysis of bystander edits and off-target rates

# Figure S10a - sanger sequencing results with bystanders

png <- readPNG("supplementary_figures/FigureS10a_sanger.png") 
plotA <- ggplot() + 
  background_image(png) +
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0), "cm"))


# Figure S10b - off-target rates

file <- read.delim(file = '../base_editing/off_targets/outputs/off_targets.txt', header = TRUE, sep = "\t")

file$name <- factor(file$name, levels = c('3047-dead', '3047-edited', '3075-dead', '3075-edited'), 
                    labels = c("Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"))

plotB <- ggplot(data = file, aes(x = name, y = rate)) +
  geom_bar(stat = "identity") +
  geom_point(size = 1) + 
  paper_theme +
  theme(axis.title.x = element_text(size = 6), 
        axis.text.x  = element_text(size = 6, angle = 40, hjust=1), 
        axis.title.y = element_text(size = 6),
        axis.text.y  = element_text(size = 6),
        plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm")) + 
  scale_y_continuous(limits = c(0, 0.008), labels = scales::percent) +
  labs(x = 'DdCBE', y = 'Off-target rate (%)') 


# collate figure

ggarrange(plotA, plotB, labels = c("a", "b"), font.label = list(size = 10), nrow = 1, widths = c(0.8, 0.2))

ggsave("supplementary_figures/FigureS10.jpeg", width = 180, height = 50, dpi = 600, units = c("mm"))
