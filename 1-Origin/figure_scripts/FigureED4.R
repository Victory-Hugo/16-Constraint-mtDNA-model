library(data.table) 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)

# EXTENDED DATA FIGURE 4

# Figure ED4a - generate linear plots of regional missense constraint in the protein genes

dir.create("extended_data_figures/rdata/")

file <- read.delim(file = '../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file <- file[grepl("MT-A|MT-C|MT-N", file$symbol), ]
file <- data.frame(file %>% mutate(protein_position = strsplit(as.character(protein_position), ","), symbol = strsplit(as.character(symbol), ",")) %>% unnest(protein_position, symbol))
file$identifier <- paste(file$symbol, file$protein_position, sep = ":")
plot_file <- file[!duplicated(file$identifier), ]

# by gene
for(gene in list("MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6")){
  print(gene)
  plot <- ggplot(plot_file[plot_file$symbol == gene,], aes(x = as.numeric(as.character(protein_position)), y = 1)) +
    geom_bar(aes(fill = as.factor(in_rc), colour = as.factor(in_rc)), stat = "identity") +
    scale_fill_manual(values = c("#377eb8", "#ff4040"), guide = FALSE) +
    scale_colour_manual(values = c("#377eb8", "#ff4040"), guide = FALSE) + 
    labs(title = gene, x = "Residue number") + 
    scale_x_continuous(breaks = c(1, max(as.numeric(as.character(file[file$symbol == gene, c("protein_position")]))) - 1)) +
    paper_theme + 
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.05),
          axis.title.x = element_text(vjust = 5),
          plot.margin = unit(c(0.25, 0.25, 0.05, 0.25), "cm")) 
  save(plot, file = sprintf("extended_data_figures/rdata/%s.Rdata", gene))
}

load("extended_data_figures/rdata/MT-ATP6.Rdata")
plot0 <- plot
load("extended_data_figures/rdata/MT-ATP8.Rdata")
plot1 <- plot
load("extended_data_figures/rdata/MT-CO1.Rdata")
plot2 <- plot
load("extended_data_figures/rdata/MT-CO2.Rdata")
plot3 <- plot
load("extended_data_figures/rdata/MT-CO3.Rdata")
plot4 <- plot
load("extended_data_figures/rdata/MT-CYB.Rdata")
plot5 <- plot
load("extended_data_figures/rdata/MT-ND1.Rdata")
plot6 <- plot
load("extended_data_figures/rdata/MT-ND2.Rdata")
plot7 <- plot
load("extended_data_figures/rdata/MT-ND3.Rdata")
plot8 <- plot
load("extended_data_figures/rdata/MT-ND4.Rdata")
plot9 <- plot
load("extended_data_figures/rdata/MT-ND4L.Rdata")
plot10 <- plot
load("extended_data_figures/rdata/MT-ND5.Rdata")
plot11 <- plot
load("extended_data_figures/rdata/MT-ND6.Rdata")
plot12 <- plot

write.table(plot_file[,c("symbol", "protein_position", "in_rc")], 
            file = 'final_figures_source_data/FigureED4a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED4b - generate linear plots of regional constraint in the rRNA genes

file <- read.delim(file = '../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file <- file[grepl("MT-R", file$symbol), ]
file$identifier <- paste(file$symbol, file$POS, sep = ":")
plot_file <- file[!duplicated(file$identifier), ]

# by gene
for(gene in list("MT-RNR1", "MT-RNR2")){
  print(gene)
  plot <- ggplot(plot_file[plot_file$symbol == gene,], aes(x = as.numeric(as.character(POS)), y = 1)) +
    geom_bar(aes(fill = as.factor(in_rc), colour = as.factor(in_rc)), stat = "identity") +
    scale_fill_manual(values = c("#984ea3", "#ff4040"), guide = FALSE) +
    scale_colour_manual(values = c("#984ea3", "#ff4040"), guide = FALSE) + 
    labs(title = gene, x = "mtDNA position") + 
    scale_x_continuous(breaks = c(min(as.numeric(as.character(file[file$symbol == gene, c("POS")]))), 
                                  max(as.numeric(as.character(file[file$symbol == gene, c("POS")]))))) +
    paper_theme + 
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.05),
          axis.title.x = element_text(vjust = 5),
          plot.margin = unit(c(0.25, 0.25, 0.05, 0.25), "cm")) 
  save(plot, file = sprintf("extended_data_figures/rdata/%s.Rdata", gene))
}

load("extended_data_figures/rdata/MT-RNR1.Rdata")
plot13 <- plot
load("extended_data_figures/rdata/MT-RNR2.Rdata")
plot14 <- plot

write.table(plot_file[,c("symbol", "POS", "in_rc")], col.names = c("symbol", "pos", "in_rc"),
            file = 'final_figures_source_data/FigureED4b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# plot as figure panel
ggarrange(
  ggarrange(plot0, plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, nrow = 7, ncol = 2), 
  ggarrange(plot13, plot14, nrow = 1, ncol = 2),
  nrow = 2, ncol = 1, labels = c("a", "b"), heights = c(7, 1), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED4.jpeg", width = 180, height = 170, units = c("mm"))



