library(forcats)
library(ggplot2)
library(ggpubr)
library(spgs)
library(stringr)

# FIGURE 1

# Figure 1a - plot the mutational signature estimated by the model, across the reference mtDNA (excluding OriB-OriH)

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# convert to pyrimidine for plotting
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
file$pyr_mut <- ifelse(file$mut == "G>T","C>A",
                     ifelse(file$mut == "G>A", "C>T",
                            ifelse(file$mut == "G>C", "C>G",
                                   ifelse(file$mut == "A>T", "T>A",
                                          ifelse(file$mut == "A>C", "T>G",
                                                 ifelse(file$mut == "A>G", "T>C", file$mut))))))
file$pyr_tri <- ifelse(file$REF == "G" | file$REF == "A", reverseComplement(file$trinucleotide, case = "upper"), as.character(file$trinucleotide))
file$strand <- factor(ifelse(file$REF == "G" | file$REF == "A", "Heavy", "Light"), levels = c("Light", "Heavy"), labels = c("Reference / Light", "Reverse complement / Heavy"))

# exclude OriB-OriH m.191-16197
for_plot <- unique(file[file$POS > 191 & file$POS < 16197, c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")]) 
  
plotA <- ggplot(data = for_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Pastel1", name = "Strand") +
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide", y = "Likelihood") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0)) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25)

write.table(for_plot, file = 'final_figures_source_data/Figure1a.tsv', 
            col.names = c("likelihood", "trinucleotide", "pyrimidine_mutation", "pyrimidine_trinucleotide", "strand"), row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 1b - plot the observed:expected ratio and 90% confidence interval for functional classes of variation in mtDNA

file <- read.delim(file='../output_files/oe/consequences_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file <- file[!grepl("lost", file$consequence), ]

file$group <- factor(ifelse(grepl("RNA", file$consequence), "RNA", ifelse(file$consequence == "intergenic", "Other", "Protein")),
                     levels = c("Protein", "RNA", "Other"))
file$consequence <- factor(file$consequence, levels = c("synonymous", "missense", "stop_gain", "tRNA", "rRNA", "intergenic"), 
                           labels = c("Synonymous", "Missense", "Stop gain", "tRNA", "rRNA", "Non-coding"))
# colors used within the manuscript for each functional class
mycolors = c('#4daf4a', '#377eb8', '#b22222', '#ff7f00', '#984ea3', '#ffcc00')

plotB <- ggplot(file, aes(y = fct_rev(consequence), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.75, colour = "#777575") +
  geom_point(aes(colour = consequence), shape = 18, size = 6, show.legend = FALSE) +
  labs(x = "observed:expected ratio") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size = 7)) +
  facet_grid(rows = vars(group), scales = "free", space = 'free') + 
  scale_color_manual(values = mycolors) 

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[,c("consequence", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")],
            file = 'final_figures_source_data/Figure1b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 1c - plot the observed:expected ratio and 90% confidence interval for disease-associated variation

file <- read.delim(file = '../output_files/oe/disease_variants_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$group <- str_split(file$classification, "\\-", simplify = T)[,2]
file$classification <- factor(str_split(file$classification, "\\-", simplify = T)[,1], 
                              levels=c("Reported", "Cfrm", "Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic"),
                              labels=c("Reported", "Confirmed", "Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic"))

plotC <- ggplot(file, aes(y = fct_rev(classification), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.85, colour = "#777575") +
  geom_point(aes(colour = group), shape = 18, size = 6, show.legend = FALSE) +
  facet_grid(rows = vars(group), scales = "free", space = 'free') + 
  labs(x = "observed:expected ratio") + 
  paper_theme +
  xlim(0, 1) +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size = 7)) +
  scale_color_manual(values = c("dark blue", "dark green"))

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[,c("classification", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")], 
            file = 'final_figures_source_data/Figure1c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 1d-f - plot the observed vs expected synonymous, missense, stop gain variation in protein genes

gene_oe <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header=TRUE, sep = "\t")

# synonymous
plotD <- ggplot(data = gene_oe[gene_oe$consequence == "synonymous",], aes(y = observed, x = expected)) +
  geom_point(size = 1.5, color = "#4daf4a") + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(limits = c(0, 800), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 800), expand = c(0, 0)) + 
  labs(x = "Expected", y = "Observed") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0.5, 0.25, 0.5), "cm")) +
  geom_text(aes(label = "Synonymous"), x = 220, y = 750, colour = "#4daf4a", size = 3.5) +
  stat_cor(method = "pearson", label.x = 25, label.y = 650, digits = 3, aes(label = ..r.label..), size = 3)

# missense
plotE <- ggplot(data = gene_oe[gene_oe$consequence == "missense" ,], aes(y = observed, x = expected)) +
  geom_point(size = 1.5, color = "#377eb8") +
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(limits = c(0, 1500), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0)) + 
  labs(x = "Expected", y = "Observed") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0.5, 0.25, 0.25), "cm")) +
  geom_text(aes(label = "Missense"), x = 320, y = 1400, colour = "#377eb8", size = 3.5)

# stop gain 
plotF <- ggplot(data = gene_oe[gene_oe$consequence == "stop_gain",], aes(y = observed, x = expected)) +
  geom_point(size = 1.5, color = "#b22222") +
  geom_abline(intercept = 0, slope = 1) + 
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  labs(x = "Expected", y = "Observed") + 
  paper_theme +
  theme(plot.margin = unit(c(0.5, 0.5, 0.25, 0.25), "cm")) +
  geom_text(aes(label = "Stop gain"), x = 21, y = 95, colour = "#b22222", size = 3.5)

# to write table
gene_oe$expected <- round(gene_oe$expected, 4)

write.table(gene_oe[grepl("synonymous|missense|stop_gain", gene_oe$consequence), c("locus", "description", "consequence", "observed", "expected")], 
            file = 'final_figures_source_data/Figure1d-f.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# compile panel
ggarrange(plotA, 
          ggarrange(plotB, plotC, nrow = 1, ncol = 2, widths = c(0.9, 1), labels = c("b", "c"), font.label = list(size = 10)), 
          ggarrange(plotD, plotE, plotF, ncol = 3, nrow = 1, labels = c("d", "e", "f"), font.label = list(size = 10)),
          nrow = 3, ncol = 1, heights = c(0.65, 0.75, 0.7), labels = c("a", "", ""), font.label = list(size = 10))

ggsave("figures/Figure1.jpeg", width = 180, height = 160, dpi = 600, units = c("mm"))


