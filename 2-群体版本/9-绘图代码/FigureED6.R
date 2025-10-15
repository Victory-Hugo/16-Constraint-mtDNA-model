library(dplyr)
library(ggplot2)
library(ggpubr)
library(png)
library(scales)
library(stringr)

# EXTENDED FIGURE 6

# Figure ED6a,d - plot tRNA disease-associated variants as stacked bar to show proportions of (a) base types and (b) domains

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# use ClinVar pathogenic and likely pathogenic variants and MITOMAP confirmed disease variants as pathogenic, and ClinVar benign variants as benign
file$group <- factor(ifelse(grepl("Cfrm", file$mitomap_status) | grepl("athogenic", file$clinvar_interp), "Pathogenic",
                            ifelse(grepl("Benign", file$clinvar_interp), "Benign", "neither")),
                     levels = c("Pathogenic", "Benign", "neither"))
  
# base types
file$RNA_base_type <- factor(file$RNA_base_type, levels = c("WC", "non-WC", "loop-or-other", "loop-or-other,WC", "WC,loop-or-other", "WC,WC"), 
                           labels = c("WC pair", "Non-WC pair", "Loop or other", "loop-or-other,WC", "WC,loop-or-other", "WC,WC"))

plotA <- ggplot(file[grepl("MT-T", file$symbol) & file$group != "neither" & file$RNA_base_type != "" & !grepl(",", file$RNA_base_type),], aes(group, fill = RNA_base_type)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(x = 'tRNA variants', y = "Proportion", fill = "Base type") + 
  paper_theme 

write.table(file[grepl("MT-T", file$symbol) & file$group != "neither" & file$RNA_base_type != "" & !grepl(",", file$RNA_base_type), 
                 c("POS", "REF", "ALT", "symbol", "RNA_base_type", "group")], 
            col.names = c("pos", "ref", "alt", "symbol", "RNA_base_type", "group"),
            file = 'final_figures_source_data/FigureED6a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# count in each category
as.data.frame(file[grepl("MT-T", file$symbol) & file$group != "neither" & file$RNA_base_type != "" & !grepl(",", file$RNA_base_type),] 
              %>% group_by(group, RNA_base_type) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))

# domain
file$tRNA_domain <- factor(file$tRNA_domain, levels = c("", "Acceptorstem", "D-stem", "D-loop", "Anticodonstem", "Anticodonloop", "Variableregion", "T-stem", "T-loop"), 
                           labels = c("", "Acceptor stem", "D-stem", "D-loop", "Anticodon stem", "Anticodon loop", "Variable region", "T-stem", "T-loop"))

plotD <- ggplot(file[grepl("MT-T", file$symbol) & file$group != "neither" & file$tRNA_domain != "",] , aes(group, fill = tRNA_domain)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(x = '\ntRNA variants', y = "Proportion", fill = "Domain") + 
  paper_theme +
  theme(legend.position = 'left')

write.table(file[grepl("MT-T", file$symbol) & file$group != "neither" & file$tRNA_domain != "", 
                 c("POS", "REF", "ALT", "symbol", "tRNA_domain", "group")], 
            col.names = c("pos", "ref", "alt", "symbol", "tRNA_domain", "group"),
            file = 'final_figures_source_data/FigureED6d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED6b - plot the observed:expected ratio and 90% CI for modified bases

file <- read.delim(file = '../output_files/oe/RNA_base_types_obs_exp.txt', header = TRUE, sep = "\t")
file$group <- factor(ifelse(grepl("tRNA", file$RNA_base_type), "tRNA", "rRNA"), levels = c("tRNA", "rRNA"))
file$RNA_base_type <- factor(str_sub(file$RNA_base_type, start = 1, end = -6), 
                             levels = c("WC", "non-WC", "loop-or-other", "WC", "non-WC", "loop-or-other", "modified", "modified", "non-modified", "non-modified"), 
                             labels = c("WC pair", "Non-WC pair", "Loop or other", "WC pair", "Non-WC pair", "Loop or other", "Modified", "Modified", "Non-modified", "Non-modified"))
# use RNA colors
mycolors = c('#ff7f00', '#984ea3')

plotB <- ggplot(file[grepl("odified", file$RNA_base_type),], aes(y = fct_rev(RNA_base_type), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.55, colour = "#777575") +
  geom_point(aes(colour = group), shape = 18, size = 7) +
  labs(x = "Observed:expected ratio") + 
  paper_theme +
  theme(axis.title.y = element_blank()) +
  xlim(0, 1) +
  scale_color_manual(values = mycolors, name = "RNA type") +
  facet_grid(rows = vars(group), scales = "free", space = 'free') +
  guides(colour = FALSE)

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[grepl("odified", file$RNA_base_type), c("RNA_base_type", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")],
            file = 'final_figures_source_data/FigureED6b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED6c - svg to with same color gradient to show location of domains

# color palette used can be extracted via hue_pal()(8)
svg <- readPNG("extended_data_figures/FigureED6c_tRNA_domains.png")  # 3D
plotC <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


# compile panel
ggarrange(
  ggarrange(plotA, plotB, labels = c("a", "b"), nrow = 1, ncol = 2, font.label = list(size = 10)),
  ggarrange(plotC, plotD, labels = c("c", "d"), nrow = 1, ncol = 2, font.label = list(size = 10)),
  nrow = 2, ncol = 1, heights = c(0.7, 1))

ggsave("extended_data_figures/FigureED6.jpeg", width = 180, height = 150, dpi = 600, units = c("mm"))
