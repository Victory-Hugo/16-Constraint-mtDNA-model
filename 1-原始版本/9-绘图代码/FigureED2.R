library(forcats)
library(ggplot2)
library(ggpubr)
library(spgs)
library(stringr)
library(tidyverse)

# EXTENDED DATA FIGURE 2

# Figure ED2a - plot the mutational signature estimated by the model, across the OriB-OriH region

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# convert to pyridimine mutation
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
file$pyr_mut <- ifelse(file$mut == "G>T","C>A",
                     ifelse(file$mut == "G>A", "C>T",
                            ifelse(file$mut == "G>C", "C>G",
                                   ifelse(file$mut == "A>T", "T>A",
                                          ifelse(file$mut == "A>C", "T>G",
                                                 ifelse(file$mut == "A>G", "T>C", file$mut))))))
file$pyr_tri <- ifelse(file$REF == "G" | file$REF == "A", reverseComplement(file$trinucleotide, case = "upper"), as.character(file$trinucleotide))
file$strand <- factor(ifelse(file$REF == "G" | file$REF == "A", "Heavy", "Light"), levels = c("Light", "Heavy"), labels = c("Reference / Light", "Reverse complement / Heavy"))

# subset to ori
ori_plot <- unique(file[file$POS < 192 | file$POS > 16196, c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")])

plotA <- ggplot(data = ori_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Pastel1", name = "Strand") +
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide in OriB-OriH", y = "Likelihood") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
        legend.position = "top",
        legend.key.size = unit(3, "mm")) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25)

write.table(ori_plot, file = 'final_figures_source_data/FigureED2a.tsv', 
            col.names = c("ori_likelihood", "ori_trinucleotide", "ori_pyrimidine_mutation", "ori_pyrimidine_trinucleotide", "ori_strand"), 
            row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED2b - disease-associated variants by consequence 

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$consequence_two <- factor(ifelse(grepl("MT-T", file$symbol), 'tRNA', 
                                      ifelse(grepl("MT-R", file$symbol), 'rRNA', 
                                             ifelse(grepl(",|lost|retained|incomplete|intergenic", file$consequence), "Other", 
                                                    as.character(file$consequence)))),
                               levels = c("synonymous_variant", "missense_variant", "stop_gained", "tRNA", "rRNA", "Other"),
                               labels = c("Synonymous", "Missense", "Stop gain", "tRNA", "rRNA", "Other"))
clinvar_included <- c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")
plot_data_clinvar <- as.data.frame(file[file$clinvar_interp %in% clinvar_included,] %>% group_by(consequence_two) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
plot_data_mitomap <- as.data.frame(file[grepl("Cfrm|Reported", file$mitomap_status),] %>% group_by(consequence_two) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))

mycolors = c('#4daf4a', '#377eb8', '#b22222', '#ff7f00', '#984ea3', 'grey')

plotB1 <- ggplot(plot_data_clinvar, aes(x = "", y = n, fill = consequence_two)) +
  geom_bar(stat = "identity", color = "white", show.legend = FALSE) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycolors) +
  labs(title = "ClinVar") +
  geom_text(aes(x = 1.63, label = ifelse(consequence_two == "rRNA", "      3%", ifelse(consequence_two == "Synonymous", "5%   ", paste0(round(freq * 100, 0), "%")))), 
            position = position_stack(vjust = .5), size = 2) +
  theme_void() +
  theme(plot.margin = unit(c(0.7, 0.1, 0.7, 0), "cm"),
        plot.title = element_text(size = 8, hjust = 0.5))

plotB2 <- ggplot(plot_data_mitomap, aes(x = "", y = n, fill = consequence_two)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycolors) +
  labs(fill = "Consequences", title = "MITOMAP") +
  geom_text(aes(x = 1.63, label = paste0(round(freq * 100, 0), "%")), 
            position = position_stack(vjust = .5), size = 2) +
  theme_void() +
  theme(plot.margin = unit(c(0.0, 0.25, 0, 0), "cm"),
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(4, "mm"),
        plot.title = element_text(size = 8, hjust = 0.5))

plot_data_clinvar$group <- "clinvar"
plot_data_mitomap$group <- "mitomap"
write.table(rbind(plot_data_clinvar, plot_data_mitomap), col.names = c("consequence", "n", "frequency", "group"),
            file = 'final_figures_source_data/FigureED2b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED2c - plot the observed:expected ratio and 90% confidence interval for subsets of VUS

file <- read.delim(file = '../output_files/oe/disease_variants_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$group <- str_split(file$classification, "\\-", simplify = T)[,2]
file$classification <- factor(str_split(file$classification, "\\-", simplify = T)[,1], 
                              levels=c("Reported", "Cfrm", "Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic"),
                              labels=c("Reported", "Confirmed", "Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic"))

# combine with vus by criteria data
vus <- read.delim(file = '../output_files/oe/vus_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
vus$group <- str_split(vus$classification, "\\-", simplify = T)[,1]
file <- rbind(file[grepl("Reported|Uncertain", file$classification),], vus[grepl("PP|BP", vus$classification),])
file$classification <- factor(file$classification, 
                              levels=c("Reported", "MITOMAP-Reported-PP3-PM2s", "MITOMAP-Reported-BP4-BS1", "Uncertain significance", "ClinVar-Uncertain significance-PP3-PM2s"),
                              labels=c("Reported (all)", "Reported, with pathogenic criteria", "Reported, with benign criteria", "Uncertain significance (all)", "Uncertain significance,\nwith pathogenic criteria"))

plotC <- ggplot(file, aes(y = fct_rev(classification), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.75, colour = "#777575") +
  geom_point(aes(colour = group), shape = 18, size = 6, show.legend = FALSE) +
  facet_grid(rows = vars(group), scales = "free", space = 'free', switch = 'y') + 
  labs(x = "observed:expected ratio") + 
  paper_theme +
  scale_y_discrete(position = "right") +
  xlim(0, 1.05) +
  theme(plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"),
        axis.text.x  = element_text(size = 6), 
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_color_manual(values = c("dark blue", "dark green"))

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[, c("classification", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")],
            file = 'final_figures_source_data/FigureED2c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED2d - plot the observed:expected ratio and 90% confidence interval for different categories of in silico prediction

file <- read.delim(file = '../output_files/oe/insilicos_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$group <- factor(str_split(file$prediction, "\\-", simplify = T)[,2], 
                     levels = c("APOGEE", "MitoTip", "HmtVar")) 
file$group2 <- ifelse(file$group == "APOGEE", "Missense", "tRNA")
file$prediction <- factor(str_split(file$prediction, "\\-", simplify = T)[,1], 
                          levels=c("polymorphic", "likely_polymorphic", "likely_pathogenic", "pathogenic", "likely benign", "possibly benign", "possibly pathogenic", "likely pathogenic", "Neutral", "Pathogenic"),
                          labels=c("Polymorphic", "Likely polymorphic", "Likely pathogenic", "Pathogenic", "Likely benign", "Possibly benign", "Possibly pathogenic ", "Likely pathogenic ", "Neutral", " Pathogenic"))

plotD <- ggplot(file, aes(y = prediction, x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.9, colour = "#777575") +
  geom_point(aes(colour = group2), shape = 18, size = 6) +
  labs(x = "observed:expected ratio") + 
  paper_theme +
  facet_grid(rows = vars(group), scales = "free", space = 'free') +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size = 7)) +
  scale_color_manual(values = c('#377eb8', '#ff7f00'), name = "Algorithm for:") +
  guides(color = FALSE) 

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[, c("prediction", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")],
            file = 'final_figures_source_data/FigureED2d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED2e - plot the observed:expected ratio and 90% confidence interval for functional classes of variation in mtDNA in HelixMTdb

file <- read.delim(file = '../output_files/oe/replication_dataset/helix_consequences_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file <- file[!grepl("lost", file$consequence), ]

file$group <- factor(ifelse(grepl("RNA", file$consequence), "RNA", 
                            ifelse(file$consequence == "intergenic", "Noncoding", 
                                   "Protein")),
                     levels = c("Protein", "RNA", "Noncoding"),
                     labels = c("Protein", "RNA", "Other"))
file$consequence <- factor(file$consequence, levels = c("synonymous", "missense", "stop_gain", "tRNA", "rRNA", "intergenic"), 
                           labels = c("Synonymous", "Missense", "Stop gain", "tRNA", "rRNA", "Intergenic"))

mycolors = c('#4daf4a', '#377eb8', '#b22222', '#ff7f00', '#984ea3', '#ffcc00')

plotE <- ggplot(file, aes(y = fct_rev(consequence), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.65, colour = "#777575") +
  geom_point(aes(colour = consequence), shape = 18, size = 6, show.legend = FALSE) +
  labs(x = "observed:expected ratio in HelixMTdb") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(rows = vars(group), scales="free", space='free') + 
  scale_color_manual(values = mycolors)

# to write table
file$obs.exp <- round(file$obs.exp, 4)

write.table(file[,c("consequence", "variant_count", "obs.exp", "lower_CI", "upper_CI", "group")],
            file = 'final_figures_source_data/FigureED2e.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# compile panel
ggarrange(plotA, 
          ggarrange(plotB1, plotB2, plotC, nrow = 1, ncol = 3, widths = c(0.55, 1, 1.75), labels = c("b", "", "c"), font.label = list(size = 10)),
          ggarrange(plotD, plotE, nrow = 1, ncol = 2, widths = c(1, 0.9), labels = c("d", "e"), font.label = list(size = 10)),
          nrow = 3, ncol = 1, heights = c(0.7, 0.7 , 1), labels = c("", "", ""), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED2.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))  
