library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(png)

# EXTENDED FIGURE 5

# Figure ED5a - plot bases in regional missense constraint data as stacked bar to show proportions

# read in file and subset to protein genes
file <- read.delim(file = '../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file <- file[grepl("MT-A|MT-C|MT-N", file$symbol),]
file$in_rc <- factor(ifelse(file$in_rc == "no" & !is.na(file$min_distance_to_rc) & as.numeric(file$min_distance_to_rc) <= 6, "proximal", as.character(file$in_rc)),
                     levels = c("yes", "proximal", "no"), 
                     labels = c("Within", "Proximal", "Outside"))

# hack for plotting
file1 <- file
file1$plot <- "Protein"

# to also include functional sites in same plot
file2 <- file
file2$plot <- "Functional\nsites"

for_plot <- rbind(file1, file2[grepl("proton", file2$other_prot_annotation) | grepl("site", file2$uniprot_annotation), ])
for_plot$plot <- factor(for_plot$plot, levels = c("Protein", "Functional\nsites"))

plotA <- ggplot(for_plot, aes(plot, fill = in_rc)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(y = "Proportion", x = "Bases encoding") + 
  paper_theme +
  theme(legend.position = "left") +
  scale_fill_manual(values = c("#ff4040","#926281", "#377eb8"), 
                    labels = c("Within", "Proximal", "Outside"), name = "Regional\nconstraint")

# remove new lines for writing
for_plot$plot <- gsub("\n", " ", for_plot$plot)

write.table(unique(for_plot[,c("POS", "symbol", "protein_position", "plot", "in_rc")]), col.names = c("pos", "symbol", "protein_position", "group", "in_rc"),
            file = 'final_figures_source_data/FigureED5a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# count bp in each category, notes for paper
as.data.frame(file %>% group_by(in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
nrow(file[unique(file$POS),])
as.data.frame(file[grepl("proton", file$other_prot_annotation) | grepl("site", file$uniprot_annotation), ] %>% group_by(in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
nrow(file[grepl("proton", file$other_prot_annotation) | grepl("site", file$uniprot_annotation), ])


# Figure ED5b - plot disease associated missense vs regional constraint data as stacked bar to show proportions

# use ClinVar pathogenic and likely pathogenic variants and MITOMAP confirmed disease variants as pathogenic, and ClinVar benign variants as benign
# and missense most severe consequence
file <- file[grepl("missense", file$consequence) & !grepl("stop_gain|start_lost|stop_lost|incomplete_terminal", file$consequence),]
file$group <- factor(ifelse(grepl("Cfrm", file$mitomap_status) | grepl("athogenic", file$clinvar_interp), "Pathogenic",
                            ifelse(grepl("Benign", file$clinvar_interp), "Benign", "neither")),
                     levels = c("Benign", "Pathogenic", "neither"))

plotB <- ggplot(file[file$group != "neither", ], aes(group, fill = in_rc)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(x = '\nMissense variants', y = "Proportion") + 
  paper_theme +
  scale_fill_manual(values = c("#ff4040","#926281", "#377eb8"), labels = c("Within", "Proximal", "Outside"), 
                    guide = FALSE) 

write.table(file[file$group != "neither", c("POS", "REF", "ALT", "symbol", "consequence", "protein_position", "group", "in_rc")], 
            col.names = c("pos", "ref", "alt", "symbol", "consequence", "protein_position", "group", "in_rc"),
            file = 'final_figures_source_data/FigureED5b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)

# count in each category
as.data.frame(file[file$group != "neither", ] %>% group_by(group, in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))


# Figure ED5c - plot regional constraint rRNA data as stacked bar to show proportions

file <- read.delim(file='../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file <- file[grepl("MT-R", file$symbol),]
file$in_rc <- factor(ifelse(file$in_rc == "no" & !is.na(file$min_distance_to_rc) & (as.numeric(file$min_distance_to_rc) <= 6), "proximal", as.character(file$in_rc)),
                     levels = c("yes", "proximal", "no"), labels = c("Within", "Proximal", "Outside"))
file$RNA_modified <- ifelse(file$RNA_modified == "", "No", "Yes")

# hack for plotting
file1 <- file
file1$plot <- "rRNA"

# to also include modified and bridge bases in same plot
file2 <- file
file2$plot <- "Modified &\nbridge bases"

for_plot <- rbind(file1, file2[file2$RNA_modified == "Yes" | file$rRNA_bridge_base == "Yes",])
for_plot$plot <- factor(for_plot$plot, levels = c("rRNA", "Modified &\nbridge bases"))

plotC <- ggplot(for_plot, aes(plot, fill = in_rc)) + 
  geom_bar(position = "fill", colour = "black") +
  labs(y = "Proportion", x = "Bases encoding") + 
  paper_theme +
  theme(legend.position = "left") +
  scale_fill_manual(values = c("#ff4040","#926281", "#377eb8"), 
                    labels = c("Within", "Proximal", "Outside"), name = "Regional\nconstraint", guide = FALSE)

# remove new lines for writing
for_plot$plot <- gsub("\n", " ", for_plot$plot)

write.table(unique(for_plot[,c("POS", "symbol", "plot", "in_rc")]), col.names = c("pos", "symbol", "group", "in_rc"),
            file = 'final_figures_source_data/FigureED5c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# count in category, notes for paper
as.data.frame(file %>% group_by(in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
nrow(file[unique(file$POS),])
as.data.frame(file[file$RNA_modified == "Yes" | file$rRNA_bridge_base == "Yes",] %>% group_by(in_rc) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
nrow(file[file$RNA_modified == "Yes" | file$rRNA_bridge_base == "Yes",])


# Figure ED5d - chimeraX of MT-CYB and regional constraint

cyb <- readPNG("extended_data_figures/FigureED5d_CYB.png")  # 3D
plotD <- ggplot() + 
  background_image(cyb) +
  theme(plot.margin = unit(c(0, 0.25, 0, 0.25), "cm"))


# Figure ED5e - chimeraX of MT-ND6 and regional constraint

nd6 <- readPNG("extended_data_figures/FigureED5e_ND6.png")  # 3D
plotE <- ggplot() + 
  background_image(nd6) +
  theme(plot.margin = unit(c(0.5, 0.25, 0.5, 0.25), "cm"))


# Figure ED5f - chimera figure showing area of regional constraint in tertiary structure

chimera <- readPNG("extended_data_figures/FigureED5f_RNR1.png")  # 3D
plotF <- ggplot() + 
  background_image(chimera) +
  theme(plot.margin = unit(c(0.25, 0.5, 0.25, 0.5), "cm"))


# Figure ED5g - regional constraint within the secondary structure of MT-RNR1 

svg <- readPNG("extended_data_figures/FigureED5g_RNR1.png") 
plotG <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))


# Figure ED5h - regional constraint within the secondary structure of MT-RNR2 

svg <- readPNG("extended_data_figures/FigureED5h_RNR2.png")  # 3D
plotH <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"))



# compile panel
ggarrange(
  ggarrange(plotA, plotB, plotC, ncol = 3, nrow = 1, labels = c("a", "b", "c"), widths = c(1, 0.675, 0.675), font.label = list(size = 10)),
  ggarrange(plotD, plotE, plotF, ncol = 3, nrow = 1, labels = c("d", "e", "f"), widths = c(1, 0.5, 0.75), font.label = list(size = 10)),
  ggarrange(plotG, plotH, ncol = 2, nrow = 1, labels = c("g", "h"), widths = c(0.75, 1), font.label = list(size = 10)),
  nrow = 3, ncol = 1, heights = c(0.25, 0.3, 0.35))

ggsave("extended_data_figures/FigureED5.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))
