library(colorspace)
library(data.table)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(stringr)

# EXTENDED DATA FIGURE 9

# Figure ED9a - stacked bar plot of disease-associated variation by MLC score quartile

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")

# divide into four bins
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))

scores$group <- factor(ifelse(grepl("Cfrm", scores$mitomap_status) | (grepl("athogenic", scores$clinvar_interp)), "Pathogenic",
                              ifelse(grepl("Benign", scores$clinvar_interp), "Benign", "neither")) , 
                       levels = c("Pathogenic", "Benign", "neither"))

# stacked bar plot pathogenic vs benign by score quartile
plotA <- ggplot(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),], 
                aes(group, fill = fct_rev(rank_bin))) + 
  geom_bar(position = "fill", colour = "black") +
  labs(y = "Proportion", fill = "MLC score") + 
  paper_theme +
  theme(axis.title.x = element_blank(), 
        legend.key.height = unit(0.45, 'cm'), 
        legend.position = "left") + 
  scale_fill_manual(values = c("#ff4124", "#ffbfaa", "#cfb1ff", "#542eff")) 

write.table(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),
                   c("POS", "REF", "ALT", "rank_bin", "group")], 
            col.names = c("pos", "ref", "alt", "MLC_bin", "group"),
            file = 'final_figures_source_data/FigureED9a.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# summarize for reference
as.data.frame(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),] 
              %>% group_by(group, rank_bin) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
as.data.frame(scores[scores$group == "Pathogenic" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence)
                     & !grepl("nr|na", scores$mitomap_plasmy) & scores$mitomap_plasmy != "",] %>% group_by(mitomap_plasmy, rank_bin) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))


# Figure ED9b - density plot by pathogenic vs benign

plotB <- ggplot(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),], 
                aes(x = as.numeric(MLC_var_score), fill =  group)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Group') + 
  paper_theme +
  theme(legend.text = element_text(margin = margin(r = 35, unit = "pt")),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, 1.5)))

write.table(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),
                   c("POS", "REF", "ALT", "MLC_var_score", "group")], 
            col.names = c("pos", "ref", "alt", "MLC_var_score", "group"),
            file = 'final_figures_source_data/FigureED9b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9c - collapse disease plasmy into two groups - plot all pathogenic

scores$plasmy <- ifelse(grepl("\\+\\/", scores$mitomap_plasmy), "At homoplasmy", "Only at heteroplasmy")
plotC <- ggplot(scores[scores$group == "Pathogenic" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence)
                       & !grepl("nr|na", scores$mitomap_plasmy) & scores$mitomap_plasmy != "",], 
                aes(x = as.numeric(MLC_var_score), fill =  plasmy)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Pathogenic with\nMITOMAP plasmy status') + 
  paper_theme +
  xlim(0, 1) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, 1.5)))

write.table(scores[scores$group == "Pathogenic" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence)
                   & !grepl("nr|na", scores$mitomap_plasmy) & scores$mitomap_plasmy != "",
                   c("POS", "REF", "ALT", "MLC_var_score", "group", "plasmy")], 
            col.names = c("pos", "ref", "alt", "MLC_var_score", "group", "plasmy"),
            file = 'final_figures_source_data/FigureED9c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# write table for curation, for next plots
scores$var <- paste("m.", scores$POS, scores$REF, ">", scores$ALT, sep = "")
#write.table(scores[grepl("Cfrm", scores$mitomap_status), c("var", "symbol", "consequence", "mitomap_status", "mitomap_plasmy", "MLC_var_score")], 
#            file = 'figures/cnfm_by_status.txt', sep = "\t", row.names = FALSE, col.names = TRUE)


# Figure ED9d - collapse disease plasmy into two groups - curated confirmed (above file modified)

curated <- read.delim(file = '../required_files/other_annotations/cnfm_by_status_curated.txt', header = TRUE, sep = "\t")
curated$plasmy <- ifelse(grepl("\\+\\/", curated$curated_status), "At homoplasmy", "Only at heteroplasmy")

plotD <- ggplot(curated[grepl("missense|transcript", curated$consequence) & !grepl("gain|lost|terminal", curated$consequence),], aes(x = as.numeric(MLC_var_score), fill = plasmy)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Confirmed pathogenic\nwith curated plasmy status') + 
  paper_theme +
  xlim(0, 1) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "right",
        legend.margin = margin(c(0, 0, 0, 1.5)))

write.table(curated[grepl("missense|transcript", curated$consequence) & !grepl("gain|lost|terminal", curated$consequence),], 
            file = 'final_figures_source_data/FigureED9d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9e - boxplot to show the MLC score distribution for indels in population databases

# read in gnomAD, filter to PASS only
gnomad <- read.delim(file = '../required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv', header = TRUE, sep = "\t")
gnomad <- gnomad[gnomad$filters == "PASS",]
gnomad$pos <- gnomad$position
gnomad$type <- ifelse((nchar(as.character(gnomad$ref)) == 1 & nchar(as.character(gnomad$alt)) == 1), "SNV", "indel")

# read in HelixMTdb, extract ref, pos, alt
helix <- read.delim(file = '../required_files/databases/HelixMTdb_20200327.tsv', header = TRUE, sep = "\t")
helix <- helix[str_count(helix$alleles, ",") == 1, ]  # remove multiallelic sites
helix$ref <-  sub(",.*", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$alt <-  sub(".*,", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$pos <- sub("chrM:", "", helix$locus)
helix$type <- ifelse((nchar(as.character(helix$ref)) == 1 & nchar(as.character(helix$alt)) == 1), "SNV", "indel")

# read in MITOMAP, filter to variants observed in genbank
mitomap <- read.delim(file = '../required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt', header = TRUE, sep = "\t")
mitomap <- mitomap[mitomap$gbcnt > 0,]
mitomap$type <- ifelse((nchar(as.character(mitomap$ref)) == 1 & nchar(as.character(mitomap$alt)) == 1) & !grepl(":", mitomap$ref) & !grepl(":", mitomap$alt), "SNV", "indel")

# merge and annotate
gnomad$db <- "gnomAD"
helix$db <- "HelixMTdb"
mitomap$db <- "MITOMAP"
indels <- rbind(gnomad[gnomad$type == "indel", c("pos", "ref", "alt", "db")],
                helix[helix$type == "indel", c("pos", "ref", "alt", "db")],
                mitomap[mitomap$type == "indel", c("pos", "ref", "alt", "db")])

# annotate with position scores
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$pos <- scores$POS
indels <- merge(indels, scores[!duplicated(scores$POS), c("pos", "MLC_pos_score")], by = "pos", all.x = T)

as.data.frame(indels %>% group_by(db) %>% summarize(median = median(MLC_pos_score), n = n()))

plotE <- ggplot(data = indels, aes(x = db, y = MLC_pos_score, fill = db)) + 
  geom_boxplot(position = "dodge") +
  labs(y = 'MLC score', x = "Indels") + 
  paper_theme +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.35), "cm")) +
  ylim(0, 1) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c("#F9F6EE", "#F9F6EE", "#F9F6EE"), guide = FALSE)

write.table(indels[order(as.numeric(indels$pos)),], 
            file = 'final_figures_source_data/FigureED9e.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9f - plot MLC position scores vs PhyloP

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))

plotF <- ggplot(data = scores[!duplicated(scores$POS),], aes(x = rank_bin, y = phyloP_score, fill = rank_bin)) + 
  geom_boxplot() +
  labs(x = 'MLC score', y = 'PhyloP score') + 
  paper_theme +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.35), "cm")) +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  scale_fill_discrete_sequential(palette = "OrYel", guide = FALSE, order = c(1, 2, 3, 4)) 

# print count per bin
table(scores[!duplicated(scores$POS),c("rank_bin")])

write.table(scores[!duplicated(scores$POS), c("POS", "phyloP_score", "rank_bin")], 
            col.names = c("pos", "phyloP", "MLC_bin"),
            file = 'final_figures_source_data/FigureED9f.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9g - overlay MLC with chimpanzee data

# read in file
file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$chimp_match <- ifelse(as.character(file$REF) == as.character(file$chimp_ref), "yes", "no")

# annotate MLC scores with chimp data
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
file$var <- paste(file$POS, file$REF, file$ALT)
scores$var <- paste(scores$POS, scores$REF, scores$ALT)
for_plot <- merge(scores, file[,c("var", "chimp_match")], by = "var")

# plot
plotG <- ggplot(for_plot, aes(x = as.numeric(POS), y = as.numeric(MLC_pos_score), label = ' - chimp conserved')) +
  geom_line(aes(color = pctRank_mean_OEUF)) + 
  scale_color_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 16569)) + 
  ylim(0, 1) + 
  labs(color = "MLC score", x = "Position", y = "MLC score") + 
  paper_theme +
  theme(axis.title.x = element_text(size = 8, vjust = 5),
    legend.key.height = unit(0.25, 'cm'), 
    legend.margin = margin(-15, 0, 0, 0),
    plot.margin = unit(c(0.25, 0.75, 0.25, 0.25), "cm")) +
  geom_point(data = for_plot, mapping = aes(x = as.numeric(POS), y = 0), size = 1, shape = 124) + # using vertical dash
  geom_point(data = for_plot[for_plot$chimp_match == "no"
                             & grepl("intergenic|non_coding_transcript|missense|incomplete_terminal", for_plot$consequence)
                             ,],
             mapping = aes(x = as.numeric(POS), y = 0), size = 1, shape = 124, color = "white") + # using vertical dash
  geom_text(mapping = aes(x = 16570, y = 0), hjust = 0, size = 2.4) + 
  coord_cartesian(clip = 'off') 

write.table(for_plot[order(for_plot$POS), c("POS", "REF", "ALT", "consequence", "MLC_pos_score", "chimp_match")], 
            col.names = c("pos", "ref", "alt", "consequence", "MLC_score", "chimp_ref_match"),
            file = 'final_figures_source_data/FigureED9g.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9h - the local constraint score distribution of base/amino acids substitutions in gnomAD

# read in and merge in scores - note filter to PASS only in gnomAD
gnomad <- read.delim(file = '../required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv', header = TRUE, sep = "\t")
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
gnomad$var <- paste(gnomad$position, gnomad$ref, gnomad$alt)
scores$var <- paste(scores$POS, scores$REF, scores$ALT)
gnomad <- merge(gnomad[gnomad$filters == "PASS",], scores[, c("var", "MLC_var_score")], by = "var", all.x = T) 

# assign variant type and category
gnomad$type <- ifelse((nchar(as.character(gnomad$ref)) == 1 & nchar(as.character(gnomad$alt)) == 1), "SNV", "indel")
gnomad$category <- factor(ifelse(gnomad$AF_hom == 0, "Heteroplasmy\nonly", 
                                 ifelse(gnomad$AF_hom > 0 & gnomad$AF_hom < 1/50000, "Homoplasmy\nAF <0.002%", 
                                        "Homoplasmy\nAF ≥0.002%")),
                          levels = c("Homoplasmy\nAF ≥0.002%", "Homoplasmy\nAF <0.002%",  "Heteroplasmy\nonly"))
as.data.frame(gnomad %>% group_by(category) %>% summarize(n()))

plotH <- ggplot(data = gnomad[gnomad$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "gnomAD SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey", "white"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.35), "cm"))

# remove new lines for writing and subset
gnomad$category <- gsub("\n", " ", gnomad$category)
gnomad <- gnomad[gnomad$type == "SNV",]

write.table(gnomad[order(gnomad$position), c("position", "ref", "alt", "filters", "AF_hom", "AF_het", "MLC_var_score", "type", "category")],
            file = 'final_figures_source_data/FigureED9h.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9i - the local constraint score distribution of base/amino acids substitutions in HelixMTdb

# read in, format and merge in scores
helix <- read.delim(file = '../required_files/databases/HelixMTdb_20200327.tsv', header = TRUE, sep = "\t")
helix <- helix[str_count(helix$alleles, ",") == 1, ]  # remove multiallelic sites
helix$ref <-  sub(",.*", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$alt <-  sub(".*,", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$POS <- sub("chrM:", "", helix$locus)
helix$var <- paste(helix$POS, helix$ref, helix$alt)
helix <- merge(helix, scores[, c("var", "MLC_var_score")], by = "var", all.x = T) 

# assign variant type and category
helix$type <- ifelse((nchar(as.character(helix$ref)) == 1 & nchar(as.character(helix$alt)) == 1), "SNV", "indel")
helix$category <- factor(ifelse(helix$AF_hom == 0, "Heteroplasmy\nonly", 
                                ifelse(helix$AF_hom > 0 & helix$AF_hom < 1/50000, "Homoplasmy\nAF <0.002%", 
                                       "Homoplasmy\nAF ≥0.002%")),
                         levels = c("Homoplasmy\nAF ≥0.002%", "Homoplasmy\nAF <0.002%",  "Heteroplasmy\nonly"))
as.data.frame(helix %>% group_by(category) %>% summarize(n()))

plotI <- ggplot(data = helix[helix$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "HelixMTdb SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey", "white"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.35), "cm"))

# remove new lines for writing and subset
helix$category <- gsub("\n", " ", helix$category)
helix <- helix[helix$type == "SNV",]

write.table(helix[order(helix$locus), c("locus", "alleles", "AF_hom", "AF_het", "MLC_var_score", "type", "category")],
            file = 'final_figures_source_data/FigureED9i.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure ED9j - the local constraint score distribution of base/amino acids substitutions in MITOMAP

# read in, format and merge in scores
mitomap <- read.delim(file = '../required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt', header = TRUE, sep = "\t")
mitomap$var <- paste(mitomap$pos, mitomap$ref, mitomap$alt)
mitomap <- merge(mitomap[!duplicated(mitomap$var),], scores[, c("var", "MLC_var_score")], by = "var", all.x = T)

# assign variant type and category, also calculate allele frequency
mitomap$type <- ifelse((nchar(as.character(mitomap$ref)) == 1 & nchar(as.character(mitomap$alt)) == 1) & !grepl(":", mitomap$ref) & !grepl(":", mitomap$alt), "SNV", "indel")
mitomap$af <- mitomap$gbcnt/56910  # number of sequences in MITOMAP at time of download
mitomap$category <- factor(ifelse(mitomap$af < 1/50000, "AF <0.002%", "AF ≥0.002%"), levels = c("AF ≥0.002%", "AF <0.002%"))
as.data.frame(mitomap %>% group_by(category) %>% summarize(n()))

plotJ <- ggplot(data = mitomap[mitomap$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "MITOMAP SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.35), "cm"))

# subset for writing
mitomap <- mitomap[mitomap$type == "SNV",]

write.table(mitomap[order(mitomap$pos), c("pos", "ref", "alt", "af", "MLC_var_score", "type", "category")],
            file = 'final_figures_source_data/FigureED9j.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# compile panel
ggarrange(
  ggarrange(ggarrange(plotA, plotE, plotF, nrow = 3, labels = c("a", "e", "f"), font.label = list(size = 10), heights = c(0.3, 0.25, 0.35)), 
            ggarrange(plotB, plotC, plotD, nrow = 3, labels = c("b", "c", "d"), font.label = list(size = 10)), 
            ncol = 2, widths = c(0.8, 1), labels = c("a", ""), font.label = list(size = 10)),
  ggarrange(plotG, labels = c("g"), font.label = list(size = 10)), 
  ggarrange(plotH, plotI, plotJ, ncol = 3, nrow = 1, labels = c("h", "i", "j"), widths = c(1, 1, 0.75), font.label = list(size = 10)),
  nrow = 3, ncol = 1, heights = c(1.15, 0.45, 0.45))

ggsave("extended_data_figures/FigureED9.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))


