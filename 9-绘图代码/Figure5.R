library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(png)
library(tidyr)


# FIGURE 5

# Figure 5a - sanger sequencing results

png <- readPNG("figures/Figure5a_sanger.png") 
plotA <- ggplot() + 
  background_image(png) +
  theme(plot.margin = unit(c(0.05, 0.25, 0.05, 0.25), "cm"),
        panel.background = element_blank())


# Figure 5b - cell confluence measurements in base edited cells

file_long <- read.delim(file = '../base_editing/functional_data/incucyte_cell_growth.tsv', header = TRUE, sep = "\t")

for_plot <- as.data.frame(file_long[!is.na(file_long$value) & !grepl("Std", file_long$variable),] %>% group_by(Elapsed, variable) %>% 
                            summarize(mean = mean(value), sd = sd(value), n = n(), sem = sd(value)/sqrt(n())))

for_plot$media <- factor(ifelse(grepl("glu", for_plot$variable), "GLU", "GAL"), levels = c("GLU", "GAL"))
for_plot$group <- ifelse(grepl("dead", for_plot$variable), "dead", "edited")
for_plot$target <- factor(substr(for_plot$variable, 2, 5), levels = c("5147", "3047", "3075"))

# define colors and linetypes
for_plot$label <- factor(paste(for_plot$group, for_plot$target, sep = "-"), 
                         levels = c("dead-5147", "edited-5147", "dead-3047", "edited-3047", "dead-3075", "edited-3075"), 
                         labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"))

colors <- c("Dead-ND2-5147" = "#4CBB17", "Dead-RNR2-3075" = "#880085", "Dead-RNR2-3047" = "#6050dc", "ND2-5147" = "#4CBB17", "RNR2-3075"  = "#880085", "RNR2-3047" = "#6050dc")
alphas <- c("Dead-ND2-5147" = 1.0, "Dead-RNR2-3047" = 1.0, "Dead-RNR2-3075" = 1.0, "ND2-5147" = 0.35,  "RNR2-3047"  = 0.35, "RNR2-3075" = 0.35)

# exclude timepoints missing measurements - 12 for 3075, 12 for 5147 and 39 for 3047
for_plot$exclude <- paste(for_plot$Elapsed, for_plot$target)
for_plot <- for_plot[!grepl("39 3047", for_plot$exclude) & !grepl("12 3075", for_plot$exclude) & !grepl("12 5147", for_plot$exclude),]

# plot
plotB <- ggplot(data=for_plot, aes(x = Elapsed, y = mean, color = label)) +
  geom_line(aes(alpha = label)) +
  geom_point(size = 1, aes(alpha = label)) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = c(0, 0), limits = c(0, 110)) +
  facet_grid(media~target, scales = "free") + 
  labs(x = 'Time (hours)', y = 'Cell confluency %', color = 'Group') +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem, alpha = label), width = .5) + 
  paper_theme +
  theme(legend.key.height = unit(0.5, 'cm'),
        panel.spacing = unit(0.35, "cm", data = NULL),
        plot.margin = unit(c(0.25, 0.4, 0.25, 0.4), "cm")) + 
  scale_alpha_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047",  "Dead-RNR2-3075", "RNR2-3075"), 
                     values = alphas, name = "DdCBE") +
  scale_color_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047",  "Dead-RNR2-3075", "RNR2-3075"), 
                     values = colors, name = "DdCBE",
                     guide = guide_legend(override.aes = list(alphas = c(1.0, 0.35))))

# to write table
for_plot$mean <- round(for_plot$mean, 4)
for_plot$sem <- round(for_plot$sem, 4)

write.table(for_plot[,c("Elapsed", "mean", "sem", "media", "group", "target", "label")], 
            col.names = c("time_elapsed_hours", "confluency_mean", "sem", "media", "group", "target", "label"), 
            file = 'final_figures_source_data/Figure5b.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 5c - seahorse measurements in base edited cells

# load file and standarize minutes
file <- read.delim(file = '../base_editing/functional_data/seahorse_respiration.tsv', header = TRUE, sep = "\t")

# compute mean and sd across experimental replicates
for_plot <- as.data.frame(file %>% group_by(group, minute) %>% summarize(mean = mean(units), sd = sd(units), n = n(), sem = sd(units)/sqrt(n())))
for_plot$group <- factor(for_plot$group, 
                         levels = c("Dead-5147", "5147", "Dead-3047", "3047", "Dead-3075", "3075"), 
                         labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"))

# define colors and linetypes
colors <- c("Dead-ND2-5147" = "#4CBB17", "Dead-RNR2-3075" = "#880085", "Dead-RNR2-3047" = "#6050dc", "ND2-5147" = "#4CBB17", "RNR2-3075"  = "#880085", "RNR2-3047" = "#6050dc")
linetypes <- c("Dead-ND2-5147" = "solid", "Dead-RNR2-3047" = "solid", "Dead-RNR2-3075" = "solid", "ND2-5147" = "dashed",  "RNR2-3047"  = "dashed", "RNR2-3075" = "dashed")

# plot
plotC <- ggplot(data = for_plot, aes(x = minute, y = mean, color = group)) +
  geom_line(aes(linetype = group)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = 'Time (minutes)', y = 'OCR (pmol/min/confluence)') +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .75) + 
  paper_theme +
  theme(legend.key.height = unit(0.5, 'cm'),
        plot.margin = unit(c(0.25, 0.4, 0.25, 0.4), "cm")) + 
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  geom_vline(xintercept = c(18, 38, 59), linetype = "dashed", colour = "dark grey", size = 0.25) +
  scale_linetype_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047",  "Dead-RNR2-3075", "RNR2-3075"), 
                        values = linetypes, name = "DdCBE") +
  scale_color_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047",  "Dead-RNR2-3075", "RNR2-3075"), 
                     values = colors, name = "DdCBE",
                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  annotate("text", x = c(18, 38, 59), y = 2.7, label = c('Oligomycin', 'FCCP', 'Rotenone +\nAntimycin A'), size = 2.5) +
  coord_cartesian(ylim = c(0, 2.5), clip = 'off')

# to write table
for_plot$mean <- round(for_plot$mean, 4)
for_plot$sem <- round(for_plot$sem, 4)

write.table(for_plot[,c("group", "minute", "mean", "sem")], col.names = c("group", "minute", "OCR_mean", "sem"), 
            file = 'final_figures_source_data/Figure5c.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# Figure 5d - relative seahorse measurements in base edited cells

file <- read.delim(file = '../base_editing/functional_data/seahorse_relative_OCR.tsv', header = TRUE, sep = "\t")

# calculate the mean values for each assay across replicates, across all dead controls - this also captures variation in dead plasmid controls
controls <- as.data.frame(file[grepl("Dead", file$group),] %>% group_by(assay) %>% summarize(ctl_mean = mean(unit), sd = sd(unit), n = n(), sem = sd(unit)/sqrt(n())))
file <- merge(file, controls[,c("ctl_mean", "assay")], by = "assay")
file$relative_unit <- file$unit/file$ctl_mean

# taking mean across biological replicates, using values normalized to the mean of all dead controls, across experiments, for each assay
for_plot <- as.data.frame(file %>% group_by(group, assay) %>% summarize(rel_mean = mean(relative_unit), sd = sd(relative_unit), n = n(), sem = sd(relative_unit)/sqrt(n())))

# as check, this is the mean of all dead controls for each assay; i.e. 1.0
as.data.frame(file[grepl("Dead", file$group),] %>% group_by(assay) %>% summarize(mean = mean(relative_unit)))

# compute p-values, relative to all dead controls 
assays = c("ATP_production", "Basal_respiration", "Maximal_respiration", "Proton_leak")
targets = c("5147", "3047", "3075")

sink("../base_editing/functional_data/relative_OCR_pvalues.txt")
for(assay in assays){
  for(target in targets){
    for(dead in targets){
      print(paste(assay, ":", target, "vs Dead-", dead, sep = " "))
      print(t.test(file[file$group == target & file$assay == assay, c("relative_unit")],
                   file[file$group == paste("Dead-", dead, sep="") & file$assay == assay, c("relative_unit")],
                   alternative = c("two.sided"),
                   conf.level = 0.95)$p.value)
    }
  }
}
sink()


# define levels and labels for plotting
for_plot$assay <- factor(for_plot$assay, 
                         levels = c("Basal_respiration", "ATP_production", "Proton_leak", "Maximal_respiration"), 
                         labels = c("Basal respiration", "ATP production", "Proton leak", "Maximal respiration"))
for_plot$group <- factor(for_plot$group, 
                         levels = c("Dead-5147", "5147", "Dead-3047", "3047", "Dead-3075", "3075"), 
                         labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"))

# define levels and labels for plotting for the source file
file$assay <- factor(file$assay, 
                     levels = c("Basal_respiration", "ATP_production", "Proton_leak", "Maximal_respiration"), 
                     labels = c("Basal respiration", "ATP production", "Proton leak", "Maximal respiration"))

# define colors and alphas for plotting
colors <- c("Dead-ND2-5147" = "#4CBB17", "Dead-RNR2-3075" = "#880085", "Dead-RNR2-3047" = "#6050dc", "ND2-5147" = "#4CBB17", "RNR2-3075"  = "#880085", "RNR2-3047" = "#6050dc")
alphas <- c("Dead-ND2-5147" = 1.0, "Dead-RNR2-3047" = 1.0, "Dead-RNR2-3075" = 1.0, "ND2-5147" = 0.35,  "RNR2-3047"  = 0.35, "RNR2-3075" = 0.35)
file$group <- factor(file$group, 
                     levels = c("Dead-5147", "5147", "Dead-3047", "3047", "Dead-3075", "3075"), 
                     labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"))

# plot
plotD <- ggplot(data = for_plot, aes(x = assay, y = rel_mean, fill = group)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(), aes(alpha = group)) +
  geom_point(data = file, aes(x = assay, y = relative_unit, color = group), position = position_jitterdodge(dodge.width=0.9, jitter.width = 0.05), size = 0.25, color = "black") +
  labs(y = 'Relative OCR') +
  geom_errorbar(aes(ymin = rel_mean - sem, ymax = rel_mean + sem), width = .6, position = position_dodge(.9), size = 0.25, color = "black") + 
  paper_theme +
  theme(axis.title.x = element_blank(), 
        legend.key.height = unit(0.4, 'cm'), 
        plot.margin = unit(c(0.2, 0.4, 0.25, 0.4), "cm")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25)) + 
  scale_fill_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"), 
                    values = colors, name = "DdCBE") +
  scale_alpha_manual(labels = c("Dead-ND2-5147", "ND2-5147", "Dead-RNR2-3047", "RNR2-3047", "Dead-RNR2-3075", "RNR2-3075"),
                     values = alphas, name = "DdCBE") +
  geom_hline(yintercept = c(1), linetype = "dashed", colour = "dark grey", size = 0.25) +
  geom_signif(y_position = rep(1.3, 8), xmin = c(0.92, 1.22, 1.92, 2.22, 2.92, 3.22, 3.62, 3.92), xmax = c(1.07, 1.37, 2.07, 2.37, 3.07, 3.37, 3.77, 4.07), annotations = c("*", "**", "*", "**", "**", "*", "*", "*"), size = 0.25, tip_length = 0.01) +
  coord_cartesian(ylim = c(0, 1.4), clip = 'off')

# to write table
for_plot$rel_mean <- round(for_plot$rel_mean, 4)
for_plot$sem <- round(for_plot$sem, 4)

write.table(for_plot[,c("group", "assay", "rel_mean", "sem")], col.names = c("group", "assay", "relative_OCR", "sem"), 
            file = 'final_figures_source_data/Figure5d.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


# collate
ggarrange(plotA, plotB, plotC, plotD, labels = c("a", "b", "c", "d"), font.label = list(size = 10), nrow = 4, heights = c(0.45, 0.45, 0.45, 0.3))

ggsave("figures/Figure5.jpeg", width = 180, height = 165, dpi = 600, units = c("mm")) # 170 height is max


