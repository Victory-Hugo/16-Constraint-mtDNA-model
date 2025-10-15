library(ggplot2)
library(tidyverse)

# SUPPLEMENTARY FIGURE 1
# plot a comparison of the likelihood of transitions across different de novo categories and maximum de novo counts per sample, using output of compare_denovo.py

file <- read.delim(file = '../output_files/denovo/Ts_likelihood_by_category.txt', header = TRUE, sep= "\t")

# format for plot - since the maximum per sample for germline is 3, remove those above this
file <- file[(file$category == "germline" & file$threshold < 4) | file$category != "germline", ]
# to collapse for plotting
file$threshold[file$threshold == 40] <- ">5"  
file$threshold <- factor(file$threshold, levels = c(1, 2, 3, 4, 5, ">5"))

file_plot <- file %>% gather(mutation, likelihood, likelihood_C.T:likelihood_A.G)
file_plot$category <- factor(file_plot$category, levels = c("germline", "germline + somatic tissue", "germline + somatic tissue + somatic cancer"),
                             labels = c("Germline", "Germline + somatic tissue", "Germline + somatic tissue + somatic cancer"))

ggplot(data = file_plot, aes(y = likelihood, x = threshold, color = mutation, group = mutation)) +
  geom_point(size = 2) + 
  geom_line() + 
  facet_grid(.~category, scales = "free", labeller = label_wrap_gen(width = 17)) + 
  labs(x = "Maximum number of mutations per sample", y = "Likelihood score") +
  scale_y_continuous(limits = c(0, 10.5)) + 
  paper_theme + 
  scale_color_discrete(labels = c("A>G", "C>T", "G>A", "T>C"), name = "Mutation type")
  
ggsave("supplementary_figures/FigureS1.jpeg", width = 180, height = 60, dpi = 600, units = c("mm"))


# compute n
denovo_raw <- read.delim(file = '../output_files/denovo/all_denovo.txt', header = TRUE, sep= "\t")
denovo <- unique(denovo_raw[,c("sample", "sample_denovo_count")])
denovo$type <- ifelse(grepl("germline", denovo$sample), "germline",
                      ifelse(grepl("somatic_tissue", denovo$sample), "somatic-tissue",
                             ifelse(grepl("somatic_cancer", denovo$sample), "somatic-cancer", "NA")))
table(denovo$type)
