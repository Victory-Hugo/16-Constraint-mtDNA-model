library(ggplot2)
library(ggpubr)

# SUPPLEMENTARY FIGURE 5

# Figure S5a - density plot showing expected distribution for various k-mer lengths, in protein genes

for (length in c(10, 15, 20, 25, 30)){
  file <- read.delim(file = sprintf('../output_files/local_constraint/%s_kmers_local_functional_constraint.txt', length), header = TRUE, sep = "\t")
  file$length <- as.factor(length)
  print(length)
  if(length == 10){
    kmers <- file
  } else {
    kmers <- rbind(kmers, file)
  }
}

plotA <- ggplot(kmers[!is.na(kmers$pctRank_OEUF_protein) & is.na(kmers$pctRank_OEUF_RNA) & is.na(kmers$pctRank_OEUF_noncoding),], aes(x = expected, fill =  length)) + 
  geom_density(alpha = 0.7) +
  labs(x = 'Expected amino acid substitutions in protein genes', y = 'Density', fill = 'k-mer length (bp)') + 
  xlim(0, 60) +
  ylim(0, 0.3) +
  paper_theme +
  geom_vline(xintercept = 10)


# Figure S5b - density plot showing expected distribution for various k-mer lengths, in RNA and non-coding loci

plotB <- ggplot(kmers[is.na(kmers$pctRank_OEUF_protein),], aes(x = expected, fill =  length)) + 
  geom_density(alpha = 0.7) +
  labs(x = 'Expected base substitutions in RNA and non-coding loci', y = 'Density', fill = 'k-mer length (bp)') + 
  xlim(0, 60) +
  ylim(0, 0.3) +
  paper_theme +
  geom_vline(xintercept = 10)


# compile figure panel
ggarrange(plotA, plotB, nrow = 2, labels = c("a", "b"), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS5.jpeg", width = 180, height = 90, dpi = 600, units = c("mm"))


# compute n

# in protein genes
kmers[!is.na(kmers$pctRank_OEUF_protein) & is.na(kmers$pctRank_OEUF_RNA) & is.na(kmers$pctRank_OEUF_noncoding),] %>% group_by(length) %>% summarize(n = n())
# in RNA genes and non-coding loci
kmers[is.na(kmers$pctRank_OEUF_protein),] %>% group_by(length) %>% summarize(n = n())
