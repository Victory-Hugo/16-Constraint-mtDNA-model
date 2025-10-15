library(circlize)
library(dplyr)
library(ggExtra)
library(ggplot2)
library(ggpubr)
library(png)
library(R.utils)
"%ni%" <- Negate("%in%")
library(scales)

# SUPPLEMENTARY FIGURE 9

# Figure S9a - circos plot showing MLC and UK Biobank data

# Adapted code from Nature Communications paper: https://github.com/ArkingLab/Heteroplasmy-and-Mortality/blob/main/circos_for_UKB.R 

# collate annotations
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$var <- paste(scores$POS, scores$REF, ">", scores$ALT, sep = "")
file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$var <- paste(file$POS, file$REF, ">", file$ALT, sep = "")
all_scores <- merge(file, scores[, c("var", "MLC_pos_score", "MLC_var_score", "pctRank_mean_OEUF")], by = "var", all.x = T)

# read in source data from Nature Communications manuscript
natcom <- read.delim(file = '../required_files/other_annotations/PMID37777527_Fig1_source_data.txt', header = TRUE, sep = "\t")
natcom$var <- paste(natcom$POS, natcom$REF, ">", natcom$ALT, sep = "")
all <- merge(all_scores, natcom[,c("var", "count_het", "count_homo")], by = "var", all.x = T)

# set missing values to 0
all[is.na(all$count_het), c("count_het")] <- "0"
all[is.na(all$count_homo), c("count_homo")] <- "0"

# collapse by site and variant type and merge
site_counts <- as.data.frame(all %>% group_by(POS, consequence) %>% summarize(count_het_site = sum(as.numeric(count_het))))
site_counts$match <- paste(site_counts$POS, site_counts$consequence)
all$match <- paste(all$POS, all$consequence)
all <- merge(all, site_counts[c("match", "count_het_site")], by = "match", all.x = T)

# the Y-axis of Tracks 1 and 2 for the circos plot is the log(number of participants with a heteroplasmy + 1), scaled from 0 to 9
all$het_track_uk <- log(as.numeric(all$count_het_site) + 1)

# set colors for plotting
all$color <- ifelse(grepl("intergenic", all$consequence), "#ffcc00", 
                    ifelse(grepl("MT-R", all$symbol), "#984ea3",
                           ifelse(grepl("MT-T", all$symbol), "#ff7f00", '#377eb8')))  # protein
all$color <- ifelse(grepl("stop_gain", all$consequence), "#b22222", 
                    ifelse(all$consequence == "synonymous_variant" | all$consequence == "synonymous_variant,synonymous_variant", "#4daf4a", as.character(all$color)))

# now, create data for "track 3" in Nature Communications -- adapted from their code

## files needed
df = data.frame(chr = "chrM", pos = seq(1, 16569))
test <- natcom
test$Median_AF_Het[is.na(test$Median_AF_Het)] <- 0

## identify mutational deserts, two ways 
## first with 5 zeros in a row
d1 = df[which(df$pos %ni% test$POS),] ## get list of zeros
nrow(d1)

d2 <- data.frame(seqToIntervals(d1$pos))
d2$diff <- d2$to - d2$from
table(d2$diff)

## make 3 categories, 0,1,2+ which correspond to 1,2,3+ mutational desert sites in a row
d2$col1 = paste(d2$from, "-", d2$to, sep="")
zero = subset(d2, diff == 0)
one = subset(d2, diff == 1)
more = subset(d2, diff > 1)

l1 <- lapply(strsplit(zero$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)
l2 <- lapply(strsplit(one$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)
l3 <- lapply(strsplit(more$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)

track3 <- d1
track3$chr <- "chrM"
track3$color <- ifelse(track3$pos %in% l2, "gray57", "gray87") # everything in l1 should be yellow
track3$color <- ifelse(track3$pos %in% l3, "gray28", track3$color)
track3$col1 <- ifelse(track3$color == "gray28", 1, 0.33)
track3$col1 <- ifelse(track3$color == "gray57", 0.66, as.double(track3$col1))


# now, generate the circos plot

circos.clear()
png("supplementary_figures/FigureS9a_circos.png", width = 4000, height = 4000, res = 600) 
all$chr <- "chrM"
circos.par("xaxis.clock.wise" = FALSE, "track.height" = 0.15, "start.degree" = 90, cell.padding = c(0.001, 0.001, 0.001, 0.001))
circos.initialize(factors=all$chr, x= all$POS)

# MLC - define a color gradient using colorRamp2
color_ramp <- colorRamp2(c(0, 0.5, 1.0), c("blue", "white", "red"))

circos.track(factors=all$chr, x = all$POS, y = all$MLC_pos_score,
             panel.fun = function(x, y) {
               circos.points(x, y, col = color_ramp(y), cex = 0.07)
             })

locations=c(1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16569)

for (i in 1:length(locations)){
  circos.axis(labels.cex = 0.8, major.at = locations)
}

# do the chimpanzee track
all$h <- 1.5
all$chimp_match <- ifelse(as.character(all$REF) == as.character(all$chimp_ref), "yes", "no")
all$color_chimp <- ifelse(
  all$chimp_match == "no" & grepl("intergenic|non_coding_transcript|missense|incomplete_terminal", all$consequence), "white", "black")
chimp <- unique(all[all$chimp_match == "no" & grepl("intergenic|non_coding_transcript|missense|incomplete_terminal", all$consequence), 
                    c("chr", "POS", "h", "color_chimp")])
scale_row <- data.frame(chr = "chrM", POS = c(1, 2), h = c(1, 0), color_chimp = "white")
meep <- rbind(chimp, scale_row)
chimp_scaled <- meep

circos.track(factors=chimp_scaled$chr, x = chimp_scaled$POS, y = chimp_scaled$h,
             panel.fun = function(x, y) {
               circos.lines(x, y, type = 'h', col = chimp_scaled$color_chimp, lwd = 0.75)
             }, bg.col = "black", bg.border = NA, track.height = 0.03)


# Track 1
track1 <- unique(all[all$count_het > 0 & (grepl("intergenic|non_coding_transcript", all$consequence) | all$consequence == "synonymous_variant" | all$consequence == "synonymous_variant,synonymous_variant"), 
                     c("chr", "POS", "het_track_uk", "color")])
scale_row <- data.frame(chr = "chrM", POS = c(1,2), het_track_uk = c(8,0), color = "white")
meep <- rbind(track1, scale_row)
track1_scaled <- meep

circos.track(factors=track1_scaled$chr, x = track1_scaled$POS, y = track1_scaled$het_track_uk,
             panel.fun = function(x, y) {
               circos.lines(x, y, type = 'h',col = track1_scaled$color, lwd = 0.75)
             })

# Track 2 (order to have stop gain last)
track2 <- unique(all[all$count_het > 0 & grepl("missense|stop_gained", all$consequence),
                     c("chr", "POS", "het_track_uk", "color")])
track2 <- track2[order(track2$color),]
meep <- rbind(track2, scale_row)
track2_scaled <- meep

circos.track(factors=track2_scaled$chr, x = track2_scaled$POS, y = track2_scaled$het_track_uk,
             panel.fun = function(x, y) {
               circos.lines(x, y, type = 'h', col = track2_scaled$color, lwd = 0.75)
             })

# Track 3 to scale
scale_row <- data.frame(chr = "chrM", pos = c(16565, 16570), col1 = c(0, 1.33), color = "white")
meep <- rbind(track3, scale_row)
track3_scaled <- meep

circos.track(factors=track3_scaled$chr, x = track3_scaled$pos, y = track3_scaled$col1,
             panel.fun = function(x, y) {
               circos.lines(x, y, type = 'h', col = track3_scaled$color, lwd = 0.5)
             }, track.height = 0.1)

dev.off()


# Figure S9b-c - histogram showing MLC distribution for UK Biobank data

# first, collapse just by site (nucleotide position) and merge for per position total counts (both het and hom) -- site counts just by position
all$total_ukb <- as.numeric(all$count_het) + as.numeric(all$count_homo)
site_counts <- as.data.frame(all %>% group_by(POS) %>% summarize(count_site = sum(as.numeric(total_ukb))))
all <- merge(all, site_counts[c("POS", "count_site")], by = "POS", all.x = T)

# variants at heteroplasmy seen only once in UKB
plotB <- ggplot(data = all[all$count_het == 1 & !grepl(",|&", all$consequence),], aes(as.numeric(MLC_var_score))) + 
  geom_histogram(binwidth = 0.025, aes(x = MLC_var_score, fill = ..x..), color = "black", size = 0.25) +
  paper_theme +
  scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1)) +
  scale_y_continuous(trans = "sqrt") +
  labs(x = "MLC variant score", y = "Variants seen once\nat heteroplasmy in UKB") +
  theme(legend.position="none",
        plot.margin = unit(c(0.25, 0.25, 0.25, 1.0), "cm"))

# sites with no variants in the UKB
plotC <- ggplot(data = all[all$count_site == 0 & !grepl(",|&", all$consequence),], aes(as.numeric(MLC_var_score))) + 
  geom_histogram(binwidth = 0.025, aes(x = MLC_var_score, fill = ..x..), color = "black", size = 0.25) +
  paper_theme +
  scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", space = "Lab", limits = c(0, 1)) +
  scale_y_continuous(trans = "sqrt") +
  labs(x = "MLC variant score", y = "Variants not in UKB\nat invariant sites") +
  theme(legend.position="none",
        plot.margin = unit(c(0.25, 1.0, 0.25, 0.25), "cm"))

# combine with circos plot
fig <- readPNG("supplementary_figures/FigureS9a_circos.png")  
plotA <- ggplot() + 
  background_image(fig) +
  theme(panel.background = element_rect(fill = 'white', color = 'white'))


# compile figure
ggarrange(plotA, 
          ggarrange(plotB, plotC, nrow = 1, labels = c("b", "c"), font.label = list(size = 10)),
          nrow = 2, heights = c(1, 0.25), labels = c("a", ""), font.label = list(size = 10))

ggsave("supplementary_figures/FigureS9.jpeg", width = 200, height = 250, dpi = 600, units = c("mm"))

