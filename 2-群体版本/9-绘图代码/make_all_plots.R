# generate all figures
dir.create("figures")
dir.create("extended_data_figures")
dir.create("supplementary_figures")
dir.create("supplementary_datasets")
dir.create("final_figures_source_data")

library(ggplot2)

# ggplot theme to use for manuscript
paper_theme <- theme(plot.title = element_text(size = 8),
                     axis.title.x = element_text(size = 8), 
                     axis.text.x  = element_text(size = 8), 
                     axis.title.y = element_text(size = 8),
                     axis.text.y  = element_text(size = 8),
                     axis.ticks.length = unit(0.15, "cm"),
                     panel.background = element_blank(), 
                     plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
                     axis.line = element_line(colour = "black"),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 8), 
                     strip.text.x = element_text(size = 8),
                     strip.text.y = element_text(size = 8))

# supplementary datasets
source("Supplementary_Datasets.R")

# need this for Figure2
source("FigureED4.R", echo = TRUE) 

# main text figures
source("Figure1.R", echo = TRUE) 
source("Figure2.R", echo = TRUE)
source("Figure3.R", echo = TRUE)
source("Figure4.R", echo = TRUE)
source("Figure5.R", echo = TRUE)

# extended data figures
# FigureED1 assembled separately
source("FigureED2.R", echo = TRUE)
source("FigureED3.R", echo = TRUE)
source("FigureED4.R", echo = TRUE) 
source("FigureED5.R", echo = TRUE) 
source("FigureED6.R", echo = TRUE) 
# FigureED7 assembled separately as svg
source("FigureED8.R", echo = TRUE)
source("FigureED9.R", echo = TRUE)
source("FigureED10.R", echo = TRUE) 

# supplementary figures
source("FigureS1.R", echo = TRUE) 
source("FigureS2.R", echo = TRUE) 
source("FigureS3.R", echo = TRUE) 
source("FigureS4.R", echo = TRUE) 
source("FigureS5.R", echo = TRUE)
# FigureS6 assembled separately as illustration
source("FigureS7.R", echo = TRUE)
source("FigureS8.R", echo = TRUE)
source("FigureS9.R", echo = TRUE)
source("FigureS10.R", echo = TRUE)

# source data
source("compile_source_data.R")