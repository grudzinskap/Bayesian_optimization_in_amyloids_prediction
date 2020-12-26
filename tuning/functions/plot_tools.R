library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(ggrepel)

size_mod <- -6

my_theme <- theme(plot.background=element_rect(fill = "transparent",
                                               color = "transparent"),
                  #panel.grid.major = element_line(color="grey", linetype = "dashed", size = 0.5),
                  panel.grid.major = element_line(color="lightgrey", linetype = "dashed", size = 0.5),
                  panel.background = element_rect(fill = "transparent", color = "black"),
                  legend.background = element_rect(fill = "NA"),
                  legend.position = "bottom",
                  legend.key = element_rect(fill = "NA", color = "NA", size = 0.5),
                  legend.key.size = unit(0.1, "inches"),
                  legend.spacing = unit(-0.25, "lines"),
                  axis.text = element_text(size=13 + size_mod),
                  axis.title.x = element_text(size=14 + size_mod, vjust = -1), 
                  axis.title.y = element_text(size=14 + size_mod, vjust = 1),
                  strip.text = element_text(size=13 + size_mod, face = "bold"),
                  legend.text = element_text(size=13 + size_mod), 
                  legend.title = element_text(size=14 + size_mod),
                  plot.title = element_text(size=20 + size_mod),
                  strip.background = element_rect(fill = "NA", color = "NA"),
                  plot.margin = unit(rep(0.02, 4), "inches"))
