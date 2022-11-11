# THIS SCRIPT CREATES ED FIGURE 08
# If you attempt to plot these figures in your R window instead of to a PDF, as written, 
# you may need to mute line "theme_clean" with a hashtah, ie. "#theme_clean". 

rm(list=ls())

library(ggplot2)
library(ggthemes)

setwd("data/figure_and_input_data/")

# 
est_2km_5pl <- read.csv("estimates_2km_5p.csv")
est_3km_5pl <- read.csv("estimates_3km_5p.csv")
est_4km_5pl <- read.csv("estimates_4km_5p.csv")

#Figure 3b
pdf("../../figures/raw/EDFigure08.pdf", width=8, height=6) 
EDFigure08 <- ggplot() +
  geom_point(data = est_2km_5pl, aes(x = .3, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_5pl, aes(x = .3, y = mc_l95, xend = .3, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_5pl, aes(x = .35, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_5pl, aes(x = .35, y = sc_l95, xend = .35, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_3km_5pl, aes(x = .65, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_3km_5pl, aes(x = .65, y = mc_l95, xend = .65, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_3km_5pl, aes(x = .7, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_3km_5pl, aes(x = .7, y = sc_l95, xend = .7, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_4km_5pl, aes(x = 1, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_4km_5pl, aes(x = 1, y = mc_l95, xend = 1, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_4km_5pl, aes(x = 1.05, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_4km_5pl, aes(x = 1.05, y = sc_l95, xend = 1.05, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  ylim(-.25,.75) +
  xlim(0.1,1.25) +
  ylab("Estimated Causal Effect") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.background = element_rect(color = "white"))
EDFigure08
dev.off()






