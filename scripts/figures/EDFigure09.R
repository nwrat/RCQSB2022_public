# THIS SCRIPT CREATES ED FIGURE09
# If you attempt to plot these figures in your R window instead of to a PDF, as written, 
# you may need to mute line "theme_clean" with a hashtah, ie. "#theme_clean". 

rm(list=ls())

library(ggplot2)
library(ggthemes)

setwd("data/figure_and_input_data/")

est_2km_0pl <- read.csv("estimates_2km_0pl_aug.csv")
est_2km_1pl <- read.csv("estimates_2km_1pl_aug.csv")
est_2km_3pl <- read.csv("estimates_2km_3pl_aug.csv")
est_2km_7pl <- read.csv("estimates_2km_7pl_aug.csv")

est_2km_0pl_noUG <- read.csv("estimates_2km_0pl_noUG_aug.csv")
est_2km_1pl_noUG <- read.csv("estimates_2km_1pl_noUG_aug.csv")
est_2km_3pl_noUG <- read.csv("estimates_2km_3pl_noUG_aug.csv")
est_2km_7pl_noUG <- read.csv("estimates_2km_7pl_noUG_aug.csv")

est_2km_5pl <- read.csv("estimates_2km_5pl_aug.csv")
est_2km_5pl_noUG <- read.csv("estimates_2km_5pl_noUG_aug.csv")


#ED Figure09
pdf("../../figures/raw/EDFigure09.pdf", width=8, height=6) 
EDFigure09 <- ggplot() +
  geom_point(data = est_2km_0pl, aes(x = .3, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_0pl, aes(x = .3, y = mc_l95, xend = .3, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_0pl, aes(x = .35, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_0pl, aes(x = .35, y = sc_l95, xend = .35, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_1pl, aes(x = .65, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_1pl, aes(x = .65, y = mc_l95, xend = .65, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_1pl, aes(x = .7, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_1pl, aes(x = .7, y = sc_l95, xend = .7, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_3pl, aes(x = 1, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_3pl, aes(x = 1, y = mc_l95, xend = 1, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_3pl, aes(x = 1.05, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_3pl, aes(x = 1.05, y = sc_l95, xend = 1.05, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_5pl, aes(x = 1.35, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_5pl, aes(x = 1.35, y = mc_l95, xend = 1.35, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_5pl, aes(x = 1.4, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_5pl, aes(x = 1.4, y = sc_l95, xend = 1.4, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_7pl, aes(x = 1.7, y = mc), color = "springgreen", size = 2) +
  geom_segment(data = est_2km_7pl, aes(x = 1.7, y = mc_l95, xend = 1.7, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_7pl, aes(x = 1.75, y = sc), color = "dodgerblue", size = 2, ) +
  geom_segment(data = est_2km_7pl, aes(x = 1.75, y = sc_l95, xend = 1.75, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  
  geom_point(data = est_2km_0pl_noUG, aes(x = 2.35, y = mc), color = "springgreen", size = 2, shape = 17) +
  geom_segment(data = est_2km_0pl_noUG, aes(x = 2.35, y = mc_l95, xend = 2.35, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_0pl_noUG, aes(x = 2.4, y = sc), color = "dodgerblue", size = 2, shape = 17) +
  geom_segment(data = est_2km_0pl_noUG, aes(x = 2.4, y = sc_l95, xend = 2.4, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_1pl_noUG, aes(x = 2.7, y = mc), color = "springgreen", size = 2, shape = 17) +
  geom_segment(data = est_2km_1pl_noUG, aes(x = 2.7, y = mc_l95, xend = 2.7, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_1pl_noUG, aes(x = 2.75, y = sc), color = "dodgerblue", size = 2, shape = 17 ) +
  geom_segment(data = est_2km_1pl_noUG, aes(x = 2.75, y = sc_l95, xend = 2.75, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_3pl_noUG, aes(x = 3.05, y = mc), color = "springgreen", size = 2, shape = 17) +
  geom_segment(data = est_2km_3pl_noUG, aes(x = 3.05, y = mc_l95, xend = 3.05, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_3pl_noUG, aes(x = 3.1, y = sc), color = "dodgerblue", size = 2, shape = 17 ) +
  geom_segment(data = est_2km_3pl_noUG, aes(x = 3.1, y = sc_l95, xend = 3.1, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_5pl_noUG, aes(x = 3.4, y = mc), color = "springgreen", size = 2, shape = 17) +
  geom_segment(data = est_2km_5pl_noUG, aes(x = 3.4, y = mc_l95, xend = 3.4, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_5pl_noUG, aes(x = 3.45, y = sc), color = "dodgerblue", size = 2, shape = 17 ) +
  geom_segment(data = est_2km_5pl_noUG, aes(x = 3.45, y = sc_l95, xend = 3.45, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  geom_point(data = est_2km_7pl_noUG, aes(x = 3.75, y = mc), color = "springgreen", size = 2, shape = 17) +
  geom_segment(data = est_2km_7pl_noUG, aes(x = 3.75, y = mc_l95, xend = 3.75, yend = mc_h95), color = "springgreen", size = .2) +
  geom_point(data = est_2km_7pl_noUG, aes(x = 3.8, y = sc), color = "dodgerblue", size = 2, shape = 17 ) +
  geom_segment(data = est_2km_7pl_noUG, aes(x = 3.8, y = sc_l95, xend = 3.8, yend = sc_h95), color = "dodgerblue", size = .2) +
  
  ylim(-.25,.75) +
  ylab("Estimated Causal Effect") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.background = element_rect(color = "white"))
EDFigure09
dev.off()






