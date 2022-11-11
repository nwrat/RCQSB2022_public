# This script creates ED Figure 5. 

rm(list=ls())

library(ggplot2)
library(dplyr)
library(fixest)

setwd("data/figure_and_input_data/")

cnn_output_ug <- read.csv("cnn_output_ug.csv")
cnn_output <- read.csv("cnn_output.csv")
beta_r2 <- read.csv("beta_r2.csv") # aggregated outcomes from linear regressions shown below. 

###########
# this creates the scatter plots and outputs shown in a-e
###########

summary(lm(zero ~ pred, data=cnn_output_ug))

pdf("../../figures/raw/EDFigure05a.pdf", width=4, height=4)
EDFigure05a <- ggplot(data = cnn_output_ug, aes(x = pred, y = zero)) +
  annotate("segment", x = -6, xend = 6, y = -6, yend = 6, colour = "red") +
  geom_point(alpha = .15) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black") +
  ylim(-6,6) +
  xlim(-6,6) +
  ylab("Penalized Values") +
  xlab("") +
  theme_light() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"),
        axis.ticks.x = element_blank(),       
        axis.text.x = element_blank()) 
EDFigure05a
dev.off()

summary(lm(one ~ pred, data=cnn_output_ug))

pdf("../../figures/raw/EDFigure05b.pdf", width=4, height=4)
EDFigure05b <- ggplot(data = cnn_output_ug, aes(x = pred, y = one)) +
  annotate("segment", x = -6, xend = 6, y = -6, yend = 6, colour = "red") +
  geom_point(alpha = .15) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black") +
  ylim(-6,6) +
  xlim(-6,6) +
  ylab("") +
  xlab("") +
  theme_light() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"),
        axis.ticks.x = element_blank(),       
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank()) 
EDFigure05b
dev.off()

summary(lm(three ~ pred, data=cnn_output_ug))

pdf("../../figures/raw/EDFigure05c.pdf", width=4, height=4)
EDFigure05c <- ggplot(data = cnn_output_ug, aes(x = pred, y = three)) +
  annotate("segment", x = -6, xend = 6, y = -6, yend = 6, colour = "red") +
  geom_point(alpha = .15) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black") +
  ylim(-6,6) +
  xlim(-6,6)+
  ylab("") +
  xlab("") +
  theme_light() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"),
        axis.ticks.x = element_blank(),       
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank()) 
EDFigure05c
dev.off()

summary(lm(five ~ pred, data=cnn_output_ug))

pdf("../../figures/raw/EDFigure05d.pdf", width=4, height=4)
EDFigure05d <- ggplot(data = cnn_output_ug, aes(x = pred, y = five)) +
  annotate("segment", x = -6, xend = 6, y = -6, yend = 6, colour = "red") +
  geom_point(alpha = .15) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black") +
  ylim(-6,6) +
  xlim(-6,6)+
  ylab("Penalized Values") +
  xlab("DHS Observed Values") +
  theme_light() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"))
EDFigure05d
dev.off()

summary(lm(seven ~ pred, data=cnn_output_ug))

pdf("../../figures/raw/EDFigure05e.pdf", width=4, height=4)
EDFigure05e <- ggplot(data = cnn_output_ug, aes(x = pred, y = seven)) +
  annotate("segment", x = -6, xend = 6, y = -6, yend = 6, colour = "red") +
  geom_point(alpha = .15) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black") +
  ylim(-6,6) +
  xlim(-6,6) +
  ylab("") +
  xlab("DHS Observed Values") +  theme_light() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank()) 
EDFigure05e
dev.off()

###########
# this creates the regression output for SSA in panel f
###########

summary(lm(zero ~ pred, data=cnn_output))
summary(lm(one ~ pred, data=cnn_output))
summary(lm(three ~ pred, data=cnn_output))
summary(lm(five ~ pred, data=cnn_output))
summary(lm(seven ~ pred, data=cnn_output))

pdf("../../figures/raw/EDFigure05f.pdf", width=4, height=4)
EDFigure05f <- ggplot(data = beta_r2) +
  annotate("segment", x = 1, xend = 1, y = .6, yend = .66, colour = "red", size = .25) +
  geom_line(aes(x = ssa_beta, y= ssa_r2), color = "gold", linetype = "dotted", size = .75) +
  geom_point(aes(x = ssa_beta, y= ssa_r2), color = "gold", size = 3, shape = 15) +
  geom_line(aes(x = ug_beta, y= ug_r2), linetype = "dotted", size = .75) +
  geom_point(aes(x = ug_beta, y= ug_r2), size = 3, shape = 15) +
  geom_text(aes(x = ug_beta, y= ug_r2, label=pen),hjust=-.5, vjust=-1, size = 5) +
  scale_y_continuous(position = "right", limits = c(.6,.68)) +
  xlim(.65,1.1) +
  ylab("r2") +
  xlab("Regression Coefficient") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), plot.background = element_rect(color = "white"))
EDFigure05f
dev.off()



