# This script creates Extended Data Figure 11. 

rm(list=ls())

library(ggplot2)
library(dplyr)
library(fixest)

setwd("data/figure_and_input_data/")

# Code for 11a 
ug_dhs_cell_ownership <- read.csv("ug_dhs_cellular_ownership.csv")

pdf("../../figures/raw/EDFigure11a.pdf", width = 6, height = 4)
EDFigure11a <- ggplot() +
  geom_line(data = ug_dhs_cell_ownership, aes(x = year, y = cluster_ave, group = treated, color = factor(treated))) +
  scale_color_manual(values = c("black", "red")) +
  ylab("% Ownership") +
  xlab("") +
  ylim(0,1)  +
  theme_clean() +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.position = "none")
EDFigure11a
dev.off()

# ED Figures 11 b and c

cell_coverage <- read.csv("cell_coverage.csv")

pdf("../../figures/raw/EDFigure11b.pdf", width = 5, height = 4)
EDFigure11b <- ggplot() +
  geom_line(data = cell_coverage, aes(x = year, y = ave_coverage, group = treated, color = factor(treated))) +
  scale_color_manual(values = c("black", "red")) +
  ylab("% Coverage") +
  xlab("") +
  ylim(0,1)  +
  theme_clean() +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.position = "none")
EDFigure11b
dev.off()


cell_quality <- read.csv("cell_quality.csv")

pdf("../../figures/raw/EDFigure11c.pdf", width = 5, height = 4)
EDFigure11c <- ggplot() +
  geom_line(data = cell_quality, aes(x = year, y = ave_coverage, group = treated, color = factor(treated))) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Quality of Coverage") +
  xlab("") +
  ylim(0,1)  +
  theme_clean() +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        axis.title = element_text(size = 18), axis.text = element_text(size = 18), 
        legend.position = "none")
EDFigure11c
dev.off()
