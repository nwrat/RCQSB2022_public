# This code creates Extended Date Figure 2, DHS surveys used in this study. 

rm(list=ls())

library(gridExtra)

setwd("data/figure_and_input_data/")

dhs_survey_list <- read.csv("WI_DHS_survey_list.csv", check.names=FALSE)

pdf("../../figures/raw/dhs_survey_list.pdf", height=18, width=5)
grid.table(dhs_survey_list)
dev.off()
