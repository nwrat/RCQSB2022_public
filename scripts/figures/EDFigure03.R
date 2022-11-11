# This code creates Extended Data Figure 3, DHS varibles used in wealth index creation. 

rm(list=ls())

library(gridExtra)

setwd("data/figure_and_input_data/")

dhs_varible_list <- read.csv("WI_DHS_variable_list.csv", check.names=FALSE)

pdf("../../figures/raw/dhs_varible_list.pdf", height=6, width=5)
grid.table(dhs_varible_list)
dev.off()
