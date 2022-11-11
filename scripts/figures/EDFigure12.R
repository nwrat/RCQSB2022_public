# This is the script for the pretrends table, Extended Data Figure 12. 

library(gridExtra)

setwd("data/figure_and_input_data/")

pretrend_table <- read.csv("pretrends_test_WI.csv")

colnames(pretrend_table) <- c("Penalty", "With UG, 95%", "With Ug, 90%", "Without UG, 95%", "Without UG, 90%")

pdf("../../figures/raw/EDFigure12.pdf")
grid.table(pretrend_table)
dev.off()
