# THIS SCRIPT CREATES FIGURE 3
# If you attempt to plot these figures in your R window instead of to a PDF, as written, 
# you may need to mute line "theme_clean" with a hashtah, ie. "#theme_clean". 

rm(list=ls())

library(ggplot2)
library(ggthemes)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

c_year_means <- read.csv("c_year_means.csv")
t_year_means <- read.csv("t_year_means.csv")
c_se <- read.csv("c_se.csv")
t_se <- read.csv("t_se.csv")
estimates_2km_5pl_2011 <- read.csv("estimates_2km_5pl_2011_aug.csv")
estimates_2km_5pl_2012 <- read.csv("estimates_2km_5pl_2012_aug.csv")
estimates_2km_5pl_2013 <- read.csv("estimates_2km_5pl_2013_aug.csv")
estimates_2km_5pl_2014 <- read.csv("estimates_2km_5pl_2014_aug.csv")
estimates_2km_5pl_2015 <- read.csv("estimates_2km_5pl_2015_aug.csv")
estimates_2km_5pl <- read.csv("estimates_2km_5pl_aug.csv")
estimates_2km_5pl_noUG <- read.csv("estimates_2km_5pl_noUG_aug.csv")
estimates_2km_5pl_raster <- read.csv("estimates_2km_5pl_raster_n500.csv")
idw_results_df <- read.csv("idw_results_df.csv")
two_unit <- read.csv("two_unit_results_df.csv")

figure03a_df <- as.data.frame(matrix(NA,11,7))
figure03a_df[,1] <- c(2006:2016)
figure03a_df[1:11,2] <- c_year_means[1:11,2]
figure03a_df[1:11,3] <- t_year_means[1:11,2]
figure03a_df[1:11,4] <- c_se[1:11,3]
figure03a_df[1:11,5] <- t_se[1:11,3]
figure03a_df[6,6] <- figure03a_df[6,3] - estimates_2km_5pl_2011[,1]
figure03a_df[7,6] <- figure03a_df[7,3] - estimates_2km_5pl_2012[,1]
figure03a_df[8,6] <- figure03a_df[8,3] - estimates_2km_5pl_2013[,1]
figure03a_df[9,6] <- figure03a_df[9,3] - estimates_2km_5pl_2014[,1]
figure03a_df[10,6] <- figure03a_df[10,3] - estimates_2km_5pl_2015[,1]
figure03a_df[11,6] <- figure03a_df[11,3] - estimates_2km_5pl[,1] # MC
figure03a_df[6,7] <- figure03a_df[6,3] -  estimates_2km_5pl_2011[,4] # sc
figure03a_df[7,7] <- figure03a_df[7,3] - estimates_2km_5pl_2012[,4]
figure03a_df[8,7] <- figure03a_df[8,3] - estimates_2km_5pl_2013[,4]
figure03a_df[9,7] <- figure03a_df[9,3] - estimates_2km_5pl_2014[,4]
figure03a_df[10,7] <- figure03a_df[10,3] - estimates_2km_5pl_2015[,4]
figure03a_df[11,7] <- figure03a_df[11,3] - estimates_2km_5pl[,4]

colnames(figure03a_df) <- c("year", "control", "treated","c_se", "t_se", "mc", "sc")

# Figure 3 a
pdf("Figure03a.pdf", width=8, height=6) 
Figure03a <- ggplot(data = figure03a_df, aes(x = year)) +
  geom_line(aes( y = mc), color = "springgreen", group = 1, size = .65, linetype = "dashed") +
  geom_point(aes( y = mc), color = "springgreen", group = 1, size = 1.5) +
  geom_line(aes( y = sc), color = "dodgerblue", group = 1, size = .65, linetype = "dashed") +
  geom_point(aes( y = sc), color = "dodgerblue", group = 1, size = 1.5) +
  geom_line(aes( y = control), color = "grey50", group = 1, size = 1) +
  geom_ribbon(aes(ymin = control - c_se, ymax = control + c_se), fill = "grey50", color = NA, alpha = .25) +
  geom_point(aes(y = control), color = "grey50", size = 1.5) +
  geom_line(aes( y = treated), color = "tomato", group = 1, size = 1) +
  geom_ribbon(aes(ymin = treated - t_se, ymax = treated + t_se), fill = "tomato", color = NA, alpha = .25) +
  geom_point(aes( y = treated), color = "tomato", size = 1.5) +
  ylim(-2.5,-.5) +
  ylab("Wealth Index") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18), plot.background = element_rect(color = "white"))
Figure03a
dev.off()

# 3b
figure03b_df <- cbind(estimates_2km_5pl, estimates_2km_5pl_noUG, estimates_2km_5pl_raster, 
                      estimates_2km_5pl_2013[,1:3], 
                      estimates_2km_5pl_2014[,1:3],
                      estimates_2km_5pl_2015[,1:3],
                      idw_results_df, two_unit) 

colnames(figure03b_df) <- c("mc_2016", "mc_2016_low", "mc_2016_high",
                            "sc_2016", "sc_2016_low", "sc_2016_high",
                            "mc_2016_wo", "mc_2016_low_wo", "mc_2016_high_wo",
                            "sc_2016_wo", "sc_2016_low_wo", "sc_2016_high_wo",
                            "mc_2016_r", "mc_2016_low_r", "mc_2016_high_r",
                            "sc_2016_r", "sc_2016_low_r", "sc_2016_high_r",
                            "mc_2013", "mc_2013_low", "mc_2013_high",
                            "mc_2014", "mc_2014_low", "mc_2014_high",
                            "mc_2015", "mc_2015_low", "mc_2015_high",
                            "idw", "idw_se",
                            "two_unit", "two_unit_se")

#Figure 3b
pdf("Figure03b.pdf", width=8, height=6) 
Figure03b <- ggplot(data = figure03b_df) +
  geom_point(aes(x = .3, y = mc_2016), color = "springgreen", size = 2) +
  geom_segment(aes(x = .3, y = mc_2016_low, xend = .3, yend = mc_2016_high), color = "springgreen", size = .2) +
  geom_point(aes(x = .35, y = sc_2016), color = "dodgerblue", size = 2) +
  geom_segment(aes(x = .35, y = sc_2016_low, xend = .35, yend = sc_2016_high), color = "dodgerblue", size = .2) +
  
  geom_point(aes(x = .65, y = mc_2016_wo), color = "springgreen", shape = 17, size = 2) +
  geom_segment(aes(x = .65, y = mc_2016_low_wo, xend = .65, yend = mc_2016_high_wo), color = "springgreen", size = .2) +
  geom_point(aes(x = .7, y = sc_2016_wo), color = "dodgerblue", shape = 17, size = 2) +
  geom_segment(aes(x = .7, y = sc_2016_low_wo, xend = .7, yend = sc_2016_high_wo), color = "dodgerblue", size = .2) +
  
  geom_point(aes(x = 1, y = mc_2016_r), color = "springgreen", shape = 15, size = 2) +
  geom_segment(aes(x = 1, y = mc_2016_low_r, xend = 1, yend = mc_2016_high_r), color = "springgreen", size = .2) +
  geom_point(aes(x = 1.05, y = sc_2016_r), color = "dodgerblue", shape = 15, size = 2) +
  geom_segment(aes(x = 1.05, y = sc_2016_low_r, xend = 1.05, yend = sc_2016_high_r), color = "dodgerblue", size = .2) +
  
  geom_point(aes(x = 1.25, y = two_unit), color = "black", shape = 3, size = 2) +
  geom_segment(aes(x = 1.25, y = two_unit-(1.96*two_unit_se), xend = 1.25, yend = two_unit+(1.96*two_unit_se)), color = "black", size = .2) +
  
  geom_point(aes(x = 1.35, y = idw), color = "black", shape = 3, size = 2) +
  geom_segment(aes(x = 1.35, y = idw-(1.96*idw_se), xend = 1.35, yend = idw+(1.96*idw_se)), color = "black", size = .2) +
  
  geom_point(aes(x = 1.6, y = mc_2013), color = "springgreen", size = 2) +
  geom_segment(aes(x = 1.6, y = mc_2013_low, xend = 1.6, yend = mc_2013_high), color = "springgreen", size = .2) +
  geom_point(aes(x = 1.7, y = mc_2014), color = "springgreen", size = 2) +
  geom_segment(aes(x = 1.7, y = mc_2014_low, xend = 1.7, yend = mc_2014_high), color = "springgreen", size = .2) +
  
  geom_point(aes(x = 1.8, y = mc_2015), color = "springgreen", size = 2) +
  geom_segment(aes(x = 1.8, y = mc_2015_low, xend = 1.8, yend = mc_2015_high), color = "springgreen", size = .2) +
  geom_point(aes(x = 1.9, y = mc_2016), color = "springgreen", size = 2) +
  geom_segment(aes(x = 1.9, y = mc_2016_low, xend = 1.9, yend = mc_2016_high), color = "springgreen", size = .2) +
  
  ylim(-.25,.75) +
  xlim(0.2,2) +
  ylab("Estimated Causal Effect") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), plot.background = element_rect(color = "white"))
Figure03b
dev.off()

