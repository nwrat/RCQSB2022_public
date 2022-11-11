# This code creates Extended Data Figure 4, a-c.   

rm(list=ls())

library(ggplot2)

setwd("data/figure_and_input_data/")

pca_output <- read.csv("WI_pca_output.csv", check.names=FALSE)
pca_output <- as.data.frame(pca_output)

# Panel A
pdf("../../figures/raw/base_versus_w_elec_scatter.pdf") 
base_versus_w_elec_scatter <- ggplot(pca_output, aes(x=base, y=base_w_elec)) +
  geom_point(alpha = .15) +
  ylab("Index w/ 'electricity' variable") +
  xlab("Base Wealth Index") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title = element_text(size = 22), axis.text = element_text(size = 22))+
  ylim(-5,6.5) +
  xlim(-5,6.5)
base_versus_w_elec_scatter
dev.off()

base_reg <- lm(base_w_elec ~ base, data = pca_output)
summary(base_reg) 

# Panel B
pdf("../../figures/raw/base_versus_elec_no_app_scatter.pdf") 
base_versus_elec_no_app_scatter <- ggplot(pca_output, aes(x=base, y=elec_no_app)) +
  geom_point(alpha = .15) +
  ylab("Index w/ elec. but w/o elec appliances") +
  xlab("Base Wealth Index") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title = element_text(size = 22), axis.text = element_text(size = 22)) +
  ylim(-5,6.5) +
  xlim(-5,6.5)
base_versus_elec_no_app_scatter
dev.off()

reg_elec_no_app <- lm(elec_no_app ~ base, data = pca_output)
summary(reg_elec_no_app) 


# Panel C
pdf("../../figures/raw/base_versus_no_elec_no_app_scatter.pdf") 
base_versus_no_elec_no_app_scatter <- ggplot(pca_output, aes(x=base, y=no_elec_no_app)) +
  geom_point(alpha = .15) +
  ylab("Index w/o elec. or elec. appliances") +
  xlab("Base Wealth Index") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", axis.title = element_text(size = 22), axis.text = element_text(size = 22)) +
  ylim(-5,6.5) +
  xlim(-5,6.5)
base_versus_no_elec_no_app_scatter
dev.off() 

base_reg_no_app <- lm(no_elec_no_app ~ base, data = pca_output)
summary(base_reg_no_app) 

