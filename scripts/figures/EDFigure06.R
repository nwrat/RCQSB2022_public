# This code creates Extended Data Figure 6 a-c. 

rm(list=ls())

library(ggplot2)
library(ggthemes)

setwd("data/figure_and_input_data/")

#Extended Date Figure 6a
dd_mc_output_1alpha <- read.csv("dd_mc_output_1a_test.csv")
dd_mc_output_15alpha <- read.csv("dd_mc_output_15a_test.csv")
dd_mc_output_2alpha <- read.csv("dd_mc_output_2a_test.csv")
dd_mc_output_25alpha <- read.csv("dd_mc_output_25a_test.csv")

pdf("../../figures/raw/ED_Fig_6a.pdf", width=6, height=6) 
ED_Fig_6a <- ggplot() +
  geom_line(data = dd_mc_output_1alpha, aes(x=Year, y=DID_diff), color = "dodgerblue1", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_15alpha, aes(x=Year, y=DID_diff), color = "dodgerblue2", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_2alpha, aes(x=Year, y=DID_diff), color = "dodgerblue3", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_25alpha, aes(x=Year, y=DID_diff), color = "dodgerblue4", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_1alpha, aes(x=Year, y=DD_np_diff), color = "dodgerblue1") +
  geom_line(data = dd_mc_output_15alpha, aes(x=Year, y=DD_np_diff), color = "dodgerblue2") +
  geom_line(data = dd_mc_output_2alpha, aes(x=Year, y=DD_np_diff), color = "dodgerblue3") +
  geom_line(data = dd_mc_output_25alpha, aes(x=Year, y=DD_np_diff), color = "dodgerblue4") +
  geom_line(data = dd_mc_output_1alpha, aes(x=Year, y=MC_diff), color = "orangered1") +
  geom_line(data = dd_mc_output_15alpha, aes(x=Year, y=MC_diff), color = "orangered2") +
  geom_line(data = dd_mc_output_2alpha, aes(x=Year, y=MC_diff), color = "orangered3") +
  geom_line(data = dd_mc_output_25alpha, aes(x=Year, y=MC_diff), color = "orangered4") +
  ylim(0,4.5) +
  ylab("Difference from True ATE") +
  xlab("Years in Timeseries") +
  #theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 18), axis.text = element_text(size = 18), plot.background = element_rect(color = "white"))
ED_Fig_6a
dev.off()

#Extended Date Figure 6b
dta_small <- read.csv("dta_small_sim.csv")
dta_med <- read.csv("dta_med_sim.csv")
dta_big <- read.csv("dta_big_sim.csv")
dta_biggest <- read.csv("dta_biggest_sim.csv")

pdf("../../figures/raw/ED_Fig_6b.pdf", width=6, height=6)
ED_Fig_6b <- ggplot() + 
  geom_density(data = dta_small, aes(x = y), color = "grey70") + 
  geom_density(data = dta_small, aes(x = y_b), color = "green1") +
  geom_density(data = dta_med, aes(x = y), color = "grey80") + 
  geom_density(data = dta_med, aes(x = y_b), color = "green2") +
  geom_density(data = dta_big, aes(x = y), color = "grey90") + 
  geom_density(data = dta_big, aes(x = y_b), color = "green3") +
  geom_density(data = dta_biggest, aes(x = y), color = "black") + 
  geom_density(data = dta_biggest, aes(x = y_b), color = "green4") +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 18), axis.text = element_text(size = 18), plot.background = element_rect(color = "white"), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank())
ED_Fig_6b
dev.off()


#Extended Date Figure 6c
dd_mc_output_1alpha_small <- read.csv("dd_mc_output_1alpha_small_sim.csv")
dd_mc_output_1alpha_med <- read.csv("dd_mc_output_1alpha_med_sim.csv")
dd_mc_output_1alpha_big <- read.csv("dd_mc_output_1alpha_big_sim.csv")
dd_mc_output_1alpha_biggest <- read.csv("dd_mc_output_1alpha_biggest_sim.csv")

pdf("../../figures/raw/ED_Fig_6c.pdf", width=6, height=6) 
ED_Fig_6c <- ggplot() +
  geom_line(data = dd_mc_output_1alpha_small, aes(x=Year, y=DID_diff), color = "dodgerblue1", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_1alpha_med, aes(x=Year, y=DID_diff), color = "dodgerblue2", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_1alpha_big, aes(x=Year, y=DID_diff), color = "dodgerblue3", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_1alpha_biggest, aes(x=Year, y=DID_diff), color = "dodgerblue4", linetype = "dashed", alpha = .8) +
  geom_line(data = dd_mc_output_1alpha_small, aes(x=Year, y=DD_np_diff), color = "dodgerblue1") +
  geom_line(data = dd_mc_output_1alpha_med, aes(x=Year, y=DD_np_diff), color = "dodgerblue2") +
  geom_line(data = dd_mc_output_1alpha_big, aes(x=Year, y=DD_np_diff), color = "dodgerblue3") +
  geom_line(data = dd_mc_output_1alpha_biggest, aes(x=Year, y=DD_np_diff), color = "dodgerblue4") +
  geom_line(data = dd_mc_output_1alpha_small, aes(x=Year, y=MC_diff), color = "orangered1") +
  geom_line(data = dd_mc_output_1alpha_med, aes(x=Year, y=MC_diff), color = "orangered2") +
  geom_line(data = dd_mc_output_1alpha_big, aes(x=Year, y=MC_diff), color = "orangered3") +
  geom_line(data = dd_mc_output_1alpha_biggest, aes(x=Year, y=MC_diff), color = "orangered4") +
  ylim(0,1) +
  ylab("Difference from True ATE") +
  xlab("Years in Timeseries") +
  #theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 18), axis.text = element_text(size = 18), plot.background = element_rect(color = "white"))
ED_Fig_6c
dev.off()




