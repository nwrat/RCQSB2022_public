# This code creates Extended Data Figure 07.  

rm(list=ls())

library(ggplot2)
library(ggthemes)

setwd("data/figure_and_input_data/")

cv_fold_1 <- read.csv("cv_fold_1.csv")
cv_fold_2 <- read.csv("cv_fold_2.csv")
cv_fold_3 <- read.csv("cv_fold_3.csv")
cv_fold_4 <- read.csv("cv_fold_4.csv")
cv_fold_5 <- read.csv("cv_fold_5.csv")
cv_fold_6 <- read.csv("cv_fold_6.csv")
cv_fold_7 <- read.csv("cv_fold_7.csv")
cv_fold_8 <- read.csv("cv_fold_8.csv")
cv_fold_9 <- read.csv("cv_fold_9.csv")
cv_fold_10 <- read.csv("cv_fold_10.csv")

# This code creates each of the 10 sub figures

pdf("../../figures/raw/EDFigure07_1cv.pdf", width=8, height=8) 
EDFigure07_1cv <- ggplot(data = cv_fold_1, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("Wealth Index") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), 
        axis.text = element_text(size = 30), plot.background = element_rect(color = "white"), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 60, hjust=1, color = "white")) 
EDFigure07_1cv
dev.off()

pdf("../../figures/raw/EDFigure07_2cv.pdf", width=8, height=8) 
EDFigure07_2cv <- ggplot(data = cv_fold_2, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.x = element_blank(),       
    axis.text.x = element_text(angle = 60, hjust=1, color = "white"),
    axis.ticks.y =  element_blank(),      
    axis.text.y = element_text(color = "white")) 
EDFigure07_2cv
dev.off()

pdf("../../figures/raw/EDFigure07_3cv.pdf", width=8, height=8) 
EDFigure07_3cv <- ggplot(data = cv_fold_3, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.x = element_blank(),       
    axis.text.x = element_text(angle = 60, hjust=1, color = "white"),
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white")) 
EDFigure07_3cv
dev.off()

pdf("../../figures/raw/EDFigure07_4cv.pdf", width=8, height=8) 
EDFigure07_4cv <- ggplot(data = cv_fold_4, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.x = element_blank(),       
    axis.text.x = element_text(angle = 60, hjust=1, color = "white"),
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white"))
EDFigure07_4cv
dev.off()

pdf("../../figures/raw/EDFigure07_5cv.pdf", width=8, height=8) 
EDFigure07_5cv <- ggplot(data = cv_fold_5, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.x = element_blank(),       
    axis.text.x = element_text(angle = 60, hjust=1, color = "white"),
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white"))
EDFigure07_5cv
dev.off()

pdf("../../figures/raw/EDFigure07_6cv.pdf", width=8, height=8) 
EDFigure07_6cv <- ggplot(data = cv_fold_6, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("Wealth Index") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.text.x = element_text(angle = 60, hjust=1))
EDFigure07_6cv
dev.off()

pdf("../../figures/raw/EDFigure07_7cv.pdf", width=8, height=8) 
EDFigure07_7cv <- ggplot(data = cv_fold_7, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white"), 
    axis.text.x = element_text(angle = 60, hjust=1)) 
EDFigure07_7cv
dev.off()

pdf("../../figures/raw/EDFigure07_8cv.pdf", width=8, height=8) 
EDFigure07_8cv <- ggplot(data = cv_fold_8, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white"),
    axis.text.x = element_text(angle = 60, hjust=1)) 
EDFigure07_8cv
dev.off()

pdf("../../figures/raw/EDFigure07_9cv.pdf", width=8, height=8) 
EDFigure07_9cv <- ggplot(data = cv_fold_9, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), plot.background = element_rect(
    color = "white"), 
    axis.ticks.y = element_blank(),       
    axis.text.y = element_text(color = "white"), 
    axis.text.x = element_text(angle = 60, hjust=1)) 
EDFigure07_9cv
dev.off()

pdf("../../figures/raw/EDFigure07_10cv.pdf", width=8, height=8) 
EDFigure07_10cv <- ggplot(data = cv_fold_10, aes(x = V1)) + 
  geom_line(aes(y = V2), color = "tomato") +
  geom_line(aes(y = V3), color = "springgreen") +
  geom_line(aes(y = V4), color = "dodgerblue") +
  geom_line(aes(y = V6), color = "black") +
  ylim(-2,-1) +
  ylab("") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", axis.title = element_text(size = 30), axis.text = element_text(size = 30), 
        plot.background = element_rect(color = "white"), 
        axis.ticks.y = element_blank(),       
        axis.text.y = element_text(color = "white"),
        axis.text.x = element_text(angle = 60, hjust=1)) 
EDFigure07_10cv
dev.off()





