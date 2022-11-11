# This script runs the analysis for Extended Data Figure 08. 
# The code mimics that from Figure 03 data prep. 

rm(list=ls())

library(tidyverse)
library(dplyr)
library(utils)
library(reshape2)
library(fixest)
library(MCPanel)
library(glmnet)
library(ggthemes)
library(ggplot2)

setwd("data/figure_and_input_data/")

tc_3km_ug <- read.csv("tc_3km_ug.csv")
tc_3km_ug_C <- subset(tc_3km_ug, treated == 0)
tc_3km_ug_T <- subset(tc_3km_ug, treated == 1)

tc_4km_ug <- read.csv("tc_4km_ug.csv")
tc_4km_ug_C <- subset(tc_4km_ug, treated == 0)
tc_4km_ug_T <- subset(tc_4km_ug, treated == 1)

###########
#  3km with UG
###########

input_1 <- read.csv("preds_5p_a.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_5p_b.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_5p_c.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_5p_d.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_5p_e.csv")
colnames(input_5)[1] <- "DHSID_year"

control_panel <- cbind(input_1, input_2, input_3, 
                       input_4, input_5)
control_panel <- control_panel[,c(1:2,4,6,8,10)]
control_panel$DHSID <- substr(control_panel$DHSID_year, 1, 14)
control_panel$year <- substr(control_panel$DHSID_year, 16, 19)
control_panel <- subset(control_panel, year < 2017)
control_panel <- subset(control_panel, year > 2005)
control_panel <- merge(tc_3km_ug_C, control_panel, by = "DHSID")

treatment_panel <- cbind(input_1, input_2, input_3, 
                         input_4, input_5)
treatment_panel <- treatment_panel[,c(1:2,4,6,8,10)]
treatment_panel$DHSID <- substr(treatment_panel$DHSID_year, 1, 14)
treatment_panel$year <- substr(treatment_panel$DHSID_year, 16, 19)
treatment_panel <- subset(treatment_panel, year < 2017)
treatment_panel <- subset(treatment_panel, year > 2005)
treatment_panel <- merge(tc_3km_ug_T, treatment_panel, by = "DHSID")

ptm <- proc.time()

set.seed(789) 

n = 100 
p = 76
output_3km_5p_2016 <- matrix(NA, 100, 5)
inside_boot_run <- matrix(NA,1,p) #  

for(m in 1:n) {
  sample_C <- control_panel 
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1)))  
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) 
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) 
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") 
  
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  c_df_to_reshape_boot <- sample_control[,c(1,12,13)]
  c_df_to_reshape_boot$id <- rep(1:888, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  t_df_to_reshape_boot <- sample_treatment[,c(1,12,13)]
  t_df_to_reshape_boot$id <- rep(1:76, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  output_3km_5p_2016[m,1] <- rowMeans(wide_TC_df_boot[11,889:964])
  
  sample_matrix <- wide_TC_df_boot  
  
  N_basic_control_test <- 964 # = N
  T_basic_control_test <- 11 # = T
  M_basic_control_test <- sample_matrix 
  M_basic_control_test <- as.matrix(M_basic_control_test) 
  
  mask_basic_control_test <- matrix(1,11,964) 
  mask_basic_control_test[6:11,889:964] <- 0 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) 
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  output_3km_5p_2016[m,2] <- mean(model_with_both_basic_control_test$est[11,889:964])
  
  output_3km_5p_2016[m,3] <- output_3km_5p_2016[m,1] - output_3km_5p_2016[m,2]
  
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(sample_matrix[1:5,1:888])) 
    Y_prod <- t(as.matrix(sample_matrix[11,1:888])) 
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) 
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min)) 
    
    post_prod <- as.matrix(t(sample_matrix[1:5,888+o]))  
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) 
    
    inside_boot_run[1,o] <- pred_prod  
  }
  
  output_3km_5p_2016[m,4] <- mean(inside_boot_run[1,])  
  output_3km_5p_2016[m,5] <- output_3km_5p_2016[m,1]  - output_3km_5p_2016[m,4] 
}


colnames(output_3km_5p_2016) <- c("sample_ave", "mc_ave", "mc_ate", "sc_t_ave", "sc_t_ate")

proc.time() - ptm   

estimates_3km_5p <- matrix(NA,1,6)
colnames(estimates_3km_5p) <- c( "mc", "mc_l95", "mc_h95",
                                  "sc", "sc_l95", "sc_h95")

estimates_3km_5p[1,1] <- round(mean(output_3km_5p_2016[,3]),3)
rank_mc <- as.data.frame(output_3km_5p_2016[,3])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE), ])
estimates_3km_5p[1,2] <- round(mean(rank_mc[2:3,]), 3)
estimates_3km_5p[1,3] <- round(mean(rank_mc[98:99,]), 3)

estimates_3km_5p[1,4] <-round(mean(output_3km_5p_2016[,5]),3)
rank_sc <- as.data.frame(output_3km_5p_2016[,5])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE), ])
estimates_3km_5p[1,5] <- round(mean(rank_sc[2:3,]), 3)
estimates_3km_5p[1,6] <- round(mean(rank_sc[98:99,]), 3)

View(estimates_3km_5p)

write.csv(estimates_3km_5p, "estimates_3km_5p.csv", row.names = FALSE)

###########
# with UG 4km results
###########

control_panel <- cbind(input_1, input_2, input_3, 
                       input_4, input_5)
control_panel <- control_panel[,c(1:2,4,6,8,10)]
control_panel$DHSID <- substr(control_panel$DHSID_year, 1, 14)
control_panel$year <- substr(control_panel$DHSID_year, 16, 19)
control_panel <- subset(control_panel, year < 2017)
control_panel <- subset(control_panel, year > 2005)
control_panel <- merge(tc_4km_ug_C, control_panel, by = "DHSID")

treatment_panel <- cbind(input_1, input_2, input_3, 
                         input_4, input_5)
treatment_panel <- treatment_panel[,c(1:2,4,6,8,10)]
treatment_panel$DHSID <- substr(treatment_panel$DHSID_year, 1, 14)
treatment_panel$year <- substr(treatment_panel$DHSID_year, 16, 19)
treatment_panel <- subset(treatment_panel, year < 2017)
treatment_panel <- subset(treatment_panel, year > 2005)
treatment_panel <- merge(tc_4km_ug_T, treatment_panel, by = "DHSID")

ptm <- proc.time()

set.seed(789) 

n = 100 
p = 76
output_4km_5p_2016 <- matrix(NA, 100, 5) 
inside_boot_run <- matrix(NA,1,p)  

for(m in 1:n) {
  sample_C <- control_panel 
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1)))  
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) 
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) 
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),]  
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") 
  
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  c_df_to_reshape_boot <- sample_control[,c(1,12,13)]
  c_df_to_reshape_boot$id <- rep(1:888, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  t_df_to_reshape_boot <- sample_treatment[,c(1,12,13)]
  t_df_to_reshape_boot$id <- rep(1:76, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  output_4km_5p_2016[m,1] <- rowMeans(wide_TC_df_boot[11,889:964])
  
  sample_matrix <- wide_TC_df_boot 
  
  N_basic_control_test <- 964 # = N
  T_basic_control_test <- 11 # = T
  M_basic_control_test <- sample_matrix 
  M_basic_control_test <- as.matrix(M_basic_control_test) 
  
  mask_basic_control_test <- matrix(1,11,964) 
  mask_basic_control_test[6:11,889:964] <- 0 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) 
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  output_4km_5p_2016[m,2] <- mean(model_with_both_basic_control_test$est[11,889:964])
  
  output_4km_5p_2016[m,3] <- output_4km_5p_2016[m,1] - output_4km_5p_2016[m,2]
  
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(sample_matrix[1:5,1:888])) 
    Y_prod <- t(as.matrix(sample_matrix[11,1:888])) 
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) 
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min)) 
    
    post_prod <- as.matrix(t(sample_matrix[1:5,888+o]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) 
    
    inside_boot_run[1,o] <- pred_prod 
  }
  
  output_4km_5p_2016[m,4] <- mean(inside_boot_run[1,]) 
  output_4km_5p_2016[m,5] <- output_4km_5p_2016[m,1]  - output_4km_5p_2016[m,4] 
}


colnames(output_4km_5p_2016) <- c("sample_ave", "mc_ave", "mc_ate", "sc_t_ave", "sc_t_ate")

proc.time() - ptm   

estimates_4km_5p <- matrix(NA,1,6)
colnames(estimates_4km_5p) <- c( "mc", "mc_l95", "mc_h95",
                                 "sc", "sc_l95", "sc_h95")

estimates_4km_5p[1,1] <- round(mean(output_4km_5p_2016[,3]),3)
rank_mc <- as.data.frame(output_4km_5p_2016[,3])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE), ])
estimates_4km_5p[1,2] <- round(mean(rank_mc[2:3,]), 3)
estimates_4km_5p[1,3] <- round(mean(rank_mc[98:99,]), 3)

estimates_4km_5p[1,4] <-round(mean(output_4km_5p_2016[,5]),3)
rank_sc <- as.data.frame(output_4km_5p_2016[,5])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE), ])
estimates_4km_5p[1,5] <- round(mean(rank_sc[2:3,]), 3)
estimates_4km_5p[1,6] <- round(mean(rank_sc[98:99,]), 3)

View(estimates_4km_5p)

write.csv(estimates_4km_5p, "estimates_4km_5p.csv", row.names = FALSE)
