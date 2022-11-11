library(tidyverse)
library(dplyr)
library(utils)
library(reshape2)
library(fixest)
library(MCPanel)
library(glmnet)
library(ggthemes)
library(ggplot2)
library(infer)


# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tc_2km_ug <- read.csv("tc_2km_ug.csv")
tc_2km_ug_C <- subset(tc_2km_ug, treated == 0)
tc_2km_ug_T <- subset(tc_2km_ug, treated == 1)

input_1 <- read.csv("preds_5p_a_noUG.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_5p_b_noUG.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_5p_c_noUG.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_5p_d_noUG.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_5p_e_noUG.csv")
colnames(input_5)[1] <- "DHSID_year"

control_panel <- cbind(input_1, input_2, input_3, 
                       input_4, input_5)
control_panel <- control_panel[,c(1:2,4,6,8,10)]
control_panel$DHSID <- substr(control_panel$DHSID_year, 1, 14)
control_panel$year <- substr(control_panel$DHSID_year, 16, 19)
control_panel <- subset(control_panel, year < 2017)
control_panel <- subset(control_panel, year > 2005)
control_panel <- merge(tc_2km_ug_C, control_panel, by = "DHSID")

treatment_panel <- cbind(input_1, input_2, input_3, 
                         input_4, input_5)
treatment_panel <- treatment_panel[,c(1:2,4,6,8,10)]
treatment_panel$DHSID <- substr(treatment_panel$DHSID_year, 1, 14)
treatment_panel$year <- substr(treatment_panel$DHSID_year, 16, 19)
treatment_panel <- subset(treatment_panel, year < 2017)
treatment_panel <- subset(treatment_panel, year > 2005)
treatment_panel <- merge(tc_2km_ug_T, treatment_panel, by = "DHSID")

###########
# 2016, 2km, noUG
###########
ptm <- proc.time()

set.seed(456)

n = 500
p = 76 
output_2km_5p_noUG <- matrix(NA, 500, 5) 
inside_boot_run <- matrix(NA,1,p) 

for(m in 1:n) {
  sample_C <- control_panel 
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1)))
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
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
  unique_orig_c <- length(unique(sample_unique_C_units$DHSID))
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") 
  unique_merge_c <- length(unique(sample_control$DHSID))
  if(!unique_orig_c == unique_merge_c){stop("Merge failed, number of units 
                                    in merged data not equal to original data frame")}
  
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  unique_orig_t <- length(unique(sample_unique_T_units$DHSID))
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  unique_merge_t <- length(unique(sample_treatment$DHSID))
  if(!unique_orig_t == unique_merge_t){stop("Merge failed, number of units 
                                    in merged data not equal to original data frame")}
  
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
  
  output_2km_5p_noUG[m,1] <- rowMeans(wide_TC_df_boot[11,889:964])
  
  sample_matrix <- wide_TC_df_boot 
  
  N_basic_control_test <- 964 
  T_basic_control_test <- 11 
  M_basic_control_test <- sample_matrix 
  M_basic_control_test <- as.matrix(M_basic_control_test) 
  
  mask_basic_control_test <- matrix(1,11,964) 
  mask_basic_control_test[6:11,889:964] <- 0 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) 
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  output_2km_5p_noUG[m,2] <- mean(model_with_both_basic_control_test$est[11,889:964])
  
  output_2km_5p_noUG[m,3] <- output_2km_5p_noUG[m,1] - output_2km_5p_noUG[m,2]
  
  X_prod <- t(as.matrix(sample_matrix[1:5,1:888])) 
  Y_prod <- t(as.matrix(sample_matrix[11,1:888])) 
  en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) 
  
  for(o in 1:p) {  
    post_prod <- as.matrix(t(sample_matrix[1:5,888+o])) 
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) 
    
    inside_boot_run[1,o] <- pred_prod 
  }
  
  output_2km_5p_noUG[m,4] <- mean(inside_boot_run[1,])  
  output_2km_5p_noUG[m,5] <- output_2km_5p_noUG[m,1]  - output_2km_5p_noUG[m,4] 
}

colnames(output_2km_5p_noUG) <- c("sample_ave", "mc_ave", "mc_ate", "sc_t_ave", "sc_t_ate")

proc.time() - ptm

# summarizing results
estimates_2km_5pl_noUG <- matrix(NA,1,6)
colnames(estimates_2km_5pl_noUG) <- c( "mc", "mc_l95", "mc_h95",
                                       "sc", "sc_l95", "sc_h95")

estimates_2km_5pl_noUG[1,1] <- round(mean(output_2km_5p_noUG[,3]),3)
rank_mc <- as.data.frame(output_2km_5p_noUG[,3])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE), ])
estimates_2km_5pl_noUG[1,2] <- as.data.frame(get_confidence_interval(rank_mc, level = 0.95, type = "percentile"))[1,1]
estimates_2km_5pl_noUG[1,3] <- as.data.frame(get_confidence_interval(rank_mc, level = 0.95, type = "percentile"))[1,2]

estimates_2km_5pl_noUG[1,4] <-round(mean(output_2km_5p_noUG[,5]),3)
rank_sc <- as.data.frame(output_2km_5p_noUG[,5])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE), ])
estimates_2km_5pl_noUG[1,5] <- as.data.frame(get_confidence_interval(rank_sc, level = 0.95, type = "percentile"))[1,1]
estimates_2km_5pl_noUG[1,6] <- as.data.frame(get_confidence_interval(rank_sc, level = 0.95, type = "percentile"))[1,2]

View(estimates_2km_5pl_noUG)

write.csv(estimates_2km_5pl_noUG, "estimates_2km_5pl_noUG_aug.csv", row.names = FALSE)




