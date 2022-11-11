# This replicates results in Fig 3 for the "full_country" scenario.  

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

randomized_control_r <- read.csv("C_5penalty_raster.csv")
randomized_treatment_r <- read.csv("T_5penalty_raster.csv")

ptm <- proc.time()

set.seed(789) 

n = 500
p = 268
output_2km_5p_raster <- matrix(NA, 500, 5)
inside_boot_run <- matrix(NA,1,p)

for(m in 1:n) {
  sample_C <- randomized_control_r
  randomized_C <- apply(sample_C[,8:12], 1, function(d) sample(d, 1))
  sample_C <- cbind(sample_C, randomized_C)
  colnames(sample_C)[14] <- "randomized"
  
  sample_T <- randomized_treatment_r
  randomized_T <- apply(sample_T[,8:12], 1, function(d) sample(d, 1))
  sample_T<- cbind(sample_T, randomized_T)
  colnames(sample_T)[14] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1])
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 3601, replace = TRUE),]
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  unique_orig_c <- length(unique(sample_unique_C_units$sample_unique_C_units))
  colnames(sample_unique_C_units) <- c("id") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "id")
  unique_merge_c <- length(unique(sample_control$id))
  if(!unique_orig_c == unique_merge_c){stop("Merge failed, number of units 
                                    in merged data not equal to original data frame")}
  
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 268, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  unique_orig_t <- length(unique(sample_unique_T_units$sample_unique_T_units))
  colnames(sample_unique_T_units) <- c("id") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "id")
  unique_merge_t <- length(unique(sample_treatment$id))
  if(!unique_orig_t == unique_merge_t){stop("Merge failed, number of units 
                                    in merged data not equal to original data frame")}
  
  l_c <- length(sample_control[,1])/11
  c_df_to_reshape_boot <- sample_control[,c(1,2,14)] 
  c_df_to_reshape_boot$id <- rep(1:l_c, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,2], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(1,2,3)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  l_t <- length(sample_treatment[,1])/11
  t_df_to_reshape_boot <- sample_treatment[,c(1,2,14)]
  t_df_to_reshape_boot$id <- rep(1:l_t, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,2], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(1,2,3)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  output_2km_5p_raster[m,1] <- rowMeans(wide_TC_df_boot[11,(l_c+1):(l_c+l_t)])
  
  sample_matrix <- wide_TC_df_boot
  
  N_basic_control_test <- l_c+l_t
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,l_c+l_t)
  mask_basic_control_test[6:11,(l_c+1):(l_c+l_t)] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  output_2km_5p_raster[m,2] <- mean(model_with_both_basic_control_test$est[11,(l_c+1):(l_c+l_t)])
  
  output_2km_5p_raster[m,3] <- output_2km_5p_raster[m,1] - output_2km_5p_raster[m,2]
  
  # starting SC section
  X_prod <- t(as.matrix(sample_matrix[1:5,1:l_c]))
  Y_prod <- t(as.matrix(sample_matrix[11,1:l_c]))
  
  en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
  
  for(o in 1:p) {  
    post_prod <- as.matrix(t(sample_matrix[1:5,l_c+o]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min))
    
    inside_boot_run[1,o] <- pred_prod
  }
  
  output_2km_5p_raster[m,4] <- mean(inside_boot_run[1,])  
  output_2km_5p_raster[m,5] <- output_2km_5p_raster[m,1]  - output_2km_5p_raster[m,4]   
}

proc.time() - ptm

colnames(output_2km_5p_raster) <- c("sample_ave", "mc_ave", "mc_ate", "sc_t_ave", "sc_t_ate")

proc.time() - ptm

# summarizing results
estimates_2km_5pl_raster <- matrix(NA,1,6)
colnames(estimates_2km_5pl_raster) <- c("mc", "mc_l95", "mc_h95",
                                        "sc", "sc_l95", "sc_h95")

estimates_2km_5pl_raster[1,1] <- round(mean(output_2km_5p_raster[,3]),3)
rank_mc <- as.data.frame(output_2km_5p_raster[,3])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE), ])
estimates_2km_5pl_raster[1,2] <- as.data.frame(get_confidence_interval(rank_mc, level = 0.95, type = "percentile"))[1,1]
estimates_2km_5pl_raster[1,3] <- as.data.frame(get_confidence_interval(rank_mc, level = 0.95, type = "percentile"))[1,2]

estimates_2km_5pl_raster[1,4] <-round(mean(output_2km_5p_raster[,5]),3)
rank_sc <- as.data.frame(output_2km_5p_raster[,5])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE), ])
estimates_2km_5pl_raster[1,5] <- as.data.frame(get_confidence_interval(rank_sc, level = 0.95, type = "percentile"))[1,1]
estimates_2km_5pl_raster[1,6] <- as.data.frame(get_confidence_interval(rank_sc, level = 0.95, type = "percentile"))[1,2]

View(estimates_2km_5pl_raster)

write.csv(estimates_2km_5pl_raster, "estimates_2km_5pl_raster_n500.csv", row.names = FALSE)


