# This script performs calculations that are in the paper but that are
# not included in the main or ED figures. 

rm(list=ls())

setwd("data/figure_and_input_data/")

# calculates the SD of the wealth index
est_2km_5pl <- read.csv("estimates_2km_5pl_aug.csv")

input_1 <- read.csv("preds_5p_a_f.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_5p_b_f.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_5p_c_f.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_5p_d_f.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_5p_e_f.csv")
colnames(input_5)[1] <- "DHSID_year"

full_5pl_panel <- cbind(input_1, input_2, input_3, 
                        input_4, input_5)
full_5pl_panel <- full_5pl_panel[,c(1:2,4,6,8,10)]
full_5pl_panel$year <- substr(full_5pl_panel$DHSID_year, 16, 19)

full_5pl_panel$meanWI <- rowMeans(full_5pl_panel[,c(2:6)])

# results
est_2km_5pl[1,1] / sd(full_5pl_panel$meanWI) #

###########
# calculates the change in WI using the full country raster 
# and calculates the sd
###########

full_108 <- read.csv("dhs_village_raster_full_108.csv")
full_240 <- read.csv("dhs_village_raster_full_240.csv")
full_372 <- read.csv("dhs_village_raster_full_372.csv")
full_504 <- read.csv("dhs_village_raster_full_504.csv")
full_636 <- read.csv("dhs_village_raster_full_636.csv")

full_108_2006 <- subset(full_108, year == 2006)
full_240_2006 <- subset(full_240, year == 2006)
full_372_2006 <- subset(full_372, year == 2006)
full_504_2006 <- subset(full_504, year == 2006)
full_636_2006 <- subset(full_636, year == 2006)

full_2006 <- rbind(full_108_2006, full_240_2006, full_372_2006,
                   full_504_2006, full_636_2006)
mean(full_2006$raster_val, na.rm = TRUE)
# 

full_108_2016 <- subset(full_108, year == 2016)
full_240_2016 <- subset(full_240, year == 2016)
full_372_2016 <- subset(full_372, year == 2016)
full_504_2016 <- subset(full_504, year == 2016)
full_636_2016 <- subset(full_636, year == 2016)

full_2016 <- rbind(full_108_2016, full_240_2016, full_372_2016,
                   full_504_2016, full_636_2016)
mean(full_2016$raster_val, na.rm = TRUE)
# 

# diff between 2006 and 2016 = 
mean(full_2016$raster_val, na.rm = TRUE) - mean(full_2006$raster_val, na.rm = TRUE)
# 

full_2006_2016 <- rbind(full_2006, full_2016)
sd(full_2006_2016$raster_val, na.rm = TRUE)
# 

# sd change
mean(full_2016$raster_val, na.rm = TRUE) - mean(full_2006$raster_val, na.rm = TRUE) / sd(full_2006_2016$raster_val, na.rm = TRUE)
# 

###########
# r2 calculations in main text
###########

cnn_output_ug <- read.csv("cnn_output_ug.csv")
cnn_output <- read.csv("cnn_output.csv")

# prediction runs
ugpred_1 <- read.csv("CI-V0-1-9_config_id_0.csv")
ugpred_2 <- read.csv("CI-V0-1-9_config_id_132.csv")
ugpred_3 <- read.csv("CI-V0-1-9_config_id_264.csv")
ugpred_4 <- read.csv("CI-V0-1-9_config_id_396.csv")
ugpred_5 <- read.csv("CI-V0-1-9_config_id_528.csv")
ugpred_6 <- read.csv("CI-V0-1-9_config_id_84.csv")
ugpred_7 <- read.csv("CI-V0-1-9_config_id_216.csv")
ugpred_8 <- read.csv("CI-V0-1-9_config_id_348.csv")
ugpred_9 <- read.csv("CI-V0-1-9_config_id_480.csv")
ugpred_10 <- read.csv("CI-V0-1-9_config_id_612.csv")
ugpred_11 <- read.csv("CI-V0-1-9_config_id_96.csv")
ugpred_12 <- read.csv("CI-V0-1-9_config_id_228.csv")
ugpred_13 <- read.csv("CI-V0-1-9_config_id_360.csv")
ugpred_14 <- read.csv("CI-V0-1-9_config_id_492.csv")
ugpred_15 <- read.csv("CI-V0-1-9_config_id_624.csv")
ugpred_16 <- read.csv("CI-V0-1-9_config_id_108.csv")
ugpred_17 <- read.csv("CI-V0-1-9_config_id_240.csv")
ugpred_18 <- read.csv("CI-V0-1-9_config_id_372.csv")
ugpred_19 <- read.csv("CI-V0-1-9_config_id_504.csv")
ugpred_20 <- read.csv("CI-V0-1-9_config_id_636.csv")
ugpred_21 <- read.csv("CI-V0-1-9_config_id_114.csv")
ugpred_22 <- read.csv("CI-V0-1-9_config_id_246.csv")
ugpred_23 <- read.csv("CI-V0-1-9_config_id_378.csv")
ugpred_24 <- read.csv("CI-V0-1-9_config_id_510.csv")
ugpred_25 <- read.csv("CI-V0-1-9_config_id_642.csv")

colnames(ugpred_1) <- c("DHSID", "run")
colnames(ugpred_2) <- c("DHSID", "run")
colnames(ugpred_3) <- c("DHSID", "run")
colnames(ugpred_4) <- c("DHSID", "run")
colnames(ugpred_5) <- c("DHSID", "run")
colnames(ugpred_6) <- c("DHSID", "run")
colnames(ugpred_7) <- c("DHSID", "run")
colnames(ugpred_8) <- c("DHSID", "run")
colnames(ugpred_9) <- c("DHSID", "run")
colnames(ugpred_10) <- c("DHSID", "run")
colnames(ugpred_11) <- c("DHSID", "run")
colnames(ugpred_12) <- c("DHSID", "run")
colnames(ugpred_13) <- c("DHSID", "run")
colnames(ugpred_14) <- c("DHSID", "run")
colnames(ugpred_15) <- c("DHSID", "run")
colnames(ugpred_16) <- c("DHSID", "run")
colnames(ugpred_17) <- c("DHSID", "run")
colnames(ugpred_18) <- c("DHSID", "run")
colnames(ugpred_19) <- c("DHSID", "run")
colnames(ugpred_20) <- c("DHSID", "run")
colnames(ugpred_21) <- c("DHSID", "run")
colnames(ugpred_22) <- c("DHSID", "run")
colnames(ugpred_23) <- c("DHSID", "run")
colnames(ugpred_24) <- c("DHSID", "run")
colnames(ugpred_25) <- c("DHSID", "run")

ugpred_1_m <- merge(cnn_output, ugpred_1, by = "DHSID" )
ugpred_2_m <- merge(cnn_output, ugpred_2, by = "DHSID" )
ugpred_3_m <- merge(cnn_output, ugpred_3, by = "DHSID" )
ugpred_4_m <- merge(cnn_output, ugpred_4, by = "DHSID" )
ugpred_5_m <- merge(cnn_output, ugpred_5, by = "DHSID" )
ugpred_6_m <- merge(cnn_output, ugpred_6, by = "DHSID" )
ugpred_7_m <- merge(cnn_output, ugpred_7, by = "DHSID" )
ugpred_8_m <- merge(cnn_output, ugpred_8, by = "DHSID" )
ugpred_9_m <- merge(cnn_output, ugpred_9, by = "DHSID" )
ugpred_10_m <- merge(cnn_output, ugpred_10, by = "DHSID" )
ugpred_11_m <- merge(cnn_output, ugpred_11, by = "DHSID" )
ugpred_12_m <- merge(cnn_output, ugpred_12, by = "DHSID" )
ugpred_13_m <- merge(cnn_output, ugpred_13, by = "DHSID" )
ugpred_14_m <- merge(cnn_output, ugpred_14, by = "DHSID" )
ugpred_15_m <- merge(cnn_output, ugpred_15, by = "DHSID" )
ugpred_16_m <- merge(cnn_output, ugpred_16, by = "DHSID" )
ugpred_17_m <- merge(cnn_output, ugpred_17, by = "DHSID" )
ugpred_18_m <- merge(cnn_output, ugpred_18, by = "DHSID" )
ugpred_19_m <- merge(cnn_output, ugpred_19, by = "DHSID" )
ugpred_20_m <- merge(cnn_output, ugpred_20, by = "DHSID" )
ugpred_21_m <- merge(cnn_output, ugpred_21, by = "DHSID" )
ugpred_22_m <- merge(cnn_output, ugpred_22, by = "DHSID" )
ugpred_23_m <- merge(cnn_output, ugpred_23, by = "DHSID" )
ugpred_24_m <- merge(cnn_output, ugpred_24, by = "DHSID" )
ugpred_25_m <- merge(cnn_output, ugpred_25, by = "DHSID" )

summary(lm(run ~ pred, data=ugpred_1_m))
summary(lm(run ~ pred, data=ugpred_2_m))
summary(lm(run ~ pred, data=ugpred_3_m))
summary(lm(run ~ pred, data=ugpred_4_m))
summary(lm(run ~ pred, data=ugpred_5_m))
summary(lm(run ~ pred, data=ugpred_6_m))
summary(lm(run ~ pred, data=ugpred_7_m))
summary(lm(run ~ pred, data=ugpred_8_m))
summary(lm(run ~ pred, data=ugpred_9_m))
summary(lm(run ~ pred, data=ugpred_10_m))
summary(lm(run ~ pred, data=ugpred_11_m))
summary(lm(run ~ pred, data=ugpred_12_m))
summary(lm(run ~ pred, data=ugpred_13_m))
summary(lm(run ~ pred, data=ugpred_14_m))
summary(lm(run ~ pred, data=ugpred_15_m))
summary(lm(run ~ pred, data=ugpred_16_m))
summary(lm(run ~ pred, data=ugpred_17_m))
summary(lm(run ~ pred, data=ugpred_18_m))
summary(lm(run ~ pred, data=ugpred_19_m))
summary(lm(run ~ pred, data=ugpred_20_m))
summary(lm(run ~ pred, data=ugpred_21_m))
summary(lm(run ~ pred, data=ugpred_22_m))
summary(lm(run ~ pred, data=ugpred_23_m))
summary(lm(run ~ pred, data=ugpred_24_m))
summary(lm(run ~ pred, data=ugpred_25_m))

###########
# calculating the r2 for the noUG predictions
###########

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

ug_pred_panel_noUG <- cbind(input_1, input_2, input_3, 
                            input_4, input_5)
ug_pred_panel_noUG <- ug_pred_panel_noUG[,c(1:2,4,6,8,10)]
ug_pred_panel_noUG$DHSID <- substr(ug_pred_panel_noUG$DHSID_year, 1, 14)
ug_pred_panel_noUG$year <- substr(ug_pred_panel_noUG$DHSID_year, 16, 19)
ug_pred_panel_noUG$survey_year <- substr(ug_pred_panel_noUG$DHSID_year, 3, 6)
ug_pred_panel_noUG <- subset(ug_pred_panel_noUG, year == survey_year)
ug_pred_panel_noUG$meanWI = rowMeans(ug_pred_panel_noUG[,2:6])

ug_pred_panel_m_noUG <- merge(cnn_output, ug_pred_panel_noUG, by = "DHSID" )
colnames(ug_pred_panel_m_noUG)[9:13] <- c("run1", "run2", "run3", "run4", "run5")
summary(lm(meanWI ~ pred, data=ug_pred_panel_m_noUG))


###########
# This following script compares predictive accuracy for treated units in 2010, the last pre-treatment year. 
###########

rm(list=ls())

library(dplyr)
library(reshape2)
library(fixest)
library(MCPanel)
library(glmnet)
library(ggplot2)
library(ggthemes)
library(infer)
library(utils)

setwd("data/figure_and_input_data/")

tc_2km_ug <- read.csv("tc_2km_ug.csv") # this is the treatment and control group with a 2km buffer
tc_2km_ug_C <- subset(tc_2km_ug, treated == 0)
tc_2km_ug_T <- subset(tc_2km_ug, treated == 1)

# Loading the wealth index predictions

input_1 <- read.csv("preds_5p_a_f.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_5p_b_f.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_5p_c_f.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_5p_d_f.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_5p_e_f.csv")
colnames(input_5)[1] <- "DHSID_year"

control_panel <- cbind(input_1, input_2, input_3, 
                       input_4, input_5)
control_panel <- control_panel[,c(1:2,4,6,8,10)]
control_panel$DHSID <- substr(control_panel$DHSID_year, 1, 14)
control_panel$year <- substr(control_panel$DHSID_year, 16, 19)

control_panel <- merge(tc_2km_ug_C, control_panel, by = "DHSID")

treatment_panel <- cbind(input_1, input_2, input_3, 
                         input_4, input_5)
treatment_panel <- treatment_panel[,c(1:2,4,6,8,10)]
treatment_panel$DHSID <- substr(treatment_panel$DHSID_year, 1, 14)
treatment_panel$year <- substr(treatment_panel$DHSID_year, 16, 19)

treatment_panel <- merge(tc_2km_ug_T, treatment_panel, by = "DHSID")

# 2010, 2km loop, pretrend CV
ptm <- proc.time()

set.seed(456)

n = 100 # Number of loops
p = 76 # Number of SC-EN iterations, equals the number of treated units.  
output_2km_5p_2010 <- matrix(NA, n, 6) # This df collects observations in each loop.
inside_boot_run <- matrix(NA,1,p) # This df collects sc-en estimates per loop. 

for(m in 1:n) {
  sample_C <- control_panel 
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # Randomly select 1 of 5 WI predictions per unit. 
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # Add randomized selection to the control df
  colnames(sample_C)[13] <- "randomized"
  
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1))) # Randomly select 1 of 5 WI predictions per unit.
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T)) # Add randomized selection to the control df
  colnames(sample_T)[13] <- "randomized"
  
  # Here i randomly select control units, with replacement, in each loop. 
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
  
  # Here i randomly select treated units, with replacement, in each loop. 
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
  
  dd_sample_TC <- rbind(sample_control, sample_treatment)
  dd_sample_TC <- subset(dd_sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  dd_sample_TC$fake_treat_year <- sapply(dd_sample_TC$year, fake_treated)
  dd_sample_TC$fake_treat_id <- dd_sample_TC$treated * dd_sample_TC$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  pre_dd <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=dd_sample_TC))
  sample_coeff <- as.data.frame(coeftable(pre_dd)[1])
  output_2km_5p_2010[m,1] <-sample_coeff[1,1]
  
  
  # Here I am reshaping the control group from long to wide. 
  c_df_to_reshape_boot <- sample_control[,c(1,12,13)]
  c_df_to_reshape_boot$id <- rep(1:888, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  # Here I am reshaping the treatment group from long to wide. 
  t_df_to_reshape_boot <- sample_treatment[,c(1,12,13)]
  t_df_to_reshape_boot$id <- rep(1:76, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  wide_TC_df_boot <- wide_TC_df_boot[1:5,]
  
  output_2km_5p_2010[m,2] <- rowMeans(wide_TC_df_boot[5,889:964]) # Collects the observed average of treated units for each loop. 
  
  sample_matrix <- wide_TC_df_boot 
  
  # Here I run the matrix completion estimates.  
  N_basic_control_test <- 964 
  T_basic_control_test <- 5 
  M_basic_control_test <- sample_matrix 
  M_basic_control_test <- as.matrix(M_basic_control_test) 
  
  mask_basic_control_test <- matrix(1,5,964) 
  mask_basic_control_test[5,889:964] <- 0 # Masking observed outcomes with zeros. 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) 
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  output_2km_5p_2010[m,3] <- mean(model_with_both_basic_control_test$est[5,889:964])
  
  output_2km_5p_2010[m,4] <- output_2km_5p_2010[m,2] - output_2km_5p_2010[m,3]
  
  # Here I run the matrix completion estimates. 
  X_prod <- t(as.matrix(sample_matrix[1:4,1:888]))
  Y_prod <- t(as.matrix(sample_matrix[5,1:888]))
  
  en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
  
  for(o in 1:p) {  
    post_prod <- as.matrix(t(sample_matrix[1:4,888+o])) # This applies the weights predicted above to each unit. 
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) 
    
    inside_boot_run[1,o] <- pred_prod 
  }
  
  output_2km_5p_2010[m,5] <- mean(inside_boot_run[1,])  
  output_2km_5p_2010[m,6] <- output_2km_5p_2010[m,2]  - output_2km_5p_2010[m,5] 
}

colnames(output_2km_5p_2010) <- c("dd", "obs_ave", "mc_ave", "mc_ate", "sc_ave", "sc_ate")

proc.time() - ptm

mean(output_2km_5p_2010[,1])
mean(output_2km_5p_2010[,4])
mean(output_2km_5p_2010[,6])

###########


