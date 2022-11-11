# This code runs cross validation and prepares DFs that creates Extended Data Figure 07.  

rm(list=ls())

library(dplyr)
library(reshape2)
library(MCPanel) 
library(glmnet)
library(fixest)


setwd("data/figure_and_input_data/")

wide_control_df <- read.csv("wide_control_df.csv")

# In each of 10 loops the 10% hold out set is predicted, and the difference between observed and predicted values are cataloged.  

# CV for fold 1 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_1 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

for(s in 1:t) {
  CV_output_1[1:6,] <- wide_control_df[6:11,801:888] 
  # MC
  sample_matrix <- wide_control_df
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_1[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_1[13:18,1:88] <- CV_output_1[1:6,1:88] - CV_output_1[7:12,1:88]
  CV_output_1[19:24,1:88] <- (CV_output_1[13:18,1:88])^2
  
  # SC-ENt
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = min(en_prod$lambda.min)) 
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_1[25:30,1:88] <- sc_CV 
  CV_output_1[31:36,1:88] <- CV_output_1[1:6,1:88]  - CV_output_1[25:30,1:88]   
  CV_output_1[37:42,1:88] <- (CV_output_1[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),]) # selects each unit individually
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_1[43,z] <-sample_coeff[1,1]
  }
  CV_output_1[49:54,1:88] <- CV_output_1[1:6,1:88]  - CV_output_1[43:48,1:88]   
  CV_output_1[55:60,1:88] <- (CV_output_1[43:48,1:88])^2
}

proc.time() - ptm

# this creates a df of mean values
# in order, it is the observed value, MC, SC-EN, DD
mean_df <- as.data.frame(matrix(NA,6,10))
mean_df[,1] <- c(2011:2016)
mean_df[1,2] <- rowMeans(CV_output_1[1,])
mean_df[2,2] <- rowMeans(CV_output_1[2,])
mean_df[3,2] <- rowMeans(CV_output_1[3,])
mean_df[4,2] <- rowMeans(CV_output_1[4,])
mean_df[5,2] <- rowMeans(CV_output_1[5,])
mean_df[6,2] <- rowMeans(CV_output_1[6,])

mean_df[1,3] <- rowMeans(CV_output_1[7,])
mean_df[2,3] <- rowMeans(CV_output_1[8,])
mean_df[3,3] <- rowMeans(CV_output_1[9,])
mean_df[4,3] <- rowMeans(CV_output_1[10,])
mean_df[5,3] <- rowMeans(CV_output_1[11,])
mean_df[6,3] <- rowMeans(CV_output_1[12,])

mean_df[1,4] <- rowMeans(CV_output_1[25,])
mean_df[2,4] <- rowMeans(CV_output_1[26,])
mean_df[3,4] <- rowMeans(CV_output_1[27,])
mean_df[4,4] <- rowMeans(CV_output_1[28,])
mean_df[5,4] <- rowMeans(CV_output_1[29,])
mean_df[6,4] <- rowMeans(CV_output_1[30,])

mean_df[1,5] <- rowMeans(CV_output_1[43,])
mean_df[2,5] <- rowMeans(CV_output_1[44,])
mean_df[3,5] <- rowMeans(CV_output_1[45,])
mean_df[4,5] <- rowMeans(CV_output_1[46,])
mean_df[5,5] <- rowMeans(CV_output_1[47,])
mean_df[6,5] <- rowMeans(CV_output_1[48,])

#this subtracts the did from the means
mean_df[1,6] <- mean_df[1,2] - mean_df[1,5]
mean_df[2,6] <- mean_df[2,2] - mean_df[2,5]
mean_df[3,6] <- mean_df[3,2] - mean_df[3,5]
mean_df[4,6] <- mean_df[4,2] - mean_df[4,5]
mean_df[5,6] <- mean_df[5,2] - mean_df[5,5]
mean_df[6,6] <- mean_df[6,2] - mean_df[6,5]

#now we are calculating the mean error
#mean diff of mc
mean_df[1,7] <- (mean_df[1,3] - mean_df[1,2])
mean_df[2,7] <- (mean_df[2,3] - mean_df[2,2])
mean_df[3,7] <- (mean_df[3,3] - mean_df[3,2])
mean_df[4,7] <- (mean_df[4,3] - mean_df[4,2])
mean_df[5,7] <- (mean_df[5,3] - mean_df[5,2])
mean_df[6,7] <- (mean_df[6,3] - mean_df[6,2])

#mean diff of sc
mean_df[1,8] <- (mean_df[1,4] - mean_df[1,2])
mean_df[2,8] <- (mean_df[2,4] - mean_df[2,2])
mean_df[3,8] <- (mean_df[3,4] - mean_df[3,2])
mean_df[4,8] <- (mean_df[4,4] - mean_df[4,2])
mean_df[5,8] <- (mean_df[5,4] - mean_df[5,2])
mean_df[6,8] <- (mean_df[6,4] - mean_df[6,2])

#mean diff of did
mean_df[1,9] <- (mean_df[1,6] - mean_df[1,2])
mean_df[2,9] <- (mean_df[2,6] - mean_df[2,2])
mean_df[3,9] <- (mean_df[3,6] - mean_df[3,2])
mean_df[4,9] <- (mean_df[4,6] - mean_df[4,2])
mean_df[5,9] <- (mean_df[5,6] - mean_df[5,2])
mean_df[6,9] <- (mean_df[6,6] - mean_df[6,2])


# Fold 2 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_2 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_2 <- wide_control_df[,c(89:888,1:88)]

for(s in 1:t) {
  CV_output_2[1:6,] <- wide_control_df_2[6:11,801:888] #this will change in future iterations. 
  
  #MC
  sample_matrix <- wide_control_df_2
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_2[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_2[13:18,1:88] <- CV_output_2[1:6,1:88] - CV_output_2[7:12,1:88]
  CV_output_2[19:24,1:88] <- (CV_output_2[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_2[25:30,1:88] <- sc_CV 
  CV_output_2[31:36,1:88] <- CV_output_2[1:6,1:88]  - CV_output_2[25:30,1:88]   
  CV_output_2[37:42,1:88] <- (CV_output_2[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_2
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_2[43,z] <-sample_coeff[1,1]
  }
  CV_output_2[49:54,1:88] <- CV_output_2[1:6,1:88]  - CV_output_2[43:48,1:88]   
  CV_output_2[55:60,1:88] <- (CV_output_2[43:48,1:88])^2
}

proc.time() - ptm

# 
mean_df_2 <- as.data.frame(matrix(NA,6,5))
mean_df_2[,1] <- c(2011:2016)
mean_df_2[1,2] <- rowMeans(CV_output_2[1,])
mean_df_2[2,2] <- rowMeans(CV_output_2[2,])
mean_df_2[3,2] <- rowMeans(CV_output_2[3,])
mean_df_2[4,2] <- rowMeans(CV_output_2[4,])
mean_df_2[5,2] <- rowMeans(CV_output_2[5,])
mean_df_2[6,2] <- rowMeans(CV_output_2[6,])

mean_df_2[1,3] <- rowMeans(CV_output_2[7,])
mean_df_2[2,3] <- rowMeans(CV_output_2[8,])
mean_df_2[3,3] <- rowMeans(CV_output_2[9,])
mean_df_2[4,3] <- rowMeans(CV_output_2[10,])
mean_df_2[5,3] <- rowMeans(CV_output_2[11,])
mean_df_2[6,3] <- rowMeans(CV_output_2[12,])

mean_df_2[1,4] <- rowMeans(CV_output_2[25,])
mean_df_2[2,4] <- rowMeans(CV_output_2[26,])
mean_df_2[3,4] <- rowMeans(CV_output_2[27,])
mean_df_2[4,4] <- rowMeans(CV_output_2[28,])
mean_df_2[5,4] <- rowMeans(CV_output_2[29,])
mean_df_2[6,4] <- rowMeans(CV_output_2[30,])

mean_df_2[1,5] <- rowMeans(CV_output_2[43,])
mean_df_2[2,5] <- rowMeans(CV_output_2[44,])
mean_df_2[3,5] <- rowMeans(CV_output_2[45,])
mean_df_2[4,5] <- rowMeans(CV_output_2[46,])
mean_df_2[5,5] <- rowMeans(CV_output_2[47,])
mean_df_2[6,5] <- rowMeans(CV_output_2[48,])

#this subtracts the did from the means
mean_df_2[1,6] <- mean_df_2[1,2] - mean_df_2[1,5]
mean_df_2[2,6] <- mean_df_2[2,2] - mean_df_2[2,5]
mean_df_2[3,6] <- mean_df_2[3,2] - mean_df_2[3,5]
mean_df_2[4,6] <- mean_df_2[4,2] - mean_df_2[4,5]
mean_df_2[5,6] <- mean_df_2[5,2] - mean_df_2[5,5]
mean_df_2[6,6] <- mean_df_2[6,2] - mean_df_2[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_2[1,7] <- (mean_df_2[1,3] - mean_df_2[1,2])
mean_df_2[2,7] <- (mean_df_2[2,3] - mean_df_2[2,2])
mean_df_2[3,7] <- (mean_df_2[3,3] - mean_df_2[3,2])
mean_df_2[4,7] <- (mean_df_2[4,3] - mean_df_2[4,2])
mean_df_2[5,7] <- (mean_df_2[5,3] - mean_df_2[5,2])
mean_df_2[6,7] <- (mean_df_2[6,3] - mean_df_2[6,2])

#mean diff of sc
mean_df_2[1,8] <- (mean_df_2[1,4] - mean_df_2[1,2])
mean_df_2[2,8] <- (mean_df_2[2,4] - mean_df_2[2,2])
mean_df_2[3,8] <- (mean_df_2[3,4] - mean_df_2[3,2])
mean_df_2[4,8] <- (mean_df_2[4,4] - mean_df_2[4,2])
mean_df_2[5,8] <- (mean_df_2[5,4] - mean_df_2[5,2])
mean_df_2[6,8] <- (mean_df_2[6,4] - mean_df_2[6,2])

#mean diff of did
mean_df_2[1,9] <- (mean_df_2[1,6] - mean_df_2[1,2])
mean_df_2[2,9] <- (mean_df_2[2,6] - mean_df_2[2,2])
mean_df_2[3,9] <- (mean_df_2[3,6] - mean_df_2[3,2])
mean_df_2[4,9] <- (mean_df_2[4,6] - mean_df_2[4,2])
mean_df_2[5,9] <- (mean_df_2[5,6] - mean_df_2[5,2])
mean_df_2[6,9] <- (mean_df_2[6,6] - mean_df_2[6,2])


# Fold 3 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_3 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_3 <- wide_control_df[,c(1:88,177:888,89:176)]

for(s in 1:t) {
  CV_output_3[1:6,] <- wide_control_df_3[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_3
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_3[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_3[13:18,1:88] <- CV_output_3[1:6,1:88] - CV_output_3[7:12,1:88]
  CV_output_3[19:24,1:88] <- (CV_output_3[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_3[25:30,1:88] <- sc_CV 
  CV_output_3[31:36,1:88] <- CV_output_3[1:6,1:88]  - CV_output_3[25:30,1:88]   
  CV_output_3[37:42,1:88] <- (CV_output_3[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_3
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_3[43,z] <-sample_coeff[1,1]
  }
  CV_output_3[49:54,1:88] <- CV_output_3[1:6,1:88]  - CV_output_3[43:48,1:88]   
  CV_output_3[55:60,1:88] <- (CV_output_3[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_3 <- as.data.frame(matrix(NA,6,5))
mean_df_3[,1] <- c(2011:2016)
mean_df_3[1,2] <- rowMeans(CV_output_3[1,])
mean_df_3[2,2] <- rowMeans(CV_output_3[2,])
mean_df_3[3,2] <- rowMeans(CV_output_3[3,])
mean_df_3[4,2] <- rowMeans(CV_output_3[4,])
mean_df_3[5,2] <- rowMeans(CV_output_3[5,])
mean_df_3[6,2] <- rowMeans(CV_output_3[6,])

mean_df_3[1,3] <- rowMeans(CV_output_3[7,])
mean_df_3[2,3] <- rowMeans(CV_output_3[8,])
mean_df_3[3,3] <- rowMeans(CV_output_3[9,])
mean_df_3[4,3] <- rowMeans(CV_output_3[10,])
mean_df_3[5,3] <- rowMeans(CV_output_3[11,])
mean_df_3[6,3] <- rowMeans(CV_output_3[12,])

mean_df_3[1,4] <- rowMeans(CV_output_3[25,])
mean_df_3[2,4] <- rowMeans(CV_output_3[26,])
mean_df_3[3,4] <- rowMeans(CV_output_3[27,])
mean_df_3[4,4] <- rowMeans(CV_output_3[28,])
mean_df_3[5,4] <- rowMeans(CV_output_3[29,])
mean_df_3[6,4] <- rowMeans(CV_output_3[30,])

mean_df_3[1,5] <- rowMeans(CV_output_3[43,])
mean_df_3[2,5] <- rowMeans(CV_output_3[44,])
mean_df_3[3,5] <- rowMeans(CV_output_3[45,])
mean_df_3[4,5] <- rowMeans(CV_output_3[46,])
mean_df_3[5,5] <- rowMeans(CV_output_3[47,])
mean_df_3[6,5] <- rowMeans(CV_output_3[48,])

#this subtracts the did from the means
mean_df_3[1,6] <- mean_df_3[1,2] - mean_df_3[1,5]
mean_df_3[2,6] <- mean_df_3[2,2] - mean_df_3[2,5]
mean_df_3[3,6] <- mean_df_3[3,2] - mean_df_3[3,5]
mean_df_3[4,6] <- mean_df_3[4,2] - mean_df_3[4,5]
mean_df_3[5,6] <- mean_df_3[5,2] - mean_df_3[5,5]
mean_df_3[6,6] <- mean_df_3[6,2] - mean_df_3[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_3[1,7] <- (mean_df_3[1,3] - mean_df_3[1,2])
mean_df_3[2,7] <- (mean_df_3[2,3] - mean_df_3[2,2])
mean_df_3[3,7] <- (mean_df_3[3,3] - mean_df_3[3,2])
mean_df_3[4,7] <- (mean_df_3[4,3] - mean_df_3[4,2])
mean_df_3[5,7] <- (mean_df_3[5,3] - mean_df_3[5,2])
mean_df_3[6,7] <- (mean_df_3[6,3] - mean_df_3[6,2])

#mean diff of sc
mean_df_3[1,8] <- (mean_df_3[1,4] - mean_df_3[1,2])
mean_df_3[2,8] <- (mean_df_3[2,4] - mean_df_3[2,2])
mean_df_3[3,8] <- (mean_df_3[3,4] - mean_df_3[3,2])
mean_df_3[4,8] <- (mean_df_3[4,4] - mean_df_3[4,2])
mean_df_3[5,8] <- (mean_df_3[5,4] - mean_df_3[5,2])
mean_df_3[6,8] <- (mean_df_3[6,4] - mean_df_3[6,2])

#mean diff of did
mean_df_3[1,9] <- (mean_df_3[1,6] - mean_df_3[1,2])
mean_df_3[2,9] <- (mean_df_3[2,6] - mean_df_3[2,2])
mean_df_3[3,9] <- (mean_df_3[3,6] - mean_df_3[3,2])
mean_df_3[4,9] <- (mean_df_3[4,6] - mean_df_3[4,2])
mean_df_3[5,9] <- (mean_df_3[5,6] - mean_df_3[5,2])
mean_df_3[6,9] <- (mean_df_3[6,6] - mean_df_3[6,2])



# here is the CV for fold 4 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_4 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_4 <- wide_control_df[,c(1:176,265:888,177:264)]

for(s in 1:t) {
  CV_output_4[1:6,] <- wide_control_df_4[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_4
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_4[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_4[13:18,1:88] <- CV_output_4[1:6,1:88] - CV_output_4[7:12,1:88]
  CV_output_4[19:24,1:88] <- (CV_output_4[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_4[25:30,1:88] <- sc_CV 
  CV_output_4[31:36,1:88] <- CV_output_4[1:6,1:88]  - CV_output_4[25:30,1:88]   
  CV_output_4[37:42,1:88] <- (CV_output_4[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_4
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_4[43,z] <-sample_coeff[1,1]
  }
  CV_output_4[49:54,1:88] <- CV_output_4[1:6,1:88]  - CV_output_4[43:48,1:88]   
  CV_output_4[55:60,1:88] <- (CV_output_4[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_4 <- as.data.frame(matrix(NA,6,5))
mean_df_4[,1] <- c(2011:2016)
mean_df_4[1,2] <- rowMeans(CV_output_4[1,])
mean_df_4[2,2] <- rowMeans(CV_output_4[2,])
mean_df_4[3,2] <- rowMeans(CV_output_4[3,])
mean_df_4[4,2] <- rowMeans(CV_output_4[4,])
mean_df_4[5,2] <- rowMeans(CV_output_4[5,])
mean_df_4[6,2] <- rowMeans(CV_output_4[6,])

mean_df_4[1,3] <- rowMeans(CV_output_4[7,])
mean_df_4[2,3] <- rowMeans(CV_output_4[8,])
mean_df_4[3,3] <- rowMeans(CV_output_4[9,])
mean_df_4[4,3] <- rowMeans(CV_output_4[10,])
mean_df_4[5,3] <- rowMeans(CV_output_4[11,])
mean_df_4[6,3] <- rowMeans(CV_output_4[12,])

mean_df_4[1,4] <- rowMeans(CV_output_4[25,])
mean_df_4[2,4] <- rowMeans(CV_output_4[26,])
mean_df_4[3,4] <- rowMeans(CV_output_4[27,])
mean_df_4[4,4] <- rowMeans(CV_output_4[28,])
mean_df_4[5,4] <- rowMeans(CV_output_4[29,])
mean_df_4[6,4] <- rowMeans(CV_output_4[30,])

mean_df_4[1,5] <- rowMeans(CV_output_4[43,])
mean_df_4[2,5] <- rowMeans(CV_output_4[44,])
mean_df_4[3,5] <- rowMeans(CV_output_4[45,])
mean_df_4[4,5] <- rowMeans(CV_output_4[46,])
mean_df_4[5,5] <- rowMeans(CV_output_4[47,])
mean_df_4[6,5] <- rowMeans(CV_output_4[48,])

#this subtracts the did from the means
mean_df_4[1,6] <- mean_df_4[1,2] - mean_df_4[1,5]
mean_df_4[2,6] <- mean_df_4[2,2] - mean_df_4[2,5]
mean_df_4[3,6] <- mean_df_4[3,2] - mean_df_4[3,5]
mean_df_4[4,6] <- mean_df_4[4,2] - mean_df_4[4,5]
mean_df_4[5,6] <- mean_df_4[5,2] - mean_df_4[5,5]
mean_df_4[6,6] <- mean_df_4[6,2] - mean_df_4[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_4[1,7] <- (mean_df_4[1,3] - mean_df_4[1,2])
mean_df_4[2,7] <- (mean_df_4[2,3] - mean_df_4[2,2])
mean_df_4[3,7] <- (mean_df_4[3,3] - mean_df_4[3,2])
mean_df_4[4,7] <- (mean_df_4[4,3] - mean_df_4[4,2])
mean_df_4[5,7] <- (mean_df_4[5,3] - mean_df_4[5,2])
mean_df_4[6,7] <- (mean_df_4[6,3] - mean_df_4[6,2])

#mean diff of sc
mean_df_4[1,8] <- (mean_df_4[1,4] - mean_df_4[1,2])
mean_df_4[2,8] <- (mean_df_4[2,4] - mean_df_4[2,2])
mean_df_4[3,8] <- (mean_df_4[3,4] - mean_df_4[3,2])
mean_df_4[4,8] <- (mean_df_4[4,4] - mean_df_4[4,2])
mean_df_4[5,8] <- (mean_df_4[5,4] - mean_df_4[5,2])
mean_df_4[6,8] <- (mean_df_4[6,4] - mean_df_4[6,2])

#mean diff of did
mean_df_4[1,9] <- (mean_df_4[1,6] - mean_df_4[1,2])
mean_df_4[2,9] <- (mean_df_4[2,6] - mean_df_4[2,2])
mean_df_4[3,9] <- (mean_df_4[3,6] - mean_df_4[3,2])
mean_df_4[4,9] <- (mean_df_4[4,6] - mean_df_4[4,2])
mean_df_4[5,9] <- (mean_df_4[5,6] - mean_df_4[5,2])
mean_df_4[6,9] <- (mean_df_4[6,6] - mean_df_4[6,2])


# Fold 5 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_5 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_5 <- wide_control_df[,c(1:264,353:888,265:352)]

for(s in 1:t) {
  CV_output_5[1:6,] <- wide_control_df_5[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_5
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_5[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_5[13:18,1:88] <- CV_output_5[1:6,1:88] - CV_output_5[7:12,1:88]
  CV_output_5[19:24,1:88] <- (CV_output_5[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_5[25:30,1:88] <- sc_CV 
  CV_output_5[31:36,1:88] <- CV_output_5[1:6,1:88]  - CV_output_5[25:30,1:88]   
  CV_output_5[37:42,1:88] <- (CV_output_5[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_5
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_5[43,z] <-sample_coeff[1,1]
  }
  CV_output_5[49:54,1:88] <- CV_output_5[1:6,1:88]  - CV_output_5[43:48,1:88]   
  CV_output_5[55:60,1:88] <- (CV_output_5[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_5 <- as.data.frame(matrix(NA,6,5))
mean_df_5[,1] <- c(2011:2016)
mean_df_5[1,2] <- rowMeans(CV_output_5[1,])
mean_df_5[2,2] <- rowMeans(CV_output_5[2,])
mean_df_5[3,2] <- rowMeans(CV_output_5[3,])
mean_df_5[4,2] <- rowMeans(CV_output_5[4,])
mean_df_5[5,2] <- rowMeans(CV_output_5[5,])
mean_df_5[6,2] <- rowMeans(CV_output_5[6,])

mean_df_5[1,3] <- rowMeans(CV_output_5[7,])
mean_df_5[2,3] <- rowMeans(CV_output_5[8,])
mean_df_5[3,3] <- rowMeans(CV_output_5[9,])
mean_df_5[4,3] <- rowMeans(CV_output_5[10,])
mean_df_5[5,3] <- rowMeans(CV_output_5[11,])
mean_df_5[6,3] <- rowMeans(CV_output_5[12,])

mean_df_5[1,4] <- rowMeans(CV_output_5[25,])
mean_df_5[2,4] <- rowMeans(CV_output_5[26,])
mean_df_5[3,4] <- rowMeans(CV_output_5[27,])
mean_df_5[4,4] <- rowMeans(CV_output_5[28,])
mean_df_5[5,4] <- rowMeans(CV_output_5[29,])
mean_df_5[6,4] <- rowMeans(CV_output_5[30,])

mean_df_5[1,5] <- rowMeans(CV_output_5[43,])
mean_df_5[2,5] <- rowMeans(CV_output_5[44,])
mean_df_5[3,5] <- rowMeans(CV_output_5[45,])
mean_df_5[4,5] <- rowMeans(CV_output_5[46,])
mean_df_5[5,5] <- rowMeans(CV_output_5[47,])
mean_df_5[6,5] <- rowMeans(CV_output_5[48,])

#this subtracts the did from the means
mean_df_5[1,6] <- mean_df_5[1,2] - mean_df_5[1,5]
mean_df_5[2,6] <- mean_df_5[2,2] - mean_df_5[2,5]
mean_df_5[3,6] <- mean_df_5[3,2] - mean_df_5[3,5]
mean_df_5[4,6] <- mean_df_5[4,2] - mean_df_5[4,5]
mean_df_5[5,6] <- mean_df_5[5,2] - mean_df_5[5,5]
mean_df_5[6,6] <- mean_df_5[6,2] - mean_df_5[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_5[1,7] <- (mean_df_5[1,3] - mean_df_5[1,2])
mean_df_5[2,7] <- (mean_df_5[2,3] - mean_df_5[2,2])
mean_df_5[3,7] <- (mean_df_5[3,3] - mean_df_5[3,2])
mean_df_5[4,7] <- (mean_df_5[4,3] - mean_df_5[4,2])
mean_df_5[5,7] <- (mean_df_5[5,3] - mean_df_5[5,2])
mean_df_5[6,7] <- (mean_df_5[6,3] - mean_df_5[6,2])

#mean diff of sc
mean_df_5[1,8] <- (mean_df_5[1,4] - mean_df_5[1,2])
mean_df_5[2,8] <- (mean_df_5[2,4] - mean_df_5[2,2])
mean_df_5[3,8] <- (mean_df_5[3,4] - mean_df_5[3,2])
mean_df_5[4,8] <- (mean_df_5[4,4] - mean_df_5[4,2])
mean_df_5[5,8] <- (mean_df_5[5,4] - mean_df_5[5,2])
mean_df_5[6,8] <- (mean_df_5[6,4] - mean_df_5[6,2])

#mean diff of did
mean_df_5[1,9] <- (mean_df_5[1,6] - mean_df_5[1,2])
mean_df_5[2,9] <- (mean_df_5[2,6] - mean_df_5[2,2])
mean_df_5[3,9] <- (mean_df_5[3,6] - mean_df_5[3,2])
mean_df_5[4,9] <- (mean_df_5[4,6] - mean_df_5[4,2])
mean_df_5[5,9] <- (mean_df_5[5,6] - mean_df_5[5,2])
mean_df_5[6,9] <- (mean_df_5[6,6] - mean_df_5[6,2])


# Fold 6 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_6 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_6 <- wide_control_df[,c(1:352,441:888,353:440)]

for(s in 1:t) {
  CV_output_6[1:6,] <- wide_control_df_6[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_6
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_6[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_6[13:18,1:88] <- CV_output_6[1:6,1:88] - CV_output_6[7:12,1:88]
  CV_output_6[19:24,1:88] <- (CV_output_6[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_6[25:30,1:88] <- sc_CV 
  CV_output_6[31:36,1:88] <- CV_output_6[1:6,1:88]  - CV_output_6[25:30,1:88]   
  CV_output_6[37:42,1:88] <- (CV_output_6[31:36,1:88])^2
  
  # DiD 
  long_c_df_boot_cv <- wide_control_df_6
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_6[43,z] <-sample_coeff[1,1]
  }
  CV_output_6[49:54,1:88] <- CV_output_6[1:6,1:88]  - CV_output_6[43:48,1:88]   
  CV_output_6[55:60,1:88] <- (CV_output_6[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_6 <- as.data.frame(matrix(NA,6,5))
mean_df_6[,1] <- c(2011:2016)
mean_df_6[1,2] <- rowMeans(CV_output_6[1,])
mean_df_6[2,2] <- rowMeans(CV_output_6[2,])
mean_df_6[3,2] <- rowMeans(CV_output_6[3,])
mean_df_6[4,2] <- rowMeans(CV_output_6[4,])
mean_df_6[5,2] <- rowMeans(CV_output_6[5,])
mean_df_6[6,2] <- rowMeans(CV_output_6[6,])

mean_df_6[1,3] <- rowMeans(CV_output_6[7,])
mean_df_6[2,3] <- rowMeans(CV_output_6[8,])
mean_df_6[3,3] <- rowMeans(CV_output_6[9,])
mean_df_6[4,3] <- rowMeans(CV_output_6[10,])
mean_df_6[5,3] <- rowMeans(CV_output_6[11,])
mean_df_6[6,3] <- rowMeans(CV_output_6[12,])

mean_df_6[1,4] <- rowMeans(CV_output_6[25,])
mean_df_6[2,4] <- rowMeans(CV_output_6[26,])
mean_df_6[3,4] <- rowMeans(CV_output_6[27,])
mean_df_6[4,4] <- rowMeans(CV_output_6[28,])
mean_df_6[5,4] <- rowMeans(CV_output_6[29,])
mean_df_6[6,4] <- rowMeans(CV_output_6[30,])

mean_df_6[1,5] <- rowMeans(CV_output_6[43,])
mean_df_6[2,5] <- rowMeans(CV_output_6[44,])
mean_df_6[3,5] <- rowMeans(CV_output_6[45,])
mean_df_6[4,5] <- rowMeans(CV_output_6[46,])
mean_df_6[5,5] <- rowMeans(CV_output_6[47,])
mean_df_6[6,5] <- rowMeans(CV_output_6[48,])

#this subtracts the did from the means
mean_df_6[1,6] <- mean_df_6[1,2] - mean_df_6[1,5]
mean_df_6[2,6] <- mean_df_6[2,2] - mean_df_6[2,5]
mean_df_6[3,6] <- mean_df_6[3,2] - mean_df_6[3,5]
mean_df_6[4,6] <- mean_df_6[4,2] - mean_df_6[4,5]
mean_df_6[5,6] <- mean_df_6[5,2] - mean_df_6[5,5]
mean_df_6[6,6] <- mean_df_6[6,2] - mean_df_6[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_6[1,7] <- (mean_df_6[1,3] - mean_df_6[1,2])
mean_df_6[2,7] <- (mean_df_6[2,3] - mean_df_6[2,2])
mean_df_6[3,7] <- (mean_df_6[3,3] - mean_df_6[3,2])
mean_df_6[4,7] <- (mean_df_6[4,3] - mean_df_6[4,2])
mean_df_6[5,7] <- (mean_df_6[5,3] - mean_df_6[5,2])
mean_df_6[6,7] <- (mean_df_6[6,3] - mean_df_6[6,2])

#mean diff of sc
mean_df_6[1,8] <- (mean_df_6[1,4] - mean_df_6[1,2])
mean_df_6[2,8] <- (mean_df_6[2,4] - mean_df_6[2,2])
mean_df_6[3,8] <- (mean_df_6[3,4] - mean_df_6[3,2])
mean_df_6[4,8] <- (mean_df_6[4,4] - mean_df_6[4,2])
mean_df_6[5,8] <- (mean_df_6[5,4] - mean_df_6[5,2])
mean_df_6[6,8] <- (mean_df_6[6,4] - mean_df_6[6,2])

#mean diff of did
mean_df_6[1,9] <- (mean_df_6[1,6] - mean_df_6[1,2])
mean_df_6[2,9] <- (mean_df_6[2,6] - mean_df_6[2,2])
mean_df_6[3,9] <- (mean_df_6[3,6] - mean_df_6[3,2])
mean_df_6[4,9] <- (mean_df_6[4,6] - mean_df_6[4,2])
mean_df_6[5,9] <- (mean_df_6[5,6] - mean_df_6[5,2])
mean_df_6[6,9] <- (mean_df_6[6,6] - mean_df_6[6,2])


# Fold 7 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_7 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_7 <- wide_control_df[,c(1:440,529:888,441:528)]

for(s in 1:t) {
  CV_output_7[1:6,] <- wide_control_df_7[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_7
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  sqrt((1/sum(1-mask_basic_control_test)) * sum(model_with_both_basic_control_test$msk_err^2))
  
  CV_output_7[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_7[13:18,1:88] <- CV_output_7[1:6,1:88] - CV_output_7[7:12,1:88]
  CV_output_7[19:24,1:88] <- (CV_output_7[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_7[25:30,1:88] <- sc_CV 
  CV_output_7[31:36,1:88] <- CV_output_7[1:6,1:88]  - CV_output_7[25:30,1:88]   
  CV_output_7[37:42,1:88] <- (CV_output_7[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_7
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_7[43,z] <-sample_coeff[1,1]
  }
  CV_output_7[49:54,1:88] <- CV_output_7[1:6,1:88]  - CV_output_7[43:48,1:88]   
  CV_output_7[55:60,1:88] <- (CV_output_7[43:48,1:88])^2
}

proc.time() - ptm

#plot df and chart
mean_df_7 <- as.data.frame(matrix(NA,6,5))
mean_df_7[,1] <- c(2011:2016)
mean_df_7[1,2] <- rowMeans(CV_output_7[1,])
mean_df_7[2,2] <- rowMeans(CV_output_7[2,])
mean_df_7[3,2] <- rowMeans(CV_output_7[3,])
mean_df_7[4,2] <- rowMeans(CV_output_7[4,])
mean_df_7[5,2] <- rowMeans(CV_output_7[5,])
mean_df_7[6,2] <- rowMeans(CV_output_7[6,])

mean_df_7[1,3] <- rowMeans(CV_output_7[7,])
mean_df_7[2,3] <- rowMeans(CV_output_7[8,])
mean_df_7[3,3] <- rowMeans(CV_output_7[9,])
mean_df_7[4,3] <- rowMeans(CV_output_7[10,])
mean_df_7[5,3] <- rowMeans(CV_output_7[11,])
mean_df_7[6,3] <- rowMeans(CV_output_7[12,])

mean_df_7[1,4] <- rowMeans(CV_output_7[25,])
mean_df_7[2,4] <- rowMeans(CV_output_7[26,])
mean_df_7[3,4] <- rowMeans(CV_output_7[27,])
mean_df_7[4,4] <- rowMeans(CV_output_7[28,])
mean_df_7[5,4] <- rowMeans(CV_output_7[29,])
mean_df_7[6,4] <- rowMeans(CV_output_7[30,])

mean_df_7[1,5] <- rowMeans(CV_output_7[43,])
mean_df_7[2,5] <- rowMeans(CV_output_7[44,])
mean_df_7[3,5] <- rowMeans(CV_output_7[45,])
mean_df_7[4,5] <- rowMeans(CV_output_7[46,])
mean_df_7[5,5] <- rowMeans(CV_output_7[47,])
mean_df_7[6,5] <- rowMeans(CV_output_7[48,])

#this subtracts the did from the means
mean_df_7[1,6] <- mean_df_7[1,2] - mean_df_7[1,5]
mean_df_7[2,6] <- mean_df_7[2,2] - mean_df_7[2,5]
mean_df_7[3,6] <- mean_df_7[3,2] - mean_df_7[3,5]
mean_df_7[4,6] <- mean_df_7[4,2] - mean_df_7[4,5]
mean_df_7[5,6] <- mean_df_7[5,2] - mean_df_7[5,5]
mean_df_7[6,6] <- mean_df_7[6,2] - mean_df_7[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_7[1,7] <- (mean_df_7[1,3] - mean_df_7[1,2])
mean_df_7[2,7] <- (mean_df_7[2,3] - mean_df_7[2,2])
mean_df_7[3,7] <- (mean_df_7[3,3] - mean_df_7[3,2])
mean_df_7[4,7] <- (mean_df_7[4,3] - mean_df_7[4,2])
mean_df_7[5,7] <- (mean_df_7[5,3] - mean_df_7[5,2])
mean_df_7[6,7] <- (mean_df_7[6,3] - mean_df_7[6,2])

#mean diff of sc
mean_df_7[1,8] <- (mean_df_7[1,4] - mean_df_7[1,2])
mean_df_7[2,8] <- (mean_df_7[2,4] - mean_df_7[2,2])
mean_df_7[3,8] <- (mean_df_7[3,4] - mean_df_7[3,2])
mean_df_7[4,8] <- (mean_df_7[4,4] - mean_df_7[4,2])
mean_df_7[5,8] <- (mean_df_7[5,4] - mean_df_7[5,2])
mean_df_7[6,8] <- (mean_df_7[6,4] - mean_df_7[6,2])

#mean diff of did
mean_df_7[1,9] <- (mean_df_7[1,6] - mean_df_7[1,2])
mean_df_7[2,9] <- (mean_df_7[2,6] - mean_df_7[2,2])
mean_df_7[3,9] <- (mean_df_7[3,6] - mean_df_7[3,2])
mean_df_7[4,9] <- (mean_df_7[4,6] - mean_df_7[4,2])
mean_df_7[5,9] <- (mean_df_7[5,6] - mean_df_7[5,2])
mean_df_7[6,9] <- (mean_df_7[6,6] - mean_df_7[6,2])



# Fold 8 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_8 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_8 <- wide_control_df[,c(1:528,617:888,529:616)]

for(s in 1:t) {
  CV_output_8[1:6,] <- wide_control_df_8[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_8
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_8[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_8[13:18,1:88] <- CV_output_8[1:6,1:88] - CV_output_8[7:12,1:88]
  CV_output_8[19:24,1:88] <- (CV_output_8[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_8[25:30,1:88] <- sc_CV 
  CV_output_8[31:36,1:88] <- CV_output_8[1:6,1:88]  - CV_output_8[25:30,1:88]   
  CV_output_8[37:42,1:88] <- (CV_output_8[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_8
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_8[43,z] <-sample_coeff[1,1]
  }
  CV_output_8[49:54,1:88] <- CV_output_8[1:6,1:88]  - CV_output_8[43:48,1:88]   
  CV_output_8[55:60,1:88] <- (CV_output_8[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_8 <- as.data.frame(matrix(NA,6,5))
mean_df_8[,1] <- c(2011:2016)
mean_df_8[1,2] <- rowMeans(CV_output_8[1,])
mean_df_8[2,2] <- rowMeans(CV_output_8[2,])
mean_df_8[3,2] <- rowMeans(CV_output_8[3,])
mean_df_8[4,2] <- rowMeans(CV_output_8[4,])
mean_df_8[5,2] <- rowMeans(CV_output_8[5,])
mean_df_8[6,2] <- rowMeans(CV_output_8[6,])

mean_df_8[1,3] <- rowMeans(CV_output_8[7,])
mean_df_8[2,3] <- rowMeans(CV_output_8[8,])
mean_df_8[3,3] <- rowMeans(CV_output_8[9,])
mean_df_8[4,3] <- rowMeans(CV_output_8[10,])
mean_df_8[5,3] <- rowMeans(CV_output_8[11,])
mean_df_8[6,3] <- rowMeans(CV_output_8[12,])

mean_df_8[1,4] <- rowMeans(CV_output_8[25,])
mean_df_8[2,4] <- rowMeans(CV_output_8[26,])
mean_df_8[3,4] <- rowMeans(CV_output_8[27,])
mean_df_8[4,4] <- rowMeans(CV_output_8[28,])
mean_df_8[5,4] <- rowMeans(CV_output_8[29,])
mean_df_8[6,4] <- rowMeans(CV_output_8[30,])

mean_df_8[1,5] <- rowMeans(CV_output_8[43,])
mean_df_8[2,5] <- rowMeans(CV_output_8[44,])
mean_df_8[3,5] <- rowMeans(CV_output_8[45,])
mean_df_8[4,5] <- rowMeans(CV_output_8[46,])
mean_df_8[5,5] <- rowMeans(CV_output_8[47,])
mean_df_8[6,5] <- rowMeans(CV_output_8[48,])

#this subtracts the did from the means
mean_df_8[1,6] <- mean_df_8[1,2] - mean_df_8[1,5]
mean_df_8[2,6] <- mean_df_8[2,2] - mean_df_8[2,5]
mean_df_8[3,6] <- mean_df_8[3,2] - mean_df_8[3,5]
mean_df_8[4,6] <- mean_df_8[4,2] - mean_df_8[4,5]
mean_df_8[5,6] <- mean_df_8[5,2] - mean_df_8[5,5]
mean_df_8[6,6] <- mean_df_8[6,2] - mean_df_8[6,5]

#
#mean diff of mc
mean_df_8[1,7] <- (mean_df_8[1,3] - mean_df_8[1,2])
mean_df_8[2,7] <- (mean_df_8[2,3] - mean_df_8[2,2])
mean_df_8[3,7] <- (mean_df_8[3,3] - mean_df_8[3,2])
mean_df_8[4,7] <- (mean_df_8[4,3] - mean_df_8[4,2])
mean_df_8[5,7] <- (mean_df_8[5,3] - mean_df_8[5,2])
mean_df_8[6,7] <- (mean_df_8[6,3] - mean_df_8[6,2])

#mean diff of sc
mean_df_8[1,8] <- (mean_df_8[1,4] - mean_df_8[1,2])
mean_df_8[2,8] <- (mean_df_8[2,4] - mean_df_8[2,2])
mean_df_8[3,8] <- (mean_df_8[3,4] - mean_df_8[3,2])
mean_df_8[4,8] <- (mean_df_8[4,4] - mean_df_8[4,2])
mean_df_8[5,8] <- (mean_df_8[5,4] - mean_df_8[5,2])
mean_df_8[6,8] <- (mean_df_8[6,4] - mean_df_8[6,2])

#mean diff of did
mean_df_8[1,9] <- (mean_df_8[1,6] - mean_df_8[1,2])
mean_df_8[2,9] <- (mean_df_8[2,6] - mean_df_8[2,2])
mean_df_8[3,9] <- (mean_df_8[3,6] - mean_df_8[3,2])
mean_df_8[4,9] <- (mean_df_8[4,6] - mean_df_8[4,2])
mean_df_8[5,9] <- (mean_df_8[5,6] - mean_df_8[5,2])
mean_df_8[6,9] <- (mean_df_8[6,6] - mean_df_8[6,2])



# Fold 9 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_9 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_9 <- wide_control_df[,c(1:616,705:888,617:704)]

for(s in 1:t) {
  CV_output_9[1:6,] <- wide_control_df_9[6:11,801:888] # 
  
  # MC
  sample_matrix <- wide_control_df_9
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_9[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_9[13:18,1:88] <- CV_output_9[1:6,1:88] - CV_output_9[7:12,1:88]
  CV_output_9[19:24,1:88] <- (CV_output_9[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_9[25:30,1:88] <- sc_CV 
  CV_output_9[31:36,1:88] <- CV_output_9[1:6,1:88]  - CV_output_9[25:30,1:88]   
  CV_output_9[37:42,1:88] <- (CV_output_9[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_9
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_9[43,z] <-sample_coeff[1,1]
  }
  CV_output_9[49:54,1:88] <- CV_output_9[1:6,1:88]  - CV_output_9[43:48,1:88]   
  CV_output_9[55:60,1:88] <- (CV_output_9[43:48,1:88])^2
}

proc.time() - ptm

# Plot df and chart
mean_df_9 <- as.data.frame(matrix(NA,6,5))
mean_df_9[,1] <- c(2011:2016)
mean_df_9[1,2] <- rowMeans(CV_output_9[1,])
mean_df_9[2,2] <- rowMeans(CV_output_9[2,])
mean_df_9[3,2] <- rowMeans(CV_output_9[3,])
mean_df_9[4,2] <- rowMeans(CV_output_9[4,])
mean_df_9[5,2] <- rowMeans(CV_output_9[5,])
mean_df_9[6,2] <- rowMeans(CV_output_9[6,])

mean_df_9[1,3] <- rowMeans(CV_output_9[7,])
mean_df_9[2,3] <- rowMeans(CV_output_9[8,])
mean_df_9[3,3] <- rowMeans(CV_output_9[9,])
mean_df_9[4,3] <- rowMeans(CV_output_9[10,])
mean_df_9[5,3] <- rowMeans(CV_output_9[11,])
mean_df_9[6,3] <- rowMeans(CV_output_9[12,])

mean_df_9[1,4] <- rowMeans(CV_output_9[25,])
mean_df_9[2,4] <- rowMeans(CV_output_9[26,])
mean_df_9[3,4] <- rowMeans(CV_output_9[27,])
mean_df_9[4,4] <- rowMeans(CV_output_9[28,])
mean_df_9[5,4] <- rowMeans(CV_output_9[29,])
mean_df_9[6,4] <- rowMeans(CV_output_9[30,])

mean_df_9[1,5] <- rowMeans(CV_output_9[43,])
mean_df_9[2,5] <- rowMeans(CV_output_9[44,])
mean_df_9[3,5] <- rowMeans(CV_output_9[45,])
mean_df_9[4,5] <- rowMeans(CV_output_9[46,])
mean_df_9[5,5] <- rowMeans(CV_output_9[47,])
mean_df_9[6,5] <- rowMeans(CV_output_9[48,])

#this subtracts the did from the means
mean_df_9[1,6] <- mean_df_9[1,2] - mean_df_9[1,5]
mean_df_9[2,6] <- mean_df_9[2,2] - mean_df_9[2,5]
mean_df_9[3,6] <- mean_df_9[3,2] - mean_df_9[3,5]
mean_df_9[4,6] <- mean_df_9[4,2] - mean_df_9[4,5]
mean_df_9[5,6] <- mean_df_9[5,2] - mean_df_9[5,5]
mean_df_9[6,6] <- mean_df_9[6,2] - mean_df_9[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_9[1,7] <- (mean_df_9[1,3] - mean_df_9[1,2])
mean_df_9[2,7] <- (mean_df_9[2,3] - mean_df_9[2,2])
mean_df_9[3,7] <- (mean_df_9[3,3] - mean_df_9[3,2])
mean_df_9[4,7] <- (mean_df_9[4,3] - mean_df_9[4,2])
mean_df_9[5,7] <- (mean_df_9[5,3] - mean_df_9[5,2])
mean_df_9[6,7] <- (mean_df_9[6,3] - mean_df_9[6,2])

#mean diff of sc
mean_df_9[1,8] <- (mean_df_9[1,4] - mean_df_9[1,2])
mean_df_9[2,8] <- (mean_df_9[2,4] - mean_df_9[2,2])
mean_df_9[3,8] <- (mean_df_9[3,4] - mean_df_9[3,2])
mean_df_9[4,8] <- (mean_df_9[4,4] - mean_df_9[4,2])
mean_df_9[5,8] <- (mean_df_9[5,4] - mean_df_9[5,2])
mean_df_9[6,8] <- (mean_df_9[6,4] - mean_df_9[6,2])

#mean diff of did
mean_df_9[1,9] <- (mean_df_9[1,6] - mean_df_9[1,2])
mean_df_9[2,9] <- (mean_df_9[2,6] - mean_df_9[2,2])
mean_df_9[3,9] <- (mean_df_9[3,6] - mean_df_9[3,2])
mean_df_9[4,9] <- (mean_df_9[4,6] - mean_df_9[4,2])
mean_df_9[5,9] <- (mean_df_9[5,6] - mean_df_9[5,2])
mean_df_9[6,9] <- (mean_df_9[6,6] - mean_df_9[6,2])


# Fold 10 ###########
ptm <- proc.time()
set.seed(789)

t = 1
v = 6
y = 88
l = 6
r = 88
CV_output_10 <- as.data.frame(matrix(NA, 60, 88))
sc_CV <- as.data.frame(matrix(NA,6,88))

wide_control_df_10 <- wide_control_df[,c(1:704,793:888,705:792)]

for(s in 1:t) {
  CV_output_10[1:6,] <- wide_control_df_10[6:11,801:888] #this will change in future iterations. 
  
  # MC
  sample_matrix <- wide_control_df_10
  
  N_basic_control_test <- 888
  T_basic_control_test <- 11
  M_basic_control_test <- sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,11,888)
  mask_basic_control_test[6:11,801:888] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  CV_output_10[7:12,] <- model_with_both_basic_control_test$est[6:11,801:888]
  CV_output_10[13:18,1:88] <- CV_output_10[1:6,1:88] - CV_output_10[7:12,1:88]
  CV_output_10[19:24,1:88] <- (CV_output_10[13:18,1:88])^2
  
  # SC
  for(u in 1:v) {  
    for(z in 1:y) {
      X_prod <- t(as.matrix(sample_matrix[1:5,1:800]))
      Y_prod <- t(as.matrix(sample_matrix[5+u,1:800]))
      
      en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
      
      post_prod <- as.matrix(t(sample_matrix[1:5,800+z]))
      pred_prod <- predict(en_prod, newx = post_prod, 
                           s = (en_prod$lambda.min))
      sc_CV[u,z] <- pred_prod
    }
  }
  CV_output_10[25:30,1:88] <- sc_CV 
  CV_output_10[31:36,1:88] <- CV_output_10[1:6,1:88]  - CV_output_10[25:30,1:88]   
  CV_output_10[37:42,1:88] <- (CV_output_10[31:36,1:88])^2
  
  # DiD
  long_c_df_boot_cv <- wide_control_df_10
  long_c_df_boot_cv$year <- c(2006:2016)
  
  long_c_df_boot_cv <- reshape(long_c_df_boot_cv,direction='long', varying = list(names(long_c_df_boot_cv)[1:888]), timevar ='year')
  long_c_df_boot_cv <- long_c_df_boot_cv[,2:3]
  long_c_df_boot_cv$year <- c(2006:2016)
  long_c_df_boot_cv$treated_unit <- 0
  long_c_df_boot_cv[8801:9768,4] <- 1
  long_c_df_boot_cv$treat_year <- c(0,0,0,0,0,1,1,1,1,1,1)
  long_c_df_boot_cv$treat_id <- long_c_df_boot_cv[,4] * long_c_df_boot_cv[,5]
  long_c_df_boot_cv$unit_id <- rep(1:888, each = 11)
  colnames(long_c_df_boot_cv)[1] <- "WI"
  
  for(z in 1:y) {
    sample_long_c_df_boot_cv <- rbind(long_c_df_boot_cv[1:8800,], long_c_df_boot_cv[(8790:8800)+(z*11),])
    
    sample_long_c_df_boot_cv_2016 <- subset(sample_long_c_df_boot_cv)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[48,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2015 <- subset(sample_long_c_df_boot_cv, year < 2016)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2015))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[47,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2014 <- subset(sample_long_c_df_boot_cv, year < 2015)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2014))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[46,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2013 <- subset(sample_long_c_df_boot_cv, year < 2014)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2013))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[45,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2012 <- subset(sample_long_c_df_boot_cv, year < 2013)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2012))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[44,z] <-sample_coeff[1,1]
    
    sample_long_c_df_boot_cv_2011 <- subset(sample_long_c_df_boot_cv, year < 2012)
    sample_did <- summary(feols(WI ~ treat_id | unit_id + year, data=sample_long_c_df_boot_cv_2011))
    sample_coeff <- as.data.frame(coeftable(sample_did)[1])
    CV_output_10[43,z] <-sample_coeff[1,1]
  }
  CV_output_10[49:54,1:88] <- CV_output_10[1:6,1:88]  - CV_output_10[43:48,1:88]   
  CV_output_10[55:60,1:88] <- (CV_output_10[43:48,1:88])^2
}

proc.time() - ptm

#plot df and chart
mean_df_10 <- as.data.frame(matrix(NA,6,5))
mean_df_10[,1] <- c(2011:2016)
mean_df_10[1,2] <- rowMeans(CV_output_10[1,])
mean_df_10[2,2] <- rowMeans(CV_output_10[2,])
mean_df_10[3,2] <- rowMeans(CV_output_10[3,])
mean_df_10[4,2] <- rowMeans(CV_output_10[4,])
mean_df_10[5,2] <- rowMeans(CV_output_10[5,])
mean_df_10[6,2] <- rowMeans(CV_output_10[6,])

mean_df_10[1,3] <- rowMeans(CV_output_10[7,])
mean_df_10[2,3] <- rowMeans(CV_output_10[8,])
mean_df_10[3,3] <- rowMeans(CV_output_10[9,])
mean_df_10[4,3] <- rowMeans(CV_output_10[10,])
mean_df_10[5,3] <- rowMeans(CV_output_10[11,])
mean_df_10[6,3] <- rowMeans(CV_output_10[12,])

mean_df_10[1,4] <- rowMeans(CV_output_10[25,])
mean_df_10[2,4] <- rowMeans(CV_output_10[26,])
mean_df_10[3,4] <- rowMeans(CV_output_10[27,])
mean_df_10[4,4] <- rowMeans(CV_output_10[28,])
mean_df_10[5,4] <- rowMeans(CV_output_10[29,])
mean_df_10[6,4] <- rowMeans(CV_output_10[30,])

mean_df_10[1,5] <- rowMeans(CV_output_10[43,])
mean_df_10[2,5] <- rowMeans(CV_output_10[44,])
mean_df_10[3,5] <- rowMeans(CV_output_10[45,])
mean_df_10[4,5] <- rowMeans(CV_output_10[46,])
mean_df_10[5,5] <- rowMeans(CV_output_10[47,])
mean_df_10[6,5] <- rowMeans(CV_output_10[48,])

#this subtracts the did from the means
mean_df_10[1,6] <- mean_df_10[1,2] - mean_df_10[1,5]
mean_df_10[2,6] <- mean_df_10[2,2] - mean_df_10[2,5]
mean_df_10[3,6] <- mean_df_10[3,2] - mean_df_10[3,5]
mean_df_10[4,6] <- mean_df_10[4,2] - mean_df_10[4,5]
mean_df_10[5,6] <- mean_df_10[5,2] - mean_df_10[5,5]
mean_df_10[6,6] <- mean_df_10[6,2] - mean_df_10[6,5]

#now we are looking at the mean erro
#mean diff of mc
mean_df_10[1,7] <- (mean_df_10[1,3] - mean_df_10[1,2])
mean_df_10[2,7] <- (mean_df_10[2,3] - mean_df_10[2,2])
mean_df_10[3,7] <- (mean_df_10[3,3] - mean_df_10[3,2])
mean_df_10[4,7] <- (mean_df_10[4,3] - mean_df_10[4,2])
mean_df_10[5,7] <- (mean_df_10[5,3] - mean_df_10[5,2])
mean_df_10[6,7] <- (mean_df_10[6,3] - mean_df_10[6,2])

#mean diff of sc
mean_df_10[1,8] <- (mean_df_10[1,4] - mean_df_10[1,2])
mean_df_10[2,8] <- (mean_df_10[2,4] - mean_df_10[2,2])
mean_df_10[3,8] <- (mean_df_10[3,4] - mean_df_10[3,2])
mean_df_10[4,8] <- (mean_df_10[4,4] - mean_df_10[4,2])
mean_df_10[5,8] <- (mean_df_10[5,4] - mean_df_10[5,2])
mean_df_10[6,8] <- (mean_df_10[6,4] - mean_df_10[6,2])

#mean diff of did
mean_df_10[1,9] <- (mean_df_10[1,6] - mean_df_10[1,2])
mean_df_10[2,9] <- (mean_df_10[2,6] - mean_df_10[2,2])
mean_df_10[3,9] <- (mean_df_10[3,6] - mean_df_10[3,2])
mean_df_10[4,9] <- (mean_df_10[4,6] - mean_df_10[4,2])
mean_df_10[5,9] <- (mean_df_10[5,6] - mean_df_10[5,2])
mean_df_10[6,9] <- (mean_df_10[6,6] - mean_df_10[6,2])

# here i am renaming all the df for the figure creation
write.csv(mean_df, "cv_fold_1.csv", row.names = FALSE)
write.csv(mean_df_2, "cv_fold_2.csv", row.names = FALSE)
write.csv(mean_df_3, "cv_fold_3.csv", row.names = FALSE)
write.csv(mean_df_4, "cv_fold_4.csv", row.names = FALSE)
write.csv(mean_df_5, "cv_fold_5.csv", row.names = FALSE)
write.csv(mean_df_6, "cv_fold_6.csv", row.names = FALSE)
write.csv(mean_df_7, "cv_fold_7.csv", row.names = FALSE)
write.csv(mean_df_8, "cv_fold_8.csv", row.names = FALSE)
write.csv(mean_df_9, "cv_fold_9.csv", row.names = FALSE)
write.csv(mean_df_10, "cv_fold_10.csv", row.names = FALSE)

# Now getting aggregate df
all_cv <- cbind(CV_output_1, CV_output_2, CV_output_3, CV_output_4, CV_output_5,
                CV_output_6, CV_output_7, CV_output_8, CV_output_9, CV_output_10)

write.csv(all_cv, "all_cv.csv")

# Now checking to see differences from the observed mean
cv_output <- all_cv[c(1:6,7:12,25:30,43:48),]
write.csv(cv_output, "cv_output.csv", row.names = FALSE)

# MC mean difference 
((rowMeans(cv_output[1,]) - rowMeans(cv_output[7,])) + (rowMeans(cv_output[2,]) - rowMeans(cv_output[8,])) +
    (rowMeans(cv_output[3,]) - rowMeans(cv_output[9,])) + (rowMeans(cv_output[4,]) - rowMeans(cv_output[10,])) +
    (rowMeans(cv_output[5,]) - rowMeans(cv_output[11,])) + (rowMeans(cv_output[6,]) - rowMeans(cv_output[12,]))) / 6

# SC-EN mean difference 
((rowMeans(cv_output[1,]) - rowMeans(cv_output[13,])) + (rowMeans(cv_output[2,]) - rowMeans(cv_output[14,])) +
    (rowMeans(cv_output[3,]) - rowMeans(cv_output[15,])) + (rowMeans(cv_output[4,]) - rowMeans(cv_output[16,])) +
    (rowMeans(cv_output[5,]) - rowMeans(cv_output[17,])) + (rowMeans(cv_output[6,]) - rowMeans(cv_output[18,]))) / 6

# DiD mean difference 
(rowMeans(cv_output[19,]) + rowMeans(cv_output[20,]) + rowMeans(cv_output[21,]) + 
    rowMeans(cv_output[22,]) + rowMeans(cv_output[23,]) + rowMeans(cv_output[24,])) / 6 

# RMSEs
all_year_rmse_mc <- sqrt((sum(all_cv[19,]) + sum(all_cv[20,]) + sum(all_cv[21,]) + sum(all_cv[22,]) + sum(all_cv[23,]) + sum(all_cv[24,])) / (6*888))
all_year_rmse_mc

all_year_rmse_sc <- sqrt((sum(all_cv[37,]) + sum(all_cv[38,]) + sum(all_cv[39,]) + sum(all_cv[40,]) + sum(all_cv[41,]) + sum(all_cv[42,])) / (6*888))
all_year_rmse_sc 

all_year_rmse_did <- sqrt((sum(all_cv[55,]) + sum(all_cv[56,]) + sum(all_cv[57,]) + sum(all_cv[58,]) + sum(all_cv[59,]) + sum(all_cv[60,])) / (6*888))
all_year_rmse_did  


