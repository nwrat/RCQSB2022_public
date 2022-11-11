#this script creates the pretrends test results
# need to remove 

rm(list=ls())

library(fixest)

setwd("")

tc_2km_ug <- read.csv("tc_2km_ug.csv")
tc_2km_ug_C <- subset(tc_2km_ug, treated == 0)
tc_2km_ug_T <- subset(tc_2km_ug, treated == 1)

tc_3km_ug <- read.csv("tc_3km_ug.csv")
tc_3km_ug_C <- subset(tc_3km_ug, treated == 0)
tc_3km_ug_T <- subset(tc_3km_ug, treated == 1)

tc_4km_ug <- read.csv("tc_4km_ug.csv")
tc_4km_ug_C <- subset(tc_4km_ug, treated == 0)
tc_4km_ug_T <- subset(tc_4km_ug, treated == 1)

pretrends_test_WI <- as.data.frame(matrix(NA,7,5))
pretrends_test_WI[1,1] <- 0
pretrends_test_WI[2,1] <- 1
pretrends_test_WI[3,1] <- 3
pretrends_test_WI[4,1] <- 5
pretrends_test_WI[5,1] <- 7.5
pretrends_test_WI[6,1] <- "5 (3km)"
pretrends_test_WI[7,1] <- "5 (4km)"

colnames(pretrends_test_WI) <- c("Penalty", "With UG 95%", "With UG 90%", 
                                 "Without UG 95%", "Without UG 90%")

###########
# zero penalty term pretrend test
###########

input_1 <- read.csv("preds_0p_a.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_0p_b.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_0p_c.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_0p_d.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_0p_e.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df

  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")

  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm


pretrends_test_WI[1,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) #6
pretrends_test_WI[1,3] <-length(which(pretrend_5q_CL1[,2] <= .1)) #7


###########
# one penalty term pretrend test
###########

input_1 <- read.csv("preds_1p_a.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_1p_b.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_1p_c.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_1p_d.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_1p_e.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[2,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) #5
pretrends_test_WI[2,3] <-length(which(pretrend_5q_CL1[,2] <= .1)) #13

###########
# three penalty term pretrend test
###########

input_1 <- read.csv("preds_3p_a.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_3p_b.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_3p_c.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_3p_d.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_3p_e.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[3,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) #4
pretrends_test_WI[3,3] <- length(which(pretrend_5q_CL1[,2] <= .1)) #12

###########
# five penalty term pretrend test
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
control_panel <- merge(tc_2km_ug_C, control_panel, by = "DHSID")

treatment_panel <- cbind(input_1, input_2, input_3, 
                         input_4, input_5)
treatment_panel <- treatment_panel[,c(1:2,4,6,8,10)]
treatment_panel$DHSID <- substr(treatment_panel$DHSID_year, 1, 14)
treatment_panel$year <- substr(treatment_panel$DHSID_year, 16, 19)
treatment_panel <- subset(treatment_panel, year < 2017)
treatment_panel <- subset(treatment_panel, year > 2005)
treatment_panel <- merge(tc_2km_ug_T, treatment_panel, by = "DHSID")

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[4,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) # 8
pretrends_test_WI[4,3] <- length(which(pretrend_5q_CL1[,2] <= .1)) # 13

###########
# seven point 5 penalty term pretrend test
###########

input_1 <- read.csv("preds_7p_a.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_7p_b.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_7p_c.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_7p_d.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_7p_e.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[5,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) #3
pretrends_test_WI[5,3] <- length(which(pretrend_5q_CL1[,2] <= .1)) #8


############
# without ug
###########

###########
# zero penalty term pretrend test
###########

input_1 <- read.csv("preds_0p_a_noUG.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_0p_b_noUG.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_0p_c_noUG.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_0p_d_noUG.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_0p_e_noUG.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[1,4] <- length(which(pretrend_5q_CL1[,2] <= .05)) #10
pretrends_test_WI[1,5] <- length(which(pretrend_5q_CL1[,2] <= .1)) #14

###########
# one penalty term pretrend test
###########

input_1 <- read.csv("preds_1p_a_noUG.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_1p_b_noUG.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_1p_c_noUG.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_1p_d_noUG.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_1p_e_noUG.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[2,4] <- length(which(pretrend_5q_CL1[,2] <= .05)) #8
pretrends_test_WI[2,5] <- length(which(pretrend_5q_CL1[,2] <= .1)) #11

###########
# three penalty term pretrend test
###########

input_1 <- read.csv("preds_3p_a_noUG.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_3p_b_noUG.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_3p_c_noUG.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_3p_d_noUG.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_3p_e_noUG.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[3,4] <- length(which(pretrend_5q_CL1[,2] <= .05)) #19
pretrends_test_WI[3,5] <- length(which(pretrend_5q_CL1[,2] <= .1)) #32

###########
# five penalty term pretrend test
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[4,4] <- length(which(pretrend_5q_CL1[,2] <= .05)) # 22
pretrends_test_WI[4,5] <- length(which(pretrend_5q_CL1[,2] <= .1)) # 29

###########
# seven point 5 penalty term pretrend test
###########

input_1 <- read.csv("preds_7p_a_noUG.csv")
colnames(input_1)[1] <- "DHSID_year"
input_2 <- read.csv("preds_7p_b_noUG.csv")
colnames(input_2)[1] <- "DHSID_year"
input_3 <- read.csv("preds_7p_c_noUG.csv")
colnames(input_3)[1] <- "DHSID_year"
input_4 <- read.csv("preds_7p_d_noUG.csv")
colnames(input_4)[1] <- "DHSID_year"
input_5 <- read.csv("preds_7p_e_noUG.csv")
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

ptm <- proc.time()

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[5,4] <- length(which(pretrend_5q_CL1[,2] <= .05)) #15
pretrends_test_WI[5,5] <- length(which(pretrend_5q_CL1[,2] <= .1)) #25


###########
# 3 and 4 km
###########

###########
# five penalty 3 km
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

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[6,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) # 6
pretrends_test_WI[6,3] <- length(which(pretrend_5q_CL1[,2] <= .1)) # 20

###########
# five penalty 4 km
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

set.seed(456) 

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
pretrend_5q_CL1 <- matrix(NA, 100, 2) # Column headers are at the bottom.

for(m in 1:n) {
  sample_C <- control_panel # import the control df
  randomized_C <- t(apply(sample_C[,7:11], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[13] <- "randomized"
  
  #same as above for the treatment
  sample_T <- treatment_panel 
  randomized_T <- t(apply(sample_T[,7:11], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[13] <- "randomized"
  
  unique_C_units <- unique(sample_C[,1]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,1]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  pretrend_5q_CL1[m,1] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  pretrend_5q_CL1[m,2] <- placebo_p[1,1]
}

proc.time() - ptm

pretrends_test_WI[7,2] <- length(which(pretrend_5q_CL1[,2] <= .05)) # 9
pretrends_test_WI[7,3] <- length(which(pretrend_5q_CL1[,2] <= .1)) # 18

pretrends_test_WI[6,4] <- "-"
pretrends_test_WI[6,5] <- "-"
pretrends_test_WI[7,4] <- "-"
pretrends_test_WI[7,5] <- "-"


write.csv(pretrends_test_WI, "pretrends_test_WI.csv", row.names = FALSE)





