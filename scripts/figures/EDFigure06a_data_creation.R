# This code creates the DF for the DiD and MC simulation in Extended Data Figure 6a
# THIS CODE CREATES THE DF FOR THE DD / MC SIMULATION IN ED FIGURE 6a

library(MCPanel)
library(dplyr)
library(lfe)
library(reshape2)
library(ggthemes)

setwd("data/figure_and_input_data/")

###########
# .1 alpha
##########
ptm <- proc.time()
set.seed(789)

r = 10 #number of different years, i do all even years 4-22. 
dd_mc_loop_output <- matrix(NA,10,9) #this capture all the final estimates
for(s in 1:r) { 
  v=50 
  inside_loop_biased <- matrix(NA,3,v) 
  for(i in 1:v) {
    nc = 100 #number of units
    ny = 4+((s-1)*2) # sets number of years in DF   
    beta = 1  #treatment effect
    sdn = 5  
    cty <- data.frame(cty = 1:nc, alpha = rnorm(nc,0,0.1), delta=runif(nc,0,10)) # this creates entity values, alpha is basically the growth rate. 
    ll <- cty$alpha>0  # captures randomly selected positive growing places
    cty$treat <- 0
    cty$treat[ll] <- round(runif(sum(ll),((ny/2)+1),ny+(.5*ny))) # this selects only positive places to be treated, and controls the number of places that can be treated
    dta <- data.frame(cty = rep(1:nc,each=ny), year = rep(1:ny,times=nc), e = rnorm(ny*nc,0,1))  #this generates country-year panel; e is random noise in outcome
    dta <- left_join(dta,cty,by="cty")
    dta$treatdummy <- 1*(dta$year >= dta$treat)  #creates actual treatment dummy in panel
    d = nc*ny   
    for(c in 1:d){
      if(dta[c,6] == 0){
        dta[c,7] <- 0
      }
    }
    dta$ay <- dta$alpha*dta$year  #change in outcome due to trend
    dta$y <- dta$delta + dta$ay + dta$treatdummy*beta + dta$e  #observed y:  base + trend + treatment effect + noise
    did_biased <- felm(y ~ treatdummy | cty + year, data=dta) #standard pooled DD
    inside_loop_biased[1,i] <- summary(did_biased)$coefficients[1,1]
    
    #non-pooled did
    dta_nonpooled <- dta
    dta_nonpooled <- subset(dta_nonpooled, year <= (ny/2) | year == ny)
    did_nonpooled <- felm(y ~ treatdummy | cty + year, data=dta_nonpooled) #non-pooled DD
    inside_loop_biased[2,i] <- summary(did_nonpooled)$coefficients[1,1]
    
    # matrix completion
    mc_panel_biased_long <- matrix(NA,nc*ny,3)
    mc_panel_biased_long[,1] <- dta[,1]
    mc_panel_biased_long[,2] <- dta[,2]
    mc_panel_biased_long[,3] <- dta[,9]
    mc_panel_biased_long <- as.data.frame(mc_panel_biased_long)
    names(mc_panel_biased_long) <- c("cty", "year", "y")
    mc_panel_biased <- dcast(mc_panel_biased_long, cty ~ year, value.var="y")
    
    mc_treatpanel_biased_long <- matrix(NA,nc*ny,3)
    mc_treatpanel_biased_long[,1] <- dta[,1]
    mc_treatpanel_biased_long[,2] <- dta[,2]
    mc_treatpanel_biased_long[,3] <- dta[,7]
    mc_treatpanel_biased_long <- as.data.frame(mc_treatpanel_biased_long)
    names(mc_treatpanel_biased_long) <- c("cty", "year", "treated")
    mc_treatpanel_biased <- dcast(mc_treatpanel_biased_long, cty ~ year, value.var="treated")
    mc_treatpanel_biased[1:100,2:(ny+1)] <- 1-mc_treatpanel_biased[1:100,2:(ny+1)] #switching zeros and ones for mc input
    
    N_basic_control_test <- 100
    T_basic_control_test <- ny
    M_basic_control_test <- t(as.matrix(mc_panel_biased[1:100,2:(ny+1)]))
    
    mask_basic_control_test <- t(as.matrix(mc_treatpanel_biased[1:100,2:(ny+1)]))
    
    model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                   to_estimate_v = 1)
    
    model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
      replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
    
    cf <- model_with_both_basic_control_test$est #counterfactuals
    cf <- as.matrix(cf)
    
    # creating a counterfactual panel to only capture treated units. 
    cf_treated <- matrix(0L,ny,100)
    x=100
    y=ny
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          cf_treated[j,k] <- cf[j,k]
        } 
      }
    }
    
    # creating an observed panel to only capture treated units. 
    obs_treated <- matrix(0L,ny,100)
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          obs_treated[j,k] <- M_basic_control_test[j,k]
        } 
      }
    }
    
    # ATE per unit
    ate_calc <- matrix(0L,ny,100)
    ate_calc[1:ny,1:100] <- obs_treated[1:ny,1:100] - cf_treated[1:ny,1:100]
    
    # mean ATE over all units in each loop
    z=ny
    ATE_ave <- matrix(0L,ny,1)
    for(l in 1:z){
      ATE_ave <- (sum(ate_calc[ny,])) / (100 - length(which(ate_calc[ny,] == 0)))
    }
    
    # ATE in ith time period, i.e. the last time period of the current loop (so years 4, 6, 8... 22) 
    inside_loop_biased[3,i] <- ATE_ave
  }
  
  dd_mc_loop_output[s,1] <- mean(inside_loop_biased[1,]) 
  dd_mc_loop_output[s,2] <- sd(inside_loop_biased[1,]) # DID SD
  dd_mc_loop_output[s,3] <- abs((mean(inside_loop_biased[1,])-1))
  
  dd_mc_loop_output[s,4] <- mean(inside_loop_biased[3,]) # MC ATE
  dd_mc_loop_output[s,5] <- sd(inside_loop_biased[3,]) # MC SD
  dd_mc_loop_output[s,6] <- abs((mean(inside_loop_biased[3,])-1))
  
  dd_mc_loop_output[s,7] <- mean(inside_loop_biased[2,]) #
  dd_mc_loop_output[s,8] <- sd(inside_loop_biased[2,]) # 
  dd_mc_loop_output[s,9] <- abs((mean(inside_loop_biased[2,])-1))
}
proc.time() - ptm  

dd_mc_loop_output <- as.data.frame(dd_mc_loop_output)
year_col <- matrix(NA,10,1)
year_col[,1] <- seq(4,22,by=2)
dd_mc_loop_output_graph <- cbind(year_col,dd_mc_loop_output)
names(dd_mc_loop_output_graph) <- c("Year", "DID_est", "DID_sd", "DID_diff", "MC_est", "MC_sd", "MC_diff", "DD_np_est", "DD_np_sd", "DD_np_diff")

write.csv(dd_mc_loop_output_graph, "dd_mc_output_1alpha.csv")

###########
# .15 alpha
##########
ptm <- proc.time()
set.seed(789)

r = 10 #number of different years, i do all even years 4-22. 
dd_mc_loop_output <- matrix(NA,10,9) #this capture all the final estimates
for(s in 1:r) { 
  v=50 
  inside_loop_biased <- matrix(NA,3,v) 
  for(i in 1:v) {
    nc = 100 #number of units
    ny = 4+((s-1)*2) # sets number of years in DF   
    beta = 1  #treatment effect
    sdn = 5  
    cty <- data.frame(cty = 1:nc, alpha = rnorm(nc,0,0.15), delta=runif(nc,0,10)) # this creates entity values, alpha is basically the growth rate. 
    ll <- cty$alpha>0  # captures randomly selected positive growing places
    cty$treat <- 0
    cty$treat[ll] <- round(runif(sum(ll),((ny/2)+1),ny+(.5*ny))) # this selects only positive places to be treated, and controls the number of places that can be treated
    dta <- data.frame(cty = rep(1:nc,each=ny), year = rep(1:ny,times=nc), e = rnorm(ny*nc,0,1))  #this generates country-year panel; e is random noise in outcome
    dta <- left_join(dta,cty,by="cty")
    dta$treatdummy <- 1*(dta$year >= dta$treat)  #creates actual treatment dummy in panel
    d = nc*ny   
    for(c in 1:d){
      if(dta[c,6] == 0){
        dta[c,7] <- 0
      }
    }
    dta$ay <- dta$alpha*dta$year  #change in outcome due to trend
    dta$y <- dta$delta + dta$ay + dta$treatdummy*beta + dta$e  #observed y:  base + trend + treatment effect + noise
    did_biased <- felm(y ~ treatdummy | cty + year, data=dta) #standard pooled DD
    inside_loop_biased[1,i] <- summary(did_biased)$coefficients[1,1]
    
    #non-pooled did
    dta_nonpooled <- dta
    dta_nonpooled <- subset(dta_nonpooled, year <= (ny/2) | year == ny)
    did_nonpooled <- felm(y ~ treatdummy | cty + year, data=dta_nonpooled) #non-pooled DD
    inside_loop_biased[2,i] <- summary(did_nonpooled)$coefficients[1,1]
    
    # matrix completion
    mc_panel_biased_long <- matrix(NA,nc*ny,3)
    mc_panel_biased_long[,1] <- dta[,1]
    mc_panel_biased_long[,2] <- dta[,2]
    mc_panel_biased_long[,3] <- dta[,9]
    mc_panel_biased_long <- as.data.frame(mc_panel_biased_long)
    names(mc_panel_biased_long) <- c("cty", "year", "y")
    mc_panel_biased <- dcast(mc_panel_biased_long, cty ~ year, value.var="y")
    
    mc_treatpanel_biased_long <- matrix(NA,nc*ny,3)
    mc_treatpanel_biased_long[,1] <- dta[,1]
    mc_treatpanel_biased_long[,2] <- dta[,2]
    mc_treatpanel_biased_long[,3] <- dta[,7]
    mc_treatpanel_biased_long <- as.data.frame(mc_treatpanel_biased_long)
    names(mc_treatpanel_biased_long) <- c("cty", "year", "treated")
    mc_treatpanel_biased <- dcast(mc_treatpanel_biased_long, cty ~ year, value.var="treated")
    mc_treatpanel_biased[1:100,2:(ny+1)] <- 1-mc_treatpanel_biased[1:100,2:(ny+1)] #switching zeros and ones for mc input
    
    N_basic_control_test <- 100
    T_basic_control_test <- ny
    M_basic_control_test <- t(as.matrix(mc_panel_biased[1:100,2:(ny+1)]))
    
    mask_basic_control_test <- t(as.matrix(mc_treatpanel_biased[1:100,2:(ny+1)]))
    
    model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                   to_estimate_v = 1)
    
    model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
      replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
    
    cf <- model_with_both_basic_control_test$est #counterfactuals
    cf <- as.matrix(cf)
    
    # creating a counterfactual panel to only capture treated units. 
    cf_treated <- matrix(0L,ny,100)
    x=100
    y=ny
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          cf_treated[j,k] <- cf[j,k]
        } 
      }
    }
    
    # creating an observed panel to only capture treated units. 
    obs_treated <- matrix(0L,ny,100)
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          obs_treated[j,k] <- M_basic_control_test[j,k]
        } 
      }
    }
    
    # ATE per unit
    ate_calc <- matrix(0L,ny,100)
    ate_calc[1:ny,1:100] <- obs_treated[1:ny,1:100] - cf_treated[1:ny,1:100]
    
    # mean ATE over all units in each loop
    z=ny
    ATE_ave <- matrix(0L,ny,1)
    for(l in 1:z){
      ATE_ave <- (sum(ate_calc[ny,])) / (100 - length(which(ate_calc[ny,] == 0)))
    }
    
    # ATE in ith time period, i.e. the last time period of the current loop (so years 4, 6, 8... 22) 
    inside_loop_biased[3,i] <- ATE_ave
  }
  
  dd_mc_loop_output[s,1] <- mean(inside_loop_biased[1,]) 
  dd_mc_loop_output[s,2] <- sd(inside_loop_biased[1,]) # DID SD
  dd_mc_loop_output[s,3] <- abs((mean(inside_loop_biased[1,])-1))
  
  dd_mc_loop_output[s,4] <- mean(inside_loop_biased[3,]) # MC ATE
  dd_mc_loop_output[s,5] <- sd(inside_loop_biased[3,]) # MC SD
  dd_mc_loop_output[s,6] <- abs((mean(inside_loop_biased[3,])-1))
  
  dd_mc_loop_output[s,7] <- mean(inside_loop_biased[2,]) #
  dd_mc_loop_output[s,8] <- sd(inside_loop_biased[2,]) # 
  dd_mc_loop_output[s,9] <- abs((mean(inside_loop_biased[2,])-1))
}
proc.time() - ptm  

dd_mc_loop_output <- as.data.frame(dd_mc_loop_output)
year_col <- matrix(NA,10,1)
year_col[,1] <- seq(4,22,by=2)
dd_mc_loop_output_graph <- cbind(year_col,dd_mc_loop_output)
names(dd_mc_loop_output_graph) <- c("Year", "DID_est", "DID_sd", "DID_diff", "MC_est", "MC_sd", "MC_diff", "DD_np_est", "DD_np_sd", "DD_np_diff")

write.csv(dd_mc_loop_output_graph, "dd_mc_output_15alpha.csv")

###########
# .2 alpha
##########
ptm <- proc.time()
set.seed(789)

r = 10 #number of different years, i do all even years 4-22. 
dd_mc_loop_output <- matrix(NA,10,9) #this capture all the final estimates
for(s in 1:r) { 
  v=50 
  inside_loop_biased <- matrix(NA,3,v) 
  for(i in 1:v) {
    nc = 100 #number of units
    ny = 4+((s-1)*2) # sets number of years in DF   
    beta = 1  #treatment effect
    sdn = 5  
    cty <- data.frame(cty = 1:nc, alpha = rnorm(nc,0,0.2), delta=runif(nc,0,10)) # this creates entity values, alpha is basically the growth rate. 
    ll <- cty$alpha>0  # captures randomly selected positive growing places
    cty$treat <- 0
    cty$treat[ll] <- round(runif(sum(ll),((ny/2)+1),ny+(.5*ny))) # this selects only positive places to be treated, and controls the number of places that can be treated
    dta <- data.frame(cty = rep(1:nc,each=ny), year = rep(1:ny,times=nc), e = rnorm(ny*nc,0,1))  #this generates country-year panel; e is random noise in outcome
    dta <- left_join(dta,cty,by="cty")
    dta$treatdummy <- 1*(dta$year >= dta$treat)  #creates actual treatment dummy in panel
    d = nc*ny   
    for(c in 1:d){
      if(dta[c,6] == 0){
        dta[c,7] <- 0
      }
    }
    dta$ay <- dta$alpha*dta$year  #change in outcome due to trend
    dta$y <- dta$delta + dta$ay + dta$treatdummy*beta + dta$e  #observed y:  base + trend + treatment effect + noise
    did_biased <- felm(y ~ treatdummy | cty + year, data=dta) #standard pooled DD
    inside_loop_biased[1,i] <- summary(did_biased)$coefficients[1,1]
    
    #non-pooled did
    dta_nonpooled <- dta
    dta_nonpooled <- subset(dta_nonpooled, year <= (ny/2) | year == ny)
    did_nonpooled <- felm(y ~ treatdummy | cty + year, data=dta_nonpooled) #non-pooled DD
    inside_loop_biased[2,i] <- summary(did_nonpooled)$coefficients[1,1]
    
    # matrix completion
    mc_panel_biased_long <- matrix(NA,nc*ny,3)
    mc_panel_biased_long[,1] <- dta[,1]
    mc_panel_biased_long[,2] <- dta[,2]
    mc_panel_biased_long[,3] <- dta[,9]
    mc_panel_biased_long <- as.data.frame(mc_panel_biased_long)
    names(mc_panel_biased_long) <- c("cty", "year", "y")
    mc_panel_biased <- dcast(mc_panel_biased_long, cty ~ year, value.var="y")
    
    mc_treatpanel_biased_long <- matrix(NA,nc*ny,3)
    mc_treatpanel_biased_long[,1] <- dta[,1]
    mc_treatpanel_biased_long[,2] <- dta[,2]
    mc_treatpanel_biased_long[,3] <- dta[,7]
    mc_treatpanel_biased_long <- as.data.frame(mc_treatpanel_biased_long)
    names(mc_treatpanel_biased_long) <- c("cty", "year", "treated")
    mc_treatpanel_biased <- dcast(mc_treatpanel_biased_long, cty ~ year, value.var="treated")
    mc_treatpanel_biased[1:100,2:(ny+1)] <- 1-mc_treatpanel_biased[1:100,2:(ny+1)] #switching zeros and ones for mc input
    
    N_basic_control_test <- 100
    T_basic_control_test <- ny
    M_basic_control_test <- t(as.matrix(mc_panel_biased[1:100,2:(ny+1)]))
    
    mask_basic_control_test <- t(as.matrix(mc_treatpanel_biased[1:100,2:(ny+1)]))
    
    model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                   to_estimate_v = 1)
    
    model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
      replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
    
    cf <- model_with_both_basic_control_test$est #counterfactuals
    cf <- as.matrix(cf)
    
    # creating a counterfactual panel to only capture treated units. 
    cf_treated <- matrix(0L,ny,100)
    x=100
    y=ny
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          cf_treated[j,k] <- cf[j,k]
        } 
      }
    }
    
    # creating an observed panel to only capture treated units. 
    obs_treated <- matrix(0L,ny,100)
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          obs_treated[j,k] <- M_basic_control_test[j,k]
        } 
      }
    }
    
    # ATE per unit
    ate_calc <- matrix(0L,ny,100)
    ate_calc[1:ny,1:100] <- obs_treated[1:ny,1:100] - cf_treated[1:ny,1:100]
    
    # mean ATE over all units in each loop
    z=ny
    ATE_ave <- matrix(0L,ny,1)
    for(l in 1:z){
      ATE_ave <- (sum(ate_calc[ny,])) / (100 - length(which(ate_calc[ny,] == 0)))
    }
    
    # ATE in ith time period, i.e. the last time period of the current loop (so years 4, 6, 8... 22) 
    inside_loop_biased[3,i] <- ATE_ave
  }
  
  dd_mc_loop_output[s,1] <- mean(inside_loop_biased[1,]) 
  dd_mc_loop_output[s,2] <- sd(inside_loop_biased[1,]) # DID SD
  dd_mc_loop_output[s,3] <- abs((mean(inside_loop_biased[1,])-1))
  
  dd_mc_loop_output[s,4] <- mean(inside_loop_biased[3,]) # MC ATE
  dd_mc_loop_output[s,5] <- sd(inside_loop_biased[3,]) # MC SD
  dd_mc_loop_output[s,6] <- abs((mean(inside_loop_biased[3,])-1))
  
  dd_mc_loop_output[s,7] <- mean(inside_loop_biased[2,]) #
  dd_mc_loop_output[s,8] <- sd(inside_loop_biased[2,]) # 
  dd_mc_loop_output[s,9] <- abs((mean(inside_loop_biased[2,])-1))
}
proc.time() - ptm  

dd_mc_loop_output <- as.data.frame(dd_mc_loop_output)
year_col <- matrix(NA,10,1)
year_col[,1] <- seq(4,22,by=2)
dd_mc_loop_output_graph <- cbind(year_col,dd_mc_loop_output)
names(dd_mc_loop_output_graph) <- c("Year", "DID_est", "DID_sd", "DID_diff", "MC_est", "MC_sd", "MC_diff", "DD_np_est", "DD_np_sd", "DD_np_diff")

write.csv(dd_mc_loop_output_graph, "dd_mc_output_2alpha.csv")


###########
# .25 alpha
##########
ptm <- proc.time()
set.seed(789)

r = 10 #number of different years, i do all even years 4-22. 
dd_mc_loop_output <- matrix(NA,10,9) #this capture all the final estimates
for(s in 1:r) { 
  v=50 
  inside_loop_biased <- matrix(NA,3,v) 
  for(i in 1:v) {
    nc = 100 #number of units
    ny = 4+((s-1)*2) # sets number of years in DF   
    beta = 1  #treatment effect
    sdn = 5  
    cty <- data.frame(cty = 1:nc, alpha = rnorm(nc,0,0.25), delta=runif(nc,0,10)) # this creates entity values, alpha is basically the growth rate. 
    ll <- cty$alpha>0  # captures randomly selected positive growing places
    cty$treat <- 0
    cty$treat[ll] <- round(runif(sum(ll),((ny/2)+1),ny+(.5*ny))) # this selects only positive places to be treated, and controls the number of places that can be treated
    dta <- data.frame(cty = rep(1:nc,each=ny), year = rep(1:ny,times=nc), e = rnorm(ny*nc,0,1))  #this generates country-year panel; e is random noise in outcome
    dta <- left_join(dta,cty,by="cty")
    dta$treatdummy <- 1*(dta$year >= dta$treat)  #creates actual treatment dummy in panel
    d = nc*ny   
    for(c in 1:d){
      if(dta[c,6] == 0){
        dta[c,7] <- 0
      }
    }
    dta$ay <- dta$alpha*dta$year  #change in outcome due to trend
    dta$y <- dta$delta + dta$ay + dta$treatdummy*beta + dta$e  #observed y:  base + trend + treatment effect + noise
    did_biased <- felm(y ~ treatdummy | cty + year, data=dta) #standard pooled DD
    inside_loop_biased[1,i] <- summary(did_biased)$coefficients[1,1]
    
    #non-pooled did
    dta_nonpooled <- dta
    dta_nonpooled <- subset(dta_nonpooled, year <= (ny/2) | year == ny)
    did_nonpooled <- felm(y ~ treatdummy | cty + year, data=dta_nonpooled) #non-pooled DD
    inside_loop_biased[2,i] <- summary(did_nonpooled)$coefficients[1,1]
    
    # matrix completion
    mc_panel_biased_long <- matrix(NA,nc*ny,3)
    mc_panel_biased_long[,1] <- dta[,1]
    mc_panel_biased_long[,2] <- dta[,2]
    mc_panel_biased_long[,3] <- dta[,9]
    mc_panel_biased_long <- as.data.frame(mc_panel_biased_long)
    names(mc_panel_biased_long) <- c("cty", "year", "y")
    mc_panel_biased <- dcast(mc_panel_biased_long, cty ~ year, value.var="y")
    
    mc_treatpanel_biased_long <- matrix(NA,nc*ny,3)
    mc_treatpanel_biased_long[,1] <- dta[,1]
    mc_treatpanel_biased_long[,2] <- dta[,2]
    mc_treatpanel_biased_long[,3] <- dta[,7]
    mc_treatpanel_biased_long <- as.data.frame(mc_treatpanel_biased_long)
    names(mc_treatpanel_biased_long) <- c("cty", "year", "treated")
    mc_treatpanel_biased <- dcast(mc_treatpanel_biased_long, cty ~ year, value.var="treated")
    mc_treatpanel_biased[1:100,2:(ny+1)] <- 1-mc_treatpanel_biased[1:100,2:(ny+1)] #switching zeros and ones for mc input
    
    N_basic_control_test <- 100
    T_basic_control_test <- ny
    M_basic_control_test <- t(as.matrix(mc_panel_biased[1:100,2:(ny+1)]))
    
    mask_basic_control_test <- t(as.matrix(mc_treatpanel_biased[1:100,2:(ny+1)]))
    
    model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                   to_estimate_v = 1)
    
    model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
      replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
    
    cf <- model_with_both_basic_control_test$est #counterfactuals
    cf <- as.matrix(cf)
    
    # creating a counterfactual panel to only capture treated units. 
    cf_treated <- matrix(0L,ny,100)
    x=100
    y=ny
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          cf_treated[j,k] <- cf[j,k]
        } 
      }
    }
    
    # creating an observed panel to only capture treated units. 
    obs_treated <- matrix(0L,ny,100)
    for(j in 1:y){
      for(k in 1:x){
        if(mask_basic_control_test[j,k]==0){
          obs_treated[j,k] <- M_basic_control_test[j,k]
        } 
      }
    }
    
    # ATE per unit
    ate_calc <- matrix(0L,ny,100)
    ate_calc[1:ny,1:100] <- obs_treated[1:ny,1:100] - cf_treated[1:ny,1:100]
    
    # mean ATE over all units in each loop
    z=ny
    ATE_ave <- matrix(0L,ny,1)
    for(l in 1:z){
      ATE_ave <- (sum(ate_calc[ny,])) / (100 - length(which(ate_calc[ny,] == 0)))
    }
    
    # ATE in ith time period, i.e. the last time period of the current loop (so years 4, 6, 8... 22) 
    inside_loop_biased[3,i] <- ATE_ave
  }
  
  dd_mc_loop_output[s,1] <- mean(inside_loop_biased[1,]) 
  dd_mc_loop_output[s,2] <- sd(inside_loop_biased[1,]) # DID SD
  dd_mc_loop_output[s,3] <- abs((mean(inside_loop_biased[1,])-1))
  
  dd_mc_loop_output[s,4] <- mean(inside_loop_biased[3,]) # MC ATE
  dd_mc_loop_output[s,5] <- sd(inside_loop_biased[3,]) # MC SD
  dd_mc_loop_output[s,6] <- abs((mean(inside_loop_biased[3,])-1))
  
  dd_mc_loop_output[s,7] <- mean(inside_loop_biased[2,]) #
  dd_mc_loop_output[s,8] <- sd(inside_loop_biased[2,]) # 
  dd_mc_loop_output[s,9] <- abs((mean(inside_loop_biased[2,])-1))
}
proc.time() - ptm  

dd_mc_loop_output <- as.data.frame(dd_mc_loop_output)
year_col <- matrix(NA,10,1)
year_col[,1] <- seq(4,22,by=2)
dd_mc_loop_output_graph <- cbind(year_col,dd_mc_loop_output)
names(dd_mc_loop_output_graph) <- c("Year", "DID_est", "DID_sd", "DID_diff", "MC_est", "MC_sd", "MC_diff", "DD_np_est", "DD_np_sd", "DD_np_diff")

write.csv(dd_mc_loop_output_graph, "dd_mc_output_25alpha.csv")
