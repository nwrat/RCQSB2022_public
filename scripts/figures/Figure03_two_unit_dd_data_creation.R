# THIS CREATES THE TWO UNIT ESTIMATE IN FIGURE 03

rm(list=ls())

library(fixest)

# Set your working directory
setwd("")

dhs_wi_ug <- read.csv("dhs_wi_ug.csv")
tc_2km <- read.csv("tc_2km_ug.csv")
tc_2km <- tc_2km[,c(2,5)]

tc_dhs <- merge( dhs_wi_ug, tc_2km, by = "DHSID")

treated_f <- function(x) { 
  if(x == 2006 | x == 2009) y <- 0
  if(x == 2011 | x == 2014 | x == 2016 ) y <- 1
  return(y)
}

tc_dhs$treat_year <- sapply(tc_dhs$year, treated_f)
tc_dhs$treat_id <- tc_dhs$treated * tc_dhs$treat_year

two_unit_results_df <- as.data.frame(matrix(NA,1,2))

two_unit_dd <- summary(feols(base_WI_no_elec_m ~ treat_id | year + treated, se = "iid", data=tc_dhs))
two_unit_results_df[1,1]  <- as.data.frame(coeftable(two_unit_dd)[1])
two_unit_results_df[1,2]  <- as.data.frame(coeftable(two_unit_dd)[2])

colnames(two_unit_results_df) <- c("dd_est", "dd_se")

write.csv(two_unit_results_df, "two_unit_results_df.csv", row.names = FALSE)








