# This is the code that creates the data for Extended Data Figure 11. 

rm(list=ls())

library(ggplot2)
library(dplyr)
library(fixest)
library(raster)
library(reshape2)

setwd("data/figure_and_input_data/")

# This is the prep code for ED Figure 11 a
ug_cell_ownership <- read.csv("ug_cell_ownership_06_16.csv")

cell_per_cluster <- ug_cell_ownership %>%
  group_by(year, DHSID) %>%
  summarise(year_ave = mean(hv243a, na.rm=TRUE)) 
cell_per_cluster <- as.data.frame(cell_per_cluster)

tc_2km <- read.csv("TC_2km_ug.csv")
tc_2km <- tc_2km[,-1]
tc_2km$year <- substr(tc_2km$DHSID, 3,6)

TC_cell <- merge(tc_2km, cell_per_cluster, by = c("year", "DHSID"))

TC_cell_ave_per_cluster <- TC_cell %>%
  group_by(year, treated) %>%
  summarise(cluster_ave = mean(year_ave, na.rm=TRUE)) 
TC_cell_ave_per_cluster <- as.data.frame(TC_cell_ave_per_cluster)

TC_cell_ave_per_cluster_C <- subset(TC_cell_ave_per_cluster, treated == 0)
TC_cell_ave_per_cluster_T <- subset(TC_cell_ave_per_cluster, treated == 1)
TC_cell_ave_per_cluster_TC <- rbind(TC_cell_ave_per_cluster_C, TC_cell_ave_per_cluster_T)

TC_cell_ave_per_cluster_TC$treat_year <- c(0,0,1,1,1,0,0,1,1,1)
TC_cell_ave_per_cluster_TC$treated_id <- TC_cell_ave_per_cluster_TC$treated * TC_cell_ave_per_cluster_TC$treat_year

# this DF is used to create ED Figure 11, a
write.csv(TC_cell_ave_per_cluster_TC, "ug_dhs_cellular_ownership.csv", row.names = FALSE)

# DiD and pretrends test

TC_cell_dd <- TC_cell

treated_year <- function(x) { 
  if(x == 2006 | x == 2009 ) y <- 0
  if(x == 2011 | x == 2014 | x == 2016) y <- 1
  return(y)
}

TC_cell_dd$treat_year <- sapply(TC_cell_dd$year, treated_year)
TC_cell_dd$treated_id <- TC_cell_dd$treated * TC_cell_dd$treat_year

summary(feols(year_ave ~ treated_id | year + treated, se = "iid", data=TC_cell_dd))

# pretrends
TC_cell_dd_pre <- subset(TC_cell_dd, year < 2011)

fake_treated_year <- function(x) { 
  if(x == 2006 ) y <- 0
  if(x == 2009) y <- 1
  return(y)
}

TC_cell_dd_pre$fake_treat_year <- sapply(TC_cell_dd_pre$year, fake_treated_year)
TC_cell_dd_pre$fake_treated_id <- TC_cell_dd_pre$treated * TC_cell_dd_pre$fake_treat_year

summary(feols(year_ave ~ fake_treated_id | year + treated, se = "iid", data=TC_cell_dd_pre))

##########
# This is the prep code for ED Figure 12 b and c
#########
tc_2km_ll <- tc_2km[,c(3,2)] # 

coverage_2010 <- raster("MCE_UG2G_2010.tif")
pr_coverage_2010 <- projectRaster(coverage_2010,crs="+proj=longlat +datum=NAD27")

coverage_2011 <- raster("MCE_UG2G_2011.tif")
pr_coverage_2011 <- projectRaster(coverage_2011,crs="+proj=longlat +datum=NAD27")

coverage_2012 <- raster("MCE_UG2G_2012.tif")
pr_coverage_2012 <- projectRaster(coverage_2012,crs="+proj=longlat +datum=NAD27")

coverage_2013 <- raster("MCE_UG2G_2013.tif")
pr_coverage_2013 <- projectRaster(coverage_2013,crs="+proj=longlat +datum=NAD27")

coverage_2014 <- raster("MCE_UG2G_2014.tif")
pr_coverage_2014 <- projectRaster(coverage_2014,crs="+proj=longlat +datum=NAD27")

coverage_2015 <- raster("MCE_UG2G_2015.tif")
pr_coverage_2015 <- projectRaster(coverage_2015,crs="+proj=longlat +datum=NAD27")

coverage_2016 <- raster("MCE_UG2G_2016.tif")
pr_coverage_2016 <- projectRaster(coverage_2016,crs="+proj=longlat +datum=NAD27")

g2_2010 <- raster::extract(pr_coverage_2010, tc_2km_ll)
g2_2010 <- as.data.frame(g2_2010)
g2_2011 <- raster::extract(pr_coverage_2011, tc_2km_ll)
g2_2011 <- as.data.frame(g2_2011)
g2_2012 <- raster::extract(pr_coverage_2012, tc_2km_ll)
g2_2012 <- as.data.frame(g2_2012)
g2_2013 <- raster::extract(pr_coverage_2013, tc_2km_ll)
g2_2013 <- as.data.frame(g2_2013)
g2_2014 <- raster::extract(pr_coverage_2014, tc_2km_ll)
g2_2014 <- as.data.frame(g2_2014)
g2_2015 <- raster::extract(pr_coverage_2015, tc_2km_ll)
g2_2015 <- as.data.frame(g2_2015)
g2_2016 <- raster::extract(pr_coverage_2016, tc_2km_ll)
g2_2016 <- as.data.frame(g2_2016)

g2_all <- cbind(tc_2km, g2_2010, g2_2011, g2_2012, g2_2013, g2_2014, g2_2015, g2_2016)

# Replace NAs (no coverage) with zeros
g2_all$g2_2010[is.na(g2_all$g2_2010)] <- 0
g2_all$g2_2011[is.na(g2_all$g2_2011)] <- 0
g2_all$g2_2012[is.na(g2_all$g2_2012)] <- 0
g2_all$g2_2013[is.na(g2_all$g2_2013)] <- 0
g2_all$g2_2014[is.na(g2_all$g2_2014)] <- 0
g2_all$g2_2015[is.na(g2_all$g2_2015)] <- 0
g2_all$g2_2016[is.na(g2_all$g2_2016)] <- 0

g2_all <- as.data.frame(g2_all)

# Make a binary indicator for cell service
some_cell <- function(x) { 
  if(x > 0) y <- 1
  if(x == 0) y <- 0
  return(y)
}
g2_all$some_2010 <- sapply(g2_all$g2_2010, some_cell)
g2_all$some_2011 <- sapply(g2_all$g2_2011, some_cell)
g2_all$some_2012 <- sapply(g2_all$g2_2012, some_cell)
g2_all$some_2013 <- sapply(g2_all$g2_2013, some_cell)
g2_all$some_2014 <- sapply(g2_all$g2_2014, some_cell)
g2_all$some_2015 <- sapply(g2_all$g2_2015, some_cell)
g2_all$some_2016 <- sapply(g2_all$g2_2016, some_cell)

# Variable cell quality is identified by "2", strong cell quality by "1"
# Recode so worse values are less than 1, making the figure more intuitive
cell_qual <- function(x) { 
  if(x > 1) y <- x/4
  if(x == 1) y <- 1
  if(x < 1) y <- 0
  return(y)
}

g2_all$recode_2010 <- sapply(g2_all$g2_2010, cell_qual)
g2_all$recode_2011 <- sapply(g2_all$g2_2011, cell_qual)
g2_all$recode_2012 <- sapply(g2_all$g2_2012, cell_qual)
g2_all$recode_2013 <- sapply(g2_all$g2_2013, cell_qual)
g2_all$recode_2014 <- sapply(g2_all$g2_2014, cell_qual)
g2_all$recode_2015 <- sapply(g2_all$g2_2015, cell_qual)
g2_all$recode_2016 <- sapply(g2_all$g2_2016, cell_qual)

# Getting yearly average by treatment group for quality
g2_TC_2010 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2010, na.rm=TRUE)) 
g2_TC_2010 <- as.data.frame(g2_TC_2010)

g2_TC_2011 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2011, na.rm=TRUE)) 
g2_TC_2011 <- as.data.frame(g2_TC_2011)

g2_TC_2012 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2012, na.rm=TRUE)) 
g2_TC_2012 <- as.data.frame(g2_TC_2012)

g2_TC_2013 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2013, na.rm=TRUE)) 
g2_TC_2013 <- as.data.frame(g2_TC_2013)

g2_TC_2014 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2014, na.rm=TRUE)) 
g2_TC_2014 <- as.data.frame(g2_TC_2014)

g2_TC_2015 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2015, na.rm=TRUE)) 
g2_TC_2015 <- as.data.frame(g2_TC_2015)

g2_TC_2016 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(recode_2016, na.rm=TRUE)) 
g2_TC_2016 <- as.data.frame(g2_TC_2016)

g2_TC_full <- rbind(g2_TC_2010, g2_TC_2011, g2_TC_2011, g2_TC_2013, g2_TC_2014, g2_TC_2015, g2_TC_2016)

g2_TC_full$year <- c("2010", "2010", "2011", "2011", "2012", "2012", "2013", "2013", "2014", "2014", 
                     "2015", "2015", "2016", "2016")

write.csv(g2_TC_full, "cell_quality.csv", row.names = FALSE)

# Getting yearly average by treatment group for coverage
g2_TCsome_2010 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2010, na.rm=TRUE)) 
g2_TCsome_2010 <- as.data.frame(g2_TCsome_2010)

g2_TCsome_2011 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2011, na.rm=TRUE)) 
g2_TCsome_2011 <- as.data.frame(g2_TCsome_2011)

g2_TCsome_2012 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2012, na.rm=TRUE)) 
g2_TCsome_2012 <- as.data.frame(g2_TCsome_2012)

g2_TCsome_2013 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2013, na.rm=TRUE)) 
g2_TCsome_2013 <- as.data.frame(g2_TCsome_2013)

g2_TCsome_2014 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2014, na.rm=TRUE)) 
g2_TCsome_2014 <- as.data.frame(g2_TCsome_2014)

g2_TCsome_2015 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2015, na.rm=TRUE)) 
g2_TCsome_2015 <- as.data.frame(g2_TCsome_2015)

g2_TCsome_2016 <- g2_all %>% 
  group_by(treated) %>%
  summarise(ave_coverage = mean(some_2016, na.rm=TRUE)) 
g2_TCsome_2016 <- as.data.frame(g2_TCsome_2016)

g2_TCsome_full <- rbind(g2_TCsome_2010, g2_TCsome_2011, g2_TCsome_2011, g2_TCsome_2013, g2_TCsome_2014, g2_TCsome_2015, g2_TCsome_2016)

g2_TCsome_full$year <- c("2010", "2010", "2011", "2011", "2012", "2012", "2013", "2013", "2014", "2014", 
                         "2015", "2015", "2016", "2016")

write.csv(g2_TCsome_full, "cell_coverage.csv", row.names = FALSE)

###########
# combining a full data set for DiD #
###########
g2_full_long_quality <- g2_all[,c(1,4,20:26)]
g2_full_long_quality <- melt(g2_full_long_quality, id.vars = c("DHSID", "treated"))
g2_full_long_quality$year <- rep(c(2010:2016), each = 964) 

treated_year <- function(x) { 
  if(x == 2010 ) y <- 0
  if(x == 2011 | x == 2012 | x == 2013 | x == 2014 | x == 2015 | x == 2016) y <- 1
  return(y)
}

g2_full_long_quality$treat_year <- sapply(g2_full_long_quality$year, treated_year)
g2_full_long_quality$treated_id <- g2_full_long_quality$treated * g2_full_long_quality$treat_year

summary(feols(value ~ treated_id | DHSID + year, data=g2_full_long_quality))

# and now for coverage
g2_full_long_coverage <- g2_all[,c(1,4,13:19)]
g2_full_long_coverage <- melt(g2_full_long_coverage, id.vars = c("DHSID", "treated"))
g2_full_long_coverage$year <- rep(c(2010:2016), each = 964) 

treated_year <- function(x) { 
  if(x == 2010 ) y <- 0
  if(x == 2011 | x == 2012 | x == 2013 | x == 2014 | x == 2015 | x == 2016) y <- 1
  return(y)
}

g2_full_long_coverage$treat_year <- sapply(g2_full_long_coverage$year, treated_year)
g2_full_long_coverage$treated_id <- g2_full_long_coverage$treated * g2_full_long_coverage$treat_year

summary(feols(value ~ treated_id | DHSID + year, data=g2_full_long_coverage))






