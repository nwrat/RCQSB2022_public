# This creates the inverse distance weighting estimate in Figure 3. 

library(ggplot2)
library(devtools)
library(sf)
library(fixest)
library(units)

setwd()  

atlas_cluster_data_w_WI_scenarios <- read.csv("atlas_cluster_data_nov_12_2020.csv", stringsAsFactors = FALSE) 
atlas_cluster_data_w_WI_scenarios <- atlas_cluster_data_w_WI_scenarios[,-1]
as.numeric(atlas_cluster_data_w_WI_scenarios$DHSID)
as.numeric(atlas_cluster_data_w_WI_scenarios$year)
atlas_cluster_data_w_WI_scenarios$survey_year <- substr(atlas_cluster_data_w_WI_scenarios$DHSID, 3, 6)
atlas_cluster_data_w_WI_scenarios <- atlas_cluster_data_w_WI_scenarios[,c(1:2, 20:22, 27)]
atlas_cluster_data_w_WI_scenarios$country <- substr(atlas_cluster_data_w_WI_scenarios$DHSID, 1, 2)
atlas_cluster_data_w_WI_scenarios <- subset(atlas_cluster_data_w_WI_scenarios, country == "UG")

ug_panel_2006_all <- subset(atlas_cluster_data_w_WI_scenarios, year == 2006)
ug_panel_2009_all <- subset(atlas_cluster_data_w_WI_scenarios, year == 2009)
ug_panel_2011_all <- subset(atlas_cluster_data_w_WI_scenarios, year == 2011)
ug_panel_2014_all <- subset(atlas_cluster_data_w_WI_scenarios, year == 2014)
ug_panel_2016_all <- subset(atlas_cluster_data_w_WI_scenarios, year == 2016)

ug_2006_rural_sf <- ug_panel_2006_all %>% 
  dplyr::select("long", "lat") %>%
  st_as_sf(coords = c("long", "lat"), 
           crs    = "+init=epsg:4326") %>%
  st_transform("+init=epsg:29902") %>%
  st_buffer(dist = 5000) %>% 
  st_transform("+init=epsg:4326")

ug_2009_rural_sf <- ug_panel_2009_all %>% 
  dplyr::select("long", "lat") %>%
  st_as_sf(coords = c("long", "lat"), 
           crs    = "+init=epsg:4326") %>%
  st_transform("+init=epsg:29902") %>%
  st_buffer(dist = 5000) %>% 
  st_transform("+init=epsg:4326")

ug_2011_rural_sf <- ug_panel_2011_all %>% 
  dplyr::select("long", "lat") %>%
  st_as_sf(coords = c("long", "lat"), 
           crs    = "+init=epsg:4326") %>%
  st_transform("+init=epsg:29902") %>%
  st_buffer(dist = 5000) %>%
  st_transform("+init=epsg:4326")

ug_2014_rural_sf <- ug_panel_2014_all %>% 
  dplyr::select("long", "lat") %>%
  st_as_sf(coords = c("long", "lat"), 
           crs    = "+init=epsg:4326") %>%
  st_transform("+init=epsg:29902") %>%
  st_buffer(dist = 5000) %>% 
  st_transform("+init=epsg:4326")

ug_2016_rural_sf <- ug_panel_2016_all %>% 
  dplyr::select("long", "lat") %>%
  st_as_sf(coords = c("long", "lat"), 
           crs    = "+init=epsg:4326") %>%
  st_transform("+init=epsg:29902") %>%
  st_buffer(dist = 5000) %>% 
  st_transform("+init=epsg:4326")

#Loading distances between units. 
distance_ug_06_06 <- read.csv("distance_ug_06_06_dec.csv")
distance_ug_06_06 <- distance_ug_06_06[,-1]

distance_ug_09_09 <- read.csv("distance_ug_09_09_dec.csv")
distance_ug_09_09 <- distance_ug_09_09[,-1]

distance_ug_11_11 <- read.csv("distance_ug_11_11_dec.csv")
distance_ug_11_11 <- distance_ug_11_11[,-1]

distance_ug_14_14 <- read.csv("distance_ug_14_14_dec.csv")
distance_ug_14_14 <- distance_ug_14_14[,-1]

distance_ug_16_16 <- read.csv("distance_ug_16_16_dec.csv")
distance_ug_16_16 <- distance_ug_16_16[,-1]

distance_ug_16_14 <- read.csv("distance_ug_16_14_dec.csv")
distance_ug_16_14 <- distance_ug_16_14[,-1]

distance_ug_16_11 <- read.csv("distance_ug_16_11_dec.csv")
distance_ug_16_11 <- distance_ug_16_11[,-1]

distance_ug_16_09 <- read.csv("distance_ug_16_09_dec.csv")
distance_ug_16_09 <- distance_ug_16_09[,-1]

distance_ug_16_06 <- read.csv("distance_ug_16_06_dec.csv")
distance_ug_16_06 <- distance_ug_16_06[,-1]

distance_ug_14_11 <- read.csv("distance_ug_14_11_dec.csv")
distance_ug_14_11 <- distance_ug_14_11[,-1]

distance_ug_14_09 <- read.csv("distance_ug_14_09_dec.csv")
distance_ug_14_09 <- distance_ug_14_09[,-1]

distance_ug_14_06 <- read.csv("distance_ug_14_06_dec.csv")
distance_ug_14_06 <- distance_ug_14_06[,-1]

distance_ug_11_09 <- read.csv("distance_ug_11_09_dec.csv")
distance_ug_11_09 <- distance_ug_11_09[,-1]

distance_ug_11_06 <- read.csv("distance_ug_11_06_dec.csv")
distance_ug_11_06 <- distance_ug_11_06[,-1]

distance_ug_09_06 <- read.csv("distance_ug_09_06_dec.csv")
distance_ug_09_06 <- distance_ug_09_06[,-1]

###########
# starting IDW estimates
###########

###########
# 2016 in 2006
###########
distance_ug_16_06 <- drop_units(distance_ug_16_06)

# loop to get the distance weighted WI for each cell
r <- 684
c <- 336
WI_weighted_matrix_16_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_16_in_06_10k[i,j] <- if(distance_ug_16_06[i,j] < 10000) {
      (ug_panel_2006_all[j,5] / (distance_ug_16_06[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_16_in_06_10k) <- sapply(WI_weighted_matrix_16_in_06_10k, is.infinite)

# loop to complete the 1/distance component (denominator of the IDW function)
r <- 684
c <- 336
WI_weighted_distance_16_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_16_in_06_10k[i,j] <- if(distance_ug_16_06[i,j] < 10000) {
      (1 / distance_ug_16_06[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_16_in_06_10k) <- sapply(WI_weighted_distance_16_in_06_10k, is.infinite)

# then sum across each row
ave_WI_weighted_cluster_16_in_06_10k <- matrix(NA,684,1)
ave_WI_weighted_cluster_16_in_06_10k[,1] <- rowSums(WI_weighted_matrix_16_in_06_10k[,1:336], na.rm = TRUE) / rowSums(WI_weighted_distance_16_in_06_10k[,1:336], na.rm = TRUE)

##########
# 2016 in 2009
##########
distance_ug_16_09 <- drop_units(distance_ug_16_09)

r <- 684
c <- 170
WI_weighted_matrix_16_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_16_in_09_10k[i,j] <- if(distance_ug_16_09[i,j] < 10000) {
      (ug_panel_2009_all[j,5] / (distance_ug_16_09[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_16_in_09_10k) <- sapply(WI_weighted_matrix_16_in_09_10k, is.infinite)

r <- 684
c <- 170
WI_weighted_distance_16_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_16_in_09_10k[i,j] <- if(distance_ug_16_09[i,j] < 10000) {
      (1 / distance_ug_16_09[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_16_in_09_10k) <- sapply(WI_weighted_distance_16_in_09_10k, is.infinite)

ave_WI_weighted_cluster_16_in_09_10k <- matrix(NA,684,1)
ave_WI_weighted_cluster_16_in_09_10k[,1] <- rowSums(WI_weighted_matrix_16_in_09_10k[,1:170], na.rm = TRUE) / rowSums(WI_weighted_distance_16_in_09_10k[,1:170], na.rm = TRUE)

###########
# 2016 in 11
###########
distance_ug_16_11 <- drop_units(distance_ug_16_11)

r <- 684
c <- 400
WI_weighted_matrix_16_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_16_in_11_10k[i,j] <- if(distance_ug_16_11[i,j] < 10000) {
      (ug_panel_2011_all[j,5] / (distance_ug_16_11[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_16_in_11_10k) <- sapply(WI_weighted_matrix_16_in_11_10k, is.infinite)

r <- 684
c <- 400
WI_weighted_distance_16_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_16_in_11_10k[i,j] <- if(distance_ug_16_11[i,j] < 10000) {
      (1 / distance_ug_16_11[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_16_in_11_10k) <- sapply(WI_weighted_distance_16_in_11_10k, is.infinite)

ave_WI_weighted_cluster_16_in_11_10k <- matrix(NA,684,1)
ave_WI_weighted_cluster_16_in_11_10k[,1] <- rowSums(WI_weighted_matrix_16_in_11_10k[,1:400], na.rm = TRUE) / rowSums(WI_weighted_distance_16_in_11_10k[,1:400], na.rm = TRUE)

###########
# 2016 in 2014
##########
distance_ug_16_14 <- drop_units(distance_ug_16_14)

r <- 684
c <- 208
WI_weighted_matrix_16_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_16_in_14_10k[i,j] <- if(distance_ug_16_14[i,j] < 10000) {
      (ug_panel_2014_all[j,5] / (distance_ug_16_14[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_16_in_14_10k) <- sapply(WI_weighted_matrix_16_in_14_10k, is.infinite)

r <- 684
c <- 208
WI_weighted_distance_16_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_16_in_14_10k[i,j] <- if(distance_ug_16_14[i,j] < 10000) {
      (1 / distance_ug_16_14[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_16_in_14_10k) <- sapply(WI_weighted_distance_16_in_14_10k, is.infinite)

ave_WI_weighted_cluster_16_in_14_10k <- matrix(NA,684,1)
ave_WI_weighted_cluster_16_in_14_10k[,1] <- rowSums(WI_weighted_matrix_16_in_14_10k[,1:208], na.rm = TRUE) / rowSums(WI_weighted_distance_16_in_14_10k[,1:208], na.rm = TRUE)

###########
# 2014 in 2006
###########
distance_ug_14_06 <- drop_units(distance_ug_14_06)

r <- 208
c <- 336
WI_weighted_matrix_14_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_14_in_06_10k[i,j] <- if(distance_ug_14_06[i,j] < 10000) {
      (ug_panel_2006_all[j,5] / (distance_ug_14_06[i,j])) } else {
        NA
      }
  }
}

is.na(WI_weighted_matrix_14_in_06_10k) <- sapply(WI_weighted_matrix_14_in_06_10k, is.infinite)

r <- 208
c <- 336
WI_weighted_distance_14_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_14_in_06_10k[i,j] <- if(distance_ug_14_06[i,j] < 10000) {
      (1 / distance_ug_14_06[i,j]) } else {
        NA
      }
  }
}

is.na(WI_weighted_distance_14_in_06_10k) <- sapply(WI_weighted_distance_14_in_06_10k, is.infinite)

ave_WI_weighted_cluster_14_in_06_10k <- matrix(NA,208,1)
ave_WI_weighted_cluster_14_in_06_10k[,1] <- rowSums(WI_weighted_matrix_14_in_06_10k[,1:336], na.rm = TRUE) / rowSums(WI_weighted_distance_14_in_06_10k[,1:336], na.rm = TRUE)

###########
# 2014 in 2009
###########
distance_ug_14_09 <- drop_units(distance_ug_14_09)

r <- 208
c <- 170
WI_weighted_matrix_14_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_14_in_09_10k[i,j] <- if(distance_ug_14_09[i,j] < 10000) {
      (ug_panel_2009_all[j,5] / (distance_ug_14_09[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_14_in_09_10k) <- sapply(WI_weighted_matrix_14_in_09_10k, is.infinite)

r <- 208
c <- 170
WI_weighted_distance_14_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_14_in_09_10k[i,j] <- if(distance_ug_14_09[i,j] < 10000) {
      (1 / distance_ug_14_09[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_14_in_09_10k) <- sapply(WI_weighted_distance_14_in_09_10k, is.infinite)

ave_WI_weighted_cluster_14_in_09_10k <- matrix(NA,208,1)
ave_WI_weighted_cluster_14_in_09_10k[,1] <- rowSums(WI_weighted_matrix_14_in_09_10k[,1:170], na.rm = TRUE) / rowSums(WI_weighted_distance_14_in_09_10k[,1:170], na.rm = TRUE)

############
# 2014 in 2011
############
distance_ug_14_11 <- drop_units(distance_ug_14_11)

r <- 208
c <- 400
WI_weighted_matrix_14_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_14_in_11_10k[i,j] <- if(distance_ug_14_11[i,j] < 10000) {
      (ug_panel_2011_all[j,5] / (distance_ug_14_11[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_14_in_11_10k) <- sapply(WI_weighted_matrix_14_in_11_10k, is.infinite)

r <- 208
c <- 400
WI_weighted_distance_14_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_14_in_11_10k[i,j] <- if(distance_ug_14_11[i,j] < 10000) {
      (1 / distance_ug_14_11[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_14_in_11_10k) <- sapply(WI_weighted_distance_14_in_11_10k, is.infinite)

ave_WI_weighted_cluster_14_in_11_10k <- matrix(NA,208,1)
ave_WI_weighted_cluster_14_in_11_10k[,1] <- rowSums(WI_weighted_matrix_14_in_11_10k[,1:400], na.rm = TRUE) / rowSums(WI_weighted_distance_14_in_11_10k[,1:400], na.rm = TRUE)

############
# 2014 points in 2016
###########
distance_ug_14_16 <- drop_units(distance_ug_16_14)
distance_ug_14_16 <- t(distance_ug_14_16)

r <- 208
c <- 684
WI_weighted_matrix_14_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_14_in_16_10k[i,j] <- if(distance_ug_14_16[i,j] < 10000) {
      (ug_panel_2016_all[j,5] / (distance_ug_14_16[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_14_in_16_10k) <- sapply(WI_weighted_matrix_14_in_16_10k, is.infinite)

r <- 208
c <- 684
WI_weighted_distance_14_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_14_in_16_10k[i,j] <- if(distance_ug_14_16[i,j] < 10000) {
      (1 / distance_ug_14_16[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_14_in_16_10k) <- sapply(WI_weighted_distance_14_in_16_10k, is.infinite)

ave_WI_weighted_cluster_14_in_16_10k <- matrix(NA,208,1)
ave_WI_weighted_cluster_14_in_16_10k[,1] <- rowSums(WI_weighted_matrix_14_in_16_10k[,1:684], na.rm = TRUE) / rowSums(WI_weighted_distance_14_in_16_10k[,1:684], na.rm = TRUE)

############
# 2011 in 2006
###########
distance_ug_11_06 <- drop_units(distance_ug_11_06)

r <- 400
c <- 336
WI_weighted_matrix_11_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_11_in_06_10k[i,j] <- if(distance_ug_11_06[i,j] < 10000) {
      (ug_panel_2006_all[j,5] / (distance_ug_11_06[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_11_in_06_10k) <- sapply(WI_weighted_matrix_11_in_06_10k, is.infinite)

r <- 400
c <- 336
WI_weighted_distance_11_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_11_in_06_10k[i,j] <- if(distance_ug_11_06[i,j] < 10000) {
      (1 / distance_ug_11_06[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_11_in_06_10k) <- sapply(WI_weighted_distance_11_in_06_10k, is.infinite)

ave_WI_weighted_cluster_11_in_06_10k <- matrix(NA,400,1)
ave_WI_weighted_cluster_11_in_06_10k[,1] <- rowSums(WI_weighted_matrix_11_in_06_10k[,1:336], na.rm = TRUE) / rowSums(WI_weighted_distance_11_in_06_10k[,1:336], na.rm = TRUE)

############
# 2011 in 2009
###########
distance_ug_11_09 <- drop_units(distance_ug_11_09)

r <- 400
c <- 170
WI_weighted_matrix_11_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_11_in_09_10k[i,j] <- if(distance_ug_11_09[i,j] < 10000) {
      (ug_panel_2009_all[j,5] / (distance_ug_11_09[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_11_in_09_10k) <- sapply(WI_weighted_matrix_11_in_09_10k, is.infinite)

r <- 400
c <- 170
WI_weighted_distance_11_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_11_in_09_10k[i,j] <- if(distance_ug_11_09[i,j] < 10000) {
      (1 / distance_ug_11_09[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_11_in_09_10k) <- sapply(WI_weighted_distance_11_in_09_10k, is.infinite)

ave_WI_weighted_cluster_11_in_09_10k <- matrix(NA,400,1)
ave_WI_weighted_cluster_11_in_09_10k[,1] <- rowSums(WI_weighted_matrix_11_in_09_10k[,1:170], na.rm = TRUE) / rowSums(WI_weighted_distance_11_in_09_10k[,1:170], na.rm = TRUE)

###########
# 2011 point in 2014
###########
distance_ug_11_14 <- drop_units(distance_ug_14_11)
distance_ug_11_14 <- t(distance_ug_11_14)

r <- 400
c <- 208
WI_weighted_matrix_11_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_11_in_14_10k[i,j] <- if(distance_ug_11_14[i,j] < 10000) {
      (ug_panel_2014_all[j,5] / (distance_ug_11_14[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_11_in_14_10k) <- sapply(WI_weighted_matrix_11_in_14_10k, is.infinite)

r <- 400
c <- 208
WI_weighted_distance_11_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_11_in_14_10k[i,j] <- if(distance_ug_11_14[i,j] < 10000) {
      (1 / distance_ug_11_14[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_11_in_14_10k) <- sapply(WI_weighted_distance_11_in_14_10k, is.infinite)

ave_WI_weighted_cluster_11_in_14_10k <- matrix(NA,400,1)
ave_WI_weighted_cluster_11_in_14_10k[,1] <- rowSums(WI_weighted_matrix_11_in_14_10k[,1:208], na.rm = TRUE) / rowSums(WI_weighted_distance_11_in_14_10k[,1:208], na.rm = TRUE)

###########
# 2011 point in 2016
###########
distance_ug_11_16 <- drop_units(distance_ug_16_11)
distance_ug_11_16 <- t(distance_ug_11_16)

r <- 400
c <- 684
WI_weighted_matrix_11_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_11_in_16_10k[i,j] <- if(distance_ug_11_16[i,j] < 10000) {
      (ug_panel_2016_all[j,5] / (distance_ug_11_16[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_11_in_16_10k) <- sapply(WI_weighted_matrix_11_in_16_10k, is.infinite)

r <- 400
c <- 684
WI_weighted_distance_11_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_11_in_16_10k[i,j] <- if(distance_ug_11_16[i,j] < 10000) {
      (1 / distance_ug_11_16[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_11_in_16_10k) <- sapply(WI_weighted_distance_11_in_16_10k, is.infinite)

ave_WI_weighted_cluster_11_in_16_10k <- matrix(NA,400,1)
ave_WI_weighted_cluster_11_in_16_10k[,1] <- rowSums(WI_weighted_matrix_11_in_16_10k[,1:684], na.rm = TRUE) / rowSums(WI_weighted_distance_11_in_16_10k[,1:684], na.rm = TRUE)

###########
# 2009 in 2006
###########
distance_ug_09_06 <- drop_units(distance_ug_09_06)

r <- 170
c <- 336
WI_weighted_matrix_09_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_09_in_06_10k[i,j] <- if(distance_ug_09_06[i,j] < 10000) {
      (ug_panel_2006_all[j,5] / (distance_ug_09_06[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_09_in_06_10k) <- sapply(WI_weighted_matrix_09_in_06_10k, is.infinite)

r <- 170
c <- 336
WI_weighted_distance_09_in_06_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_09_in_06_10k[i,j] <- if(distance_ug_09_06[i,j] < 10000) {
      (1 / distance_ug_09_06[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_09_in_06_10k) <- sapply(WI_weighted_distance_09_in_06_10k, is.infinite)

ave_WI_weighted_cluster_09_in_06_10k <- matrix(NA,170,1)
ave_WI_weighted_cluster_09_in_06_10k[,1] <- rowSums(WI_weighted_matrix_09_in_06_10k[,1:336], na.rm = TRUE) / rowSums(WI_weighted_distance_09_in_06_10k[,1:336], na.rm = TRUE)

###########
# 2009 in 2011
###########
distance_ug_09_11 <- drop_units(distance_ug_11_09)
distance_ug_09_11 <- t(distance_ug_09_11)

r <- 170
c <- 400
WI_weighted_matrix_09_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_09_in_11_10k[i,j] <- if(distance_ug_09_11[i,j] < 10000) {
      (ug_panel_2011_all[j,5] / (distance_ug_09_11[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_09_in_11_10k) <- sapply(WI_weighted_matrix_09_in_11_10k, is.infinite)

r <- 170
c <- 400
WI_weighted_distance_09_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_09_in_11_10k[i,j] <- if(distance_ug_09_11[i,j] < 10000) {
      (1 / distance_ug_09_11[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_09_in_11_10k) <- sapply(WI_weighted_distance_09_in_11_10k, is.infinite)

ave_WI_weighted_cluster_09_in_11_10k <- matrix(NA,170,1)
ave_WI_weighted_cluster_09_in_11_10k[,1] <- rowSums(WI_weighted_matrix_09_in_11_10k[,1:400], na.rm = TRUE) / rowSums(WI_weighted_distance_09_in_11_10k[,1:400], na.rm = TRUE)

###########
# 2009 in 2014
distance_ug_09_14 <- drop_units(distance_ug_14_09)
distance_ug_09_14 <- t(distance_ug_09_14)

r <- 170
c <- 208
WI_weighted_matrix_09_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_09_in_14_10k[i,j] <- if(distance_ug_09_14[i,j] < 10000) {
      (ug_panel_2014_all[j,5] / (distance_ug_09_14[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_09_in_14_10k) <- sapply(WI_weighted_matrix_09_in_14_10k, is.infinite)

r <- 170
c <- 208
WI_weighted_distance_09_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_09_in_14_10k[i,j] <- if(distance_ug_09_14[i,j] < 10000) {
      (1 / distance_ug_09_14[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_09_in_14_10k) <- sapply(WI_weighted_distance_09_in_14_10k, is.infinite)

ave_WI_weighted_cluster_09_in_14_10k <- matrix(NA,170,1)
ave_WI_weighted_cluster_09_in_14_10k[,1] <- rowSums(WI_weighted_matrix_09_in_14_10k[,1:208], na.rm = TRUE) / rowSums(WI_weighted_distance_09_in_14_10k[,1:208], na.rm = TRUE)

###########
# 2009 in 2016
###########
distance_ug_09_16 <- drop_units(distance_ug_16_09)
distance_ug_09_16 <- t(distance_ug_09_16)

r <- 170
c <- 684
WI_weighted_matrix_09_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_09_in_16_10k[i,j] <- if(distance_ug_09_16[i,j] < 10000) {
      (ug_panel_2016_all[j,5] / (distance_ug_09_16[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_09_in_16_10k) <- sapply(WI_weighted_matrix_09_in_16_10k, is.infinite)

r <- 170
c <- 684
WI_weighted_distance_09_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_09_in_16_10k[i,j] <- if(distance_ug_09_16[i,j] < 10000) {
      (1 / distance_ug_09_16[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_09_in_16_10k) <- sapply(WI_weighted_distance_09_in_16_10k, is.infinite)

ave_WI_weighted_cluster_09_in_16_10k <- matrix(NA,170,1)
ave_WI_weighted_cluster_09_in_16_10k[,1] <- rowSums(WI_weighted_matrix_09_in_16_10k[,1:684], na.rm = TRUE) / rowSums(WI_weighted_distance_09_in_16_10k[,1:684], na.rm = TRUE)

###########
# 2006 in 2009
###########
distance_ug_06_09 <- drop_units(distance_ug_09_06)
distance_ug_06_09 <- t(distance_ug_06_09)

r <- 336
c <- 170
WI_weighted_matrix_06_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_06_in_09_10k[i,j] <- if(distance_ug_06_09[i,j] < 10000) {
      (ug_panel_2009_all[j,5] / (distance_ug_06_09[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_06_in_09_10k) <- sapply(WI_weighted_matrix_06_in_09_10k, is.infinite)

r <- 336
c <- 170
WI_weighted_distance_06_in_09_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_06_in_09_10k[i,j] <- if(distance_ug_06_09[i,j] < 10000) {
      (1 / distance_ug_06_09[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_06_in_09_10k) <- sapply(WI_weighted_distance_06_in_09_10k, is.infinite)

ave_WI_weighted_cluster_06_in_09_10k <- matrix(NA,336,1)
ave_WI_weighted_cluster_06_in_09_10k[,1] <- rowSums(WI_weighted_matrix_06_in_09_10k[,1:170], na.rm = TRUE) / rowSums(WI_weighted_distance_06_in_09_10k[,1:170], na.rm = TRUE)

###########
# 2006 point in 2011
###########
distance_ug_06_11 <- drop_units(distance_ug_11_06)
distance_ug_06_11 <- t(distance_ug_06_11)

r <- 336
c <- 400
WI_weighted_matrix_06_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_06_in_11_10k[i,j] <- if(distance_ug_06_11[i,j] < 10000) {
      (ug_panel_2011_all[j,5] / (distance_ug_06_11[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_06_in_11_10k) <- sapply(WI_weighted_matrix_06_in_11_10k, is.infinite)

r <- 336
c <- 400
WI_weighted_distance_06_in_11_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_06_in_11_10k[i,j] <- if(distance_ug_06_11[i,j] < 10000) {
      (1 / distance_ug_06_11[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_06_in_11_10k) <- sapply(WI_weighted_distance_06_in_11_10k, is.infinite)

ave_WI_weighted_cluster_06_in_11_10k <- matrix(NA,336,1)
ave_WI_weighted_cluster_06_in_11_10k[,1] <- rowSums(WI_weighted_matrix_06_in_11_10k[,1:400], na.rm = TRUE) / rowSums(WI_weighted_distance_06_in_11_10k[,1:400], na.rm = TRUE)

###########
# 2006 in 2014
###########
distance_ug_06_14 <- drop_units(distance_ug_14_06)
distance_ug_06_14 <- t(distance_ug_06_14)

r <- 336
c <- 208
WI_weighted_matrix_06_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_06_in_14_10k[i,j] <- if(distance_ug_06_14[i,j] < 10000) {
      (ug_panel_2014_all[j,5] / (distance_ug_06_14[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_06_in_14_10k) <- sapply(WI_weighted_matrix_06_in_14_10k, is.infinite)

r <- 336
c <- 208
WI_weighted_distance_06_in_14_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_06_in_14_10k[i,j] <- if(distance_ug_06_14[i,j] < 10000) {
      (1 / distance_ug_06_14[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_06_in_14_10k) <- sapply(WI_weighted_distance_06_in_14_10k, is.infinite)

ave_WI_weighted_cluster_06_in_14_10k <- matrix(NA,336,1)
ave_WI_weighted_cluster_06_in_14_10k[,1] <- rowSums(WI_weighted_matrix_06_in_14_10k[,1:208], na.rm = TRUE) / rowSums(WI_weighted_distance_06_in_14_10k[,1:208], na.rm = TRUE)

###########
# 2006 in 2016
###########
distance_ug_06_16 <- drop_units(distance_ug_16_06)
distance_ug_06_16 <- t(distance_ug_06_16)

r <- 336
c <- 684
WI_weighted_matrix_06_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_matrix_06_in_16_10k[i,j] <- if(distance_ug_06_16[i,j] < 10000) {
      (ug_panel_2016_all[j,5] / (distance_ug_06_16[i,j])) } else {
        NA
      }
  }
}
is.na(WI_weighted_matrix_06_in_16_10k) <- sapply(WI_weighted_matrix_06_in_16_10k, is.infinite)

r <- 336
c <- 684
WI_weighted_distance_06_in_16_10k  <- matrix(NA,r,c)
for(i in 1:r){
  for(j in 1:c ) {
    WI_weighted_distance_06_in_16_10k[i,j] <- if(distance_ug_06_16[i,j] < 10000) {
      (1 / distance_ug_06_16[i,j]) } else {
        NA
      }
  }
}
is.na(WI_weighted_distance_06_in_16_10k) <- sapply(WI_weighted_distance_06_in_16_10k, is.infinite)

ave_WI_weighted_cluster_06_in_16_10k <- matrix(NA,336,1)
ave_WI_weighted_cluster_06_in_16_10k[,1] <- rowSums(WI_weighted_matrix_06_in_16_10k[,1:684], na.rm = TRUE) / rowSums(WI_weighted_distance_06_in_16_10k[,1:684], na.rm = TRUE)


ug_2016_all_preds <- cbind(ug_panel_2016_all, ave_WI_weighted_cluster_16_in_06_10k,
                           ave_WI_weighted_cluster_16_in_09_10k, ave_WI_weighted_cluster_16_in_11_10k, 
                           ave_WI_weighted_cluster_16_in_14_10k)

ug_2014_all_preds <- cbind(ug_panel_2014_all, ave_WI_weighted_cluster_14_in_06_10k,
                           ave_WI_weighted_cluster_14_in_09_10k, ave_WI_weighted_cluster_14_in_11_10k, 
                           ave_WI_weighted_cluster_14_in_16_10k)

ug_2011_all_preds <- cbind(ug_panel_2011_all, ave_WI_weighted_cluster_11_in_06_10k,
                           ave_WI_weighted_cluster_11_in_09_10k, ave_WI_weighted_cluster_11_in_14_10k, 
                           ave_WI_weighted_cluster_11_in_16_10k)

ug_2009_all_preds <- cbind(ug_panel_2009_all, ave_WI_weighted_cluster_09_in_06_10k,
                           ave_WI_weighted_cluster_09_in_11_10k, ave_WI_weighted_cluster_09_in_14_10k, 
                           ave_WI_weighted_cluster_09_in_16_10k)

ug_2006_all_preds <- cbind(ug_panel_2006_all, ave_WI_weighted_cluster_06_in_09_10k,
                           ave_WI_weighted_cluster_06_in_11_10k, ave_WI_weighted_cluster_06_in_14_10k, 
                           ave_WI_weighted_cluster_06_in_16_10k)

###########
# Now I am merging with a 2km treat and control #
###########

TC_2u2y <- read.csv("TC_2u2y.csv")
hold_tc <- TC_2u2y
hold_tc <- hold_tc[,c(2,8)]

# for 2016 
all_2016_treat_and_control <- merge(hold_tc, ug_2016_all_preds, by = "DHSID")
all_2016_treat_and_control_noNA <- na.omit(all_2016_treat_and_control)

# for 2014 
all_2014_treat_and_control <- merge(hold_tc, ug_2014_all_preds, by = "DHSID")
all_2014_treat_and_control_noNA <- na.omit(all_2014_treat_and_control)

# for 2011 
all_2011_treat_and_control <- merge(hold_tc, ug_2011_all_preds, by = "DHSID")
all_2011_treat_and_control_noNA <- na.omit(all_2011_treat_and_control)

# for 2009 
all_2009_treat_and_control <- merge(hold_tc, ug_2009_all_preds, by = "DHSID")
all_2009_treat_and_control_noNA <- na.omit(all_2009_treat_and_control)

# for 2006 
all_2006_treat_and_control <- merge(hold_tc, ug_2006_all_preds, by = "DHSID")
all_2006_treat_and_control_noNA <- na.omit(all_2006_treat_and_control)

###########
##### Here I am merging the df's and doing DiD ####
###########

sub_2016_2016 <- all_2016_treat_and_control_noNA[,c(1:5,6)]
colnames(sub_2016_2016)[6] <- "WI"
sub_2016_2016[,3] <- 2016

sub_2016_2014 <- all_2016_treat_and_control_noNA[,c(1:5,12)]
colnames(sub_2016_2014)[6] <- "WI"
sub_2016_2014[,3] <- 2014

sub_2016_2011 <- all_2016_treat_and_control_noNA[,c(1:5,11)]
colnames(sub_2016_2011)[6] <- "WI"
sub_2016_2011[,3] <- 2011

sub_2016_2009 <- all_2016_treat_and_control_noNA[,c(1:5,10)]
colnames(sub_2016_2009)[6] <- "WI"
sub_2016_2009[,3] <- 2009

sub_2016_2006 <- all_2016_treat_and_control_noNA[,c(1:5,9)]
colnames(sub_2016_2006)[6] <- "WI"
sub_2016_2006[,3] <- 2006

#splitting 2014
sub_2014_2016 <- all_2014_treat_and_control_noNA[,c(1:5,6)]
colnames(sub_2014_2016)[6] <- "WI"
sub_2014_2016[,3] <- 2016

sub_2014_2014 <- all_2014_treat_and_control_noNA[,c(1:5,12)]
colnames(sub_2014_2014)[6] <- "WI"
sub_2014_2014[,3] <- 2014

sub_2014_2011 <- all_2014_treat_and_control_noNA[,c(1:5,11)]
colnames(sub_2014_2011)[6] <- "WI"
sub_2014_2011[,3] <- 2011

sub_2014_2009 <- all_2014_treat_and_control_noNA[,c(1:5,10)]
colnames(sub_2014_2009)[6] <- "WI"
sub_2014_2009[,3] <- 2009

sub_2014_2006 <- all_2014_treat_and_control_noNA[,c(1:5,9)]
colnames(sub_2014_2006)[6] <- "WI"
sub_2014_2006[,3] <- 2006

#splitting 2011
sub_2011_2016 <- all_2011_treat_and_control_noNA[,c(1:5,6)]
colnames(sub_2011_2016)[6] <- "WI"
sub_2011_2016[,3] <- 2016

sub_2011_2014 <- all_2011_treat_and_control_noNA[,c(1:5,12)]
colnames(sub_2011_2014)[6] <- "WI"
sub_2011_2014[,3] <- 2014

sub_2011_2011 <- all_2011_treat_and_control_noNA[,c(1:5,11)]
colnames(sub_2011_2011)[6] <- "WI"
sub_2011_2011[,3] <- 2011

sub_2011_2009 <- all_2011_treat_and_control_noNA[,c(1:5,10)]
colnames(sub_2011_2009)[6] <- "WI"
sub_2011_2009[,3] <- 2009

sub_2011_2006 <- all_2011_treat_and_control_noNA[,c(1:5,9)]
colnames(sub_2011_2006)[6] <- "WI"
sub_2011_2006[,3] <- 2006

#splitting 2009
sub_2009_2016 <- all_2009_treat_and_control_noNA[,c(1:5,6)]
colnames(sub_2009_2016)[6] <- "WI"
sub_2009_2016[,3] <- 2016

sub_2009_2014 <- all_2009_treat_and_control_noNA[,c(1:5,12)]
colnames(sub_2009_2014)[6] <- "WI"
sub_2009_2014[,3] <- 2014

sub_2009_2011 <- all_2009_treat_and_control_noNA[,c(1:5,11)]
colnames(sub_2009_2011)[6] <- "WI"
sub_2009_2011[,3] <- 2011

sub_2009_2009 <- all_2009_treat_and_control_noNA[,c(1:5,10)]
colnames(sub_2009_2009)[6] <- "WI"
sub_2009_2009[,3] <- 2009

sub_2009_2006 <- all_2009_treat_and_control_noNA[,c(1:5,9)]
colnames(sub_2009_2006)[6] <- "WI"
sub_2009_2006[,3] <- 2006

#splitting 2006
sub_2006_2016 <- all_2006_treat_and_control_noNA[,c(1:5,6)]
colnames(sub_2006_2016)[6] <- "WI"
sub_2006_2016[,3] <- 2016

sub_2006_2014 <- all_2006_treat_and_control_noNA[,c(1:5,12)]
colnames(sub_2006_2014)[6] <- "WI"
sub_2006_2014[,3] <- 2014

sub_2006_2011 <- all_2006_treat_and_control_noNA[,c(1:5,11)]
colnames(sub_2006_2011)[6] <- "WI"
sub_2006_2011[,3] <- 2011

sub_2006_2009 <- all_2006_treat_and_control_noNA[,c(1:5,10)]
colnames(sub_2006_2009)[6] <- "WI"
sub_2006_2009[,3] <- 2009

sub_2006_2006 <- all_2006_treat_and_control_noNA[,c(1:5,9)]
colnames(sub_2006_2006)[6] <- "WI"
sub_2006_2006[,3] <- 2006

long_sub_to_aggregated_noNA_b <- rbind(sub_2016_2016, sub_2016_2014, sub_2016_2011, sub_2016_2009, sub_2016_2006,
                                       sub_2014_2016, sub_2014_2014, sub_2014_2011, sub_2014_2009, sub_2014_2006,
                                       sub_2011_2016, sub_2011_2014, sub_2011_2011, sub_2011_2009, sub_2011_2006,
                                       sub_2009_2016, sub_2009_2014, sub_2009_2011, sub_2009_2009, sub_2009_2006,
                                       sub_2006_2016, sub_2006_2014, sub_2006_2011, sub_2006_2009, sub_2006_2006)

treated_fun <- function(x) { 
  if(x == 2006 | x == 2009 ) y <- 0
  if(x == 2011 |x == 2014 | x == 2016) y <- 1
  return(y)
}

long_sub_to_aggregated_noNA_b$treat_year <- sapply(long_sub_to_aggregated_noNA_b$year, treated_fun)
long_sub_to_aggregated_noNA_b$treat_ind <- long_sub_to_aggregated_noNA_b$treated_unit * long_sub_to_aggregated_noNA_b$treat_year

idw_results_df <- as.data.frame(matrix(NA,1,2))

idw_dd <- summary(feols(WI ~ treat_ind | DHSID + year, data=long_sub_to_aggregated_noNA_b))
idw_results_df[1,1]  <- as.data.frame(coeftable(idw_dd)[1])
idw_results_df[1,2]  <- as.data.frame(coeftable(idw_dd)[2])

colnames(idw_results_df) <- c("idw_est", "idw_se")
write.csv(idw_results_df, "idw_results_df.csv", row.names = FALSE)









