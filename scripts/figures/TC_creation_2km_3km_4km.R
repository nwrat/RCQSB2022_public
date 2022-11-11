# This script creates the treatment and control groups for 2km, 3km and 4km buffer. 

rm(list=ls())

library(dplyr)
library(sf)
library(reshape2)

setwd("data/figure_and_input_data/")

ug_dist_2010 <- st_read(dsn = "Ug_Aug_2010_shapefile.shp", layer = "Ug_Aug_2010_shapefile", stringsAsFactors = F)
ug_dist_2013  <- st_read(dsn = "UG_2013_Shapefile_Aug30.shp", layer = "UG_2013_Shapefile_Aug30", stringsAsFactors = F)
ug_dist_2016  <- st_read(dsn = "Distribution_Lines_2016_Operational.shp", layer = "Distribution_Lines_2016_Operational", stringsAsFactors = F)

dhs_survey_locations <- read.csv("WI_pca_output.csv", check.names=FALSE)

dhs_survey_locations_ug <- dhs_survey_locations
dhs_survey_locations_ug$country  <- substr(dhs_survey_locations_ug$DHSID, 1, 2)

dhs_survey_locations_ug <- subset(dhs_survey_locations_ug, country == "UG") #1,798
dhs_survey_locations_ug_ll <- dhs_survey_locations_ug[,c(4,3)]

sf::sf_use_s2(FALSE)
ug_dhs_ll_sf <- st_as_sf(dhs_survey_locations_ug_ll, coords = c("long","lat"), crs = ("+init=epsg:4326"))

ug_grid_2010_shape_2 <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 2000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2013_shape_2 <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 2000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2016_shape_2 <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 2000) %>%
  st_transform("+init=epsg:4326")

village_intersects_2013_line_2 <- st_intersects(ug_dhs_ll_sf, ug_grid_2013_shape_2, sparse = FALSE) 
true_index_village_13_2 <- which(apply(village_intersects_2013_line_2,1,any))
true_index_village_13_2_df <- as.data.frame(true_index_village_13_2)
length(true_index_village_13_2_df[,1])
true_matrix_village_13_2 <- matrix(NA,735,2)
true_matrix_village_13_2[,1] <- true_index_village_13_2
#then 2010, then next line

# 2010 intersection
village_intersects_2010_line_2 <- st_intersects(ug_dhs_ll_sf, ug_grid_2010_shape_2, sparse = FALSE) 
true_index_village_10_2 <- which(apply(village_intersects_2010_line_2,1,any))
true_index_village_10_2_df <- as.data.frame(true_index_village_10_2)
length(true_index_village_10_2_df[,1]) #648
true_matrix_village_10_2 <- matrix(NA,735,2) #same as 2013
true_matrix_village_10_2[1:648,1] <- true_index_village_10_2
true_matrix_village_10_2[649:735,2] <- 0  

true_matrix_village_13_2[,2] <- true_matrix_village_10_2[,1]
true_matrix_village_13_2 <- as.data.frame(true_matrix_village_13_2)

# connected in 2016
village_intersects_2016_line_2 <- st_intersects(ug_dhs_ll_sf, ug_grid_2016_shape_2, sparse = FALSE) #this is the line that gives you treated
true_index_village_16_2 <- which(apply(village_intersects_2016_line_2,1,any))
true_index_village_16_2_df <- as.data.frame(true_index_village_16_2)
length(true_index_village_16_2_df[,1]) #899
true_matrix_village_16_2 <- matrix(NA,899,2)
true_matrix_village_16_2[1:899,1] <- true_index_village_16_2

# non-electrified communities in 2016 
control_for_village_WI_2 <- matrix(NA,1798,4)
colnames(control_for_village_WI_2) <- c("row_num", "name", "lat", "long")
control_for_village_WI_2[,1] <- rep(1:1798) 
control_for_village_WI_2[,2] <- as.character(dhs_survey_locations_ug[,1])
control_for_village_WI_2[,3] <- dhs_survey_locations_ug[,3]
control_for_village_WI_2[,4] <- dhs_survey_locations_ug[,4]

# connected communities in 2016, treated
connected_village_WI_2 <- control_for_village_WI_2[control_for_village_WI_2[,1] %in% true_matrix_village_16_2[,1], ]
connected_village_WI_2 <- as.data.frame(connected_village_WI_2)


control_for_village_WI_2_list <- control_for_village_WI_2[!control_for_village_WI_2[,1] %in% true_matrix_village_16_2[,1], ]
control_for_village_WI_2_list <- as.data.frame(control_for_village_WI_2_list)

treatment_villages_2 <- matrix(NA,735+648,1)
treatment_villages_2[1:735,] <- true_matrix_village_13_2[1:735,1]
treatment_villages_2[736:1383,] <- true_matrix_village_13_2[1:648,2]
treatment_villages_2_list <- as.data.frame(as.numeric(names(which(table(treatment_villages_2)==1))))
colnames(treatment_villages_2_list) <- c("treated_units")

treatment_for_village_WI_2 <- matrix(NA,1798,4)
colnames(treatment_for_village_WI_2) <- c("row_num", "name", "lat", "long")
treatment_for_village_WI_2[,1] <- rep(1:1798) 
treatment_for_village_WI_2[,2] <- as.character(dhs_survey_locations_ug[,1])
treatment_for_village_WI_2[,3] <- dhs_survey_locations_ug[,3]
treatment_for_village_WI_2[,4] <- dhs_survey_locations_ug[,4]

treatment_for_village_WI_2 <- treatment_for_village_WI_2[treatment_for_village_WI_2[,1] %in% treatment_villages_2_list[,1], ]
treatment_for_village_WI_2 <- as.data.frame(treatment_for_village_WI_2)

tc_2km_test <- rbind(control_for_village_WI_2_list, treatment_for_village_WI_2)
tc_2km_test <- tc_2km_test[!(duplicated(tc_2km_test) | duplicated(tc_2km_test, fromLast = TRUE)), ]

tc_2km_test$treated <- 0
tc_2km_test[889:964,5] <- 1
colnames(tc_2km_test)[2] <- "DHSID"

write.csv(tc_2km_test, "tc_2km_ug.csv", row.names = FALSE)

###########
# 3km treatment and control creation
###########
ug_grid_2010_shape_3km <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2013_shape_3km <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2016_shape_3km <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

village_intersects_2013_line_3 <- st_intersects(ug_dhs_ll_sf, ug_grid_2013_shape_3km, sparse = FALSE) 
true_index_village_13_3 <- which(apply(village_intersects_2013_line_3,1,any))
true_index_village_13_3_df <- as.data.frame(true_index_village_13_3)
length(true_index_village_13_3_df[,1])
true_matrix_village_13_3 <- matrix(NA,858,2)
true_matrix_village_13_3[,1] <- true_index_village_13_3

village_intersects_2010_line_3 <- st_intersects(ug_dhs_ll_sf, ug_grid_2010_shape_3km, sparse = FALSE) 
true_index_village_10_3 <- which(apply(village_intersects_2010_line_3,1,any))
true_index_village_10_3_df <- as.data.frame(true_index_village_10_3)
length(true_index_village_10_3_df[,1]) #758
true_matrix_village_10_3 <- matrix(NA,858,3) #same as 2013
true_matrix_village_10_3[1:758,1] <- true_index_village_10_3
true_matrix_village_10_3[759:853,2] <- 0  

true_matrix_village_13_3[,2] <- true_matrix_village_10_3[,1]
true_matrix_village_13_3 <- as.data.frame(true_matrix_village_13_3)

village_intersects_2016_line_3 <- st_intersects(ug_dhs_ll_sf, ug_grid_2016_shape_3km, sparse = FALSE) #this is the line that gives you treated
true_index_village_16_3 <- which(apply(village_intersects_2016_line_3,1,any))
true_index_village_16_3_df <- as.data.frame(true_index_village_16_3)
length(true_index_village_16_3_df[,1]) #
true_matrix_village_16_3 <- matrix(NA,1053,2)
true_matrix_village_16_3[1:1053,1] <- true_index_village_16_3

control_for_village_WI_3 <- matrix(NA,1798,4)
colnames(control_for_village_WI_3) <- c("row_num", "name", "lat", "long")
control_for_village_WI_3[,1] <- rep(1:1798) 
control_for_village_WI_3[,2] <- as.character(dhs_survey_locations_ug[,1])
control_for_village_WI_3[,3] <- dhs_survey_locations_ug[,3]
control_for_village_WI_3[,4] <- dhs_survey_locations_ug[,4]

connected_village_WI_3 <- control_for_village_WI_3[control_for_village_WI_3[,1] %in% true_matrix_village_16_3[,1], ]
connected_village_WI_3 <- as.data.frame(connected_village_WI_3)

control_for_village_WI_3_list <- control_for_village_WI_3[!control_for_village_WI_3[,1] %in% true_matrix_village_16_3[,1], ]
control_for_village_WI_3_list <- as.data.frame(control_for_village_WI_3_list)

treatment_villages_3 <- matrix(NA,858+758,1)
treatment_villages_3[1:858,] <- true_matrix_village_13_3[1:858,1]
treatment_villages_3[859:1616,] <- true_matrix_village_13_3[1:758,2]
treatment_villages_3_list <- as.data.frame(as.numeric(names(which(table(treatment_villages_3)==1))))
colnames(treatment_villages_3_list) <- c("treated_units")

treatment_for_village_WI_3 <- matrix(NA,1798,4)
colnames(treatment_for_village_WI_3) <- c("row_num", "name", "lat", "long")
treatment_for_village_WI_3[,1] <- rep(1:1798) 
treatment_for_village_WI_3[,2] <- as.character(dhs_survey_locations_ug[,1])
treatment_for_village_WI_3[,3] <- dhs_survey_locations_ug[,3]
treatment_for_village_WI_3[,4] <- dhs_survey_locations_ug[,4]

treatment_for_village_WI_3 <- treatment_for_village_WI_3[treatment_for_village_WI_3[,1] %in% treatment_villages_3_list[,1], ]
treatment_for_village_WI_3 <- as.data.frame(treatment_for_village_WI_3)

tc_3km <- rbind(control_for_village_WI_3_list, treatment_for_village_WI_3)
tc_3km <- tc_3km[!(duplicated(tc_3km) | duplicated(tc_3km, fromLast = TRUE)), ]

tc_3km$treated <- 0
tc_3km[736:825,5] <- 1
colnames(tc_3km)[2] <- "DHSID"

write.csv(tc_3km, "tc_3km_ug.csv", row.names = FALSE)

###########
# 4km treatment and control creation
###########
ug_grid_2010_shape_4km <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>% #this is the version it was made in...?
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326") #this is the one i want at the end

ug_grid_2013_shape_4km <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2016_shape_4km <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326")

village_intersects_2013_line_4 <- st_intersects(ug_dhs_ll_sf, ug_grid_2013_shape_4km, sparse = FALSE) 
true_index_village_13_4 <- which(apply(village_intersects_2013_line_4,1,any))
true_index_village_13_4_df <- as.data.frame(true_index_village_13_4)
length(true_index_village_13_4_df[,1])
true_matrix_village_13_4 <- matrix(NA,955,2)
true_matrix_village_13_4[,1] <- true_index_village_13_4

village_intersects_2010_line_4 <- st_intersects(ug_dhs_ll_sf, ug_grid_2010_shape_4km, sparse = FALSE) 
true_index_village_10_4 <- which(apply(village_intersects_2010_line_4,1,any))
true_index_village_10_4_df <- as.data.frame(true_index_village_10_4)
length(true_index_village_10_4_df[,1]) #
true_matrix_village_10_4 <- matrix(NA,955,3) #same as 2013
true_matrix_village_10_4[1:850,1] <- true_index_village_10_4
true_matrix_village_10_4[851:955,2] <- 0  

true_matrix_village_13_4[,2] <- true_matrix_village_10_4[,1]
true_matrix_village_13_4 <- as.data.frame(true_matrix_village_13_4)

village_intersects_2016_line_4 <- st_intersects(ug_dhs_ll_sf, ug_grid_2016_shape_4km, sparse = FALSE) #this is the line that gives you treated
true_index_village_16_4 <- which(apply(village_intersects_2016_line_4,1,any))
true_index_village_16_4_df <- as.data.frame(true_index_village_16_4)
length(true_index_village_16_4_df[,1]) #
true_matrix_village_16_4 <- matrix(NA,1173,2)
true_matrix_village_16_4[1:1173,1] <- true_index_village_16_4

control_for_village_WI_4 <- matrix(NA,1798,4)
colnames(control_for_village_WI_4) <- c("row_num", "name", "lat", "long")
control_for_village_WI_4[,1] <- rep(1:1798) 
control_for_village_WI_4[,2] <- as.character(dhs_survey_locations_ug[,1])
control_for_village_WI_4[,3] <- dhs_survey_locations_ug[,3]
control_for_village_WI_4[,4] <- dhs_survey_locations_ug[,4]

connected_village_WI_4 <- control_for_village_WI_4[control_for_village_WI_4[,1] %in% true_matrix_village_16_4[,1], ]
connected_village_WI_4 <- as.data.frame(connected_village_WI_4)

control_for_village_WI_4_list <- control_for_village_WI_4[!control_for_village_WI_4[,1] %in% true_matrix_village_16_4[,1], ]
control_for_village_WI_4_list <- as.data.frame(control_for_village_WI_4_list)

treatment_villages_4 <- matrix(NA,955+850,1)
treatment_villages_4[1:955,] <- true_matrix_village_13_4[1:955,1]
treatment_villages_4[956:1805,] <- true_matrix_village_13_4[1:850,2]
treatment_villages_4_list <- as.data.frame(as.numeric(names(which(table(treatment_villages_4)==1))))
colnames(treatment_villages_4_list) <- c("treated_units")

treatment_for_village_WI_4 <- matrix(NA,1798,4)
colnames(treatment_for_village_WI_4) <- c("row_num", "name", "lat", "long")
treatment_for_village_WI_4[,1] <- rep(1:1798) 
treatment_for_village_WI_4[,2] <- as.character(dhs_survey_locations_ug[,1])
treatment_for_village_WI_4[,3] <- dhs_survey_locations_ug[,3]
treatment_for_village_WI_4[,4] <- dhs_survey_locations_ug[,4]

treatment_for_village_WI_4 <- treatment_for_village_WI_4[treatment_for_village_WI_4[,1] %in% treatment_villages_4_list[,1], ]
treatment_for_village_WI_4 <- as.data.frame(treatment_for_village_WI_4)

tc_4km <- rbind(control_for_village_WI_4_list, treatment_for_village_WI_4)
tc_4km <- tc_4km[!(duplicated(tc_4km) | duplicated(tc_4km, fromLast = TRUE)), ]

tc_4km$treated <- 0
tc_4km[614:706,5] <- 1
colnames(tc_4km)[2] <- "DHSID"

write.csv(tc_4km, "tc_4km_ug.csv", row.names = FALSE)













