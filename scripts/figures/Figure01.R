# This script creates Figure 1, a - c. 

library(rgdal)
library(ggplot2)
library(sf)

setwd("data/figure_and_input_data/")

###########
#Figure 1a
###########
africa_elec_polygon <- readOGR(dsn = "africa_polygon_w_elec_access.shp", layer = "africa_polygon_w_elec_access", stringsAsFactors = F)
africa_elec_polygon_sf <- st_as_sf(africa_elec_polygon, coords = c("long","lat"), crs = ("+init=epsg:4326"))

pdf("../../figures/raw/WI_Fig_1a.pdf")   
WI_Fig_1a <-ggplot() +
  geom_sf(data = africa_elec_polygon_sf, aes(fill = factor(elec_index))) + 
  scale_fill_manual(values = c("-1" = "grey80", "0"= "white",  
                               "1" = "lightblue1", "2" = "lightblue3", "3" = "skyblue3", 
                               "4" = "steelblue1", "5" = "royalblue", "6" = "navy"), name = "", 
                    labels = c("NA", "0", "< 1m", "1m - 5m", "5m - 10m", "10m - 20m", "20m - 50m", "> 50m")) +
  geom_path(data=UGA_border, aes(y=lat, x=long, group = group), color = "orangered", size = 1) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

WI_Fig_1a
dev.off()

###########
#Figure 1b
###########
UGA_border <- readOGR(dsn = "gadm36_UGA_0.shp", layer = "gadm36_UGA_0", stringsAsFactors = F)
KEN_border <- readOGR(dsn = "gadm36_KEN_0.shp", layer = "gadm36_KEN_0", stringsAsFactors = F)
TZA_border <- readOGR(dsn = "gadm36_TZA_0.shp", layer = "gadm36_TZA_0", stringsAsFactors = F)
RWA_border <- readOGR(dsn = "gadm36_RWA_0.shp", layer = "gadm36_RWA_0", stringsAsFactors = F)
SSD_border <- readOGR(dsn = "gadm36_SSD_0.shp", layer = "gadm36_SSD_0", stringsAsFactors = F)
COD_border <- readOGR(dsn = "gadm36_COD_0.shp", layer = "gadm36_COD_0", stringsAsFactors = F)
ug_water_bodies <- readOGR(dsn = "Ug_Waterbodies.shp", layer = "Ug_Waterbodies", stringsAsFactors = F)
pca_output <- read.csv("wi_pca_output.csv")
pca_output$country = substr(pca_output$DHSID, 1, 2)
pca_output_ug <- subset(pca_output, country == "UG")
pca_output_ug_2016 <- subset(pca_output_ug, year == 2016)
  
ug_dist_2010 <- st_read(dsn = "Ug_Aug_2010_shapefile.shp", layer = "Ug_Aug_2010_shapefile", stringsAsFactors = F)
ug_grid_2010_shape <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 1000) %>%
  st_transform("+init=epsg:4326")

ug_dist_2013 <- st_read(dsn = "UG_2013_Shapefile_Aug30.shp", layer = "UG_2013_Shapefile_Aug30", stringsAsFactors = F)
ug_grid_2013_shape <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 1000) %>%
  st_transform("+init=epsg:4326")

ug_dist_2016 <- st_read(dsn = "Distribution_Lines_2016_Operational.shp", layer = "Distribution_Lines_2016_Operational", stringsAsFactors = F)
ug_grid_2016_shape <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 1000) %>%
  st_transform("+init=epsg:4326")

pdf("../../figures/raw/WI_Fig_1b.pdf") 
WI_Fig_1b <- ggplot() +
  geom_polygon(data=UGA_border, aes(x=long, y=lat, group = group), fill = "grey70") +
  geom_polygon(data=ug_water_bodies, aes(x=long, y=lat, group = group), fill = "lightcyan") +
  geom_polygon(data=KEN_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=TZA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=RWA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=SSD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=COD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_path(data=UGA_border, aes(y=lat, x=long, group = group), size = .25) +
  geom_sf(data = ug_grid_2016_shape, colour = "purple1", fill = "purple1", size = .05) + #before i just had color set and not fill
  geom_sf(data = ug_grid_2013_shape, fill = "red", color = "red", size = .05) +
  geom_sf(data = ug_grid_2010_shape, colour = "dodgerblue1", fill = "dodgerblue1",  size = .1) +
  geom_point(data = pca_output_ug_2016, aes(x = long, y = lat), size = .25) +
  coord_sf(xlim = c(29.55, 34.85), ylim = c(-1.3, 4.05)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),       
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),       
        axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
WI_Fig_1b
dev.off()

###########
#Figure 1c
###########
dhs_survey_locations <- read.csv("WI_pca_output.csv", check.names=FALSE)
dhs_country_polygon <- readOGR(dsn = "dhs_survey_polygon_file.shp", layer = "dhs_survey_polygon_file", stringsAsFactors = F)
dhs_country_polygon_sf <- st_as_sf(dhs_country_polygon, coords = c("long","lat"), crs = ("+init=epsg:4326"))

pdf("../../figures/raw/WI_Fig_1c.pdf") 
WI_Fig_1c <- ggplot() +
  geom_sf(data = africa_elec_polygon_sf, aes(), fill = "grey70", color = "white") +
  geom_sf(data = dhs_country_polygon_sf, aes(), fill = "grey70", color = "black") +
  geom_point(data=dhs_survey_locations, 
             aes(x=long, y = lat, color = -base), size = .25, show.legend = FALSE) + 
  scale_color_gradientn(colours = rainbow(5)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

WI_Fig_1c
dev.off()














