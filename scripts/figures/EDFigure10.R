# This code creates Extended Data Figure 10. 

rm(list=ls())

library(ggplot2)
library(raster)
library(rgdal)
library(sf)
 
setwd("data/figure_and_input_data/")

UGA_border <- readOGR(dsn = "gadm36_UGA_0.shp", layer = "gadm36_UGA_0", stringsAsFactors = F)
ug_water_bodies <- readOGR(dsn = "Ug_Waterbodies.shp", layer = "Ug_Waterbodies", stringsAsFactors = F)
COD_border <- readOGR(dsn = "gadm36_COD_0.shp", layer = "gadm36_COD_0", stringsAsFactors = F)
KEN_border <- readOGR(dsn = "gadm36_KEN_0.shp", layer = "gadm36_KEN_0", stringsAsFactors = F)
TZA_border <- readOGR(dsn = "gadm36_TZA_0.shp", layer = "gadm36_TZA_0", stringsAsFactors = F)
RWA_border <- readOGR(dsn = "gadm36_RWA_0.shp", layer = "gadm36_RWA_0", stringsAsFactors = F)
SSD_border <- readOGR(dsn = "gadm36_SSD_0.shp", layer = "gadm36_SSD_0", stringsAsFactors = F)

roads_UG_2012 <- st_read("roads_UG_2012.shp",  stringsAsFactors = FALSE)
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

# Main Figure
pdf("../../figures/raw/EDFigure10.pdf") 
EDFigure10 <- ggplot() +
  geom_polygon(data=UGA_border, aes(x=long, y=lat, group = group), fill = "white") +
  geom_polygon(data=ug_water_bodies, aes(x=long, y=lat, group = group), fill = "lightcyan") +
  geom_polygon(data=KEN_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=TZA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=RWA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=SSD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=COD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_path(data=UGA_border, aes(y=lat, x=long, group = group), color = "grey50", size = .5) +
  geom_sf(data = roads_UG_2012, color = "grey30", size = .25) +
  geom_point(data = pca_output_ug_2016, aes(x = long, y = lat), size = .25) +
  geom_sf(data = ug_grid_2016_shape, size = .05) + 
  geom_sf(data = ug_grid_2013_shape, fill = "red", color = "red", size = .05) +
  geom_sf(data = ug_grid_2010_shape, size = .1) +
  coord_sf(xlim = c(29.55, 34.85), ylim = c(-1.3, 4.05)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))
EDFigure10
dev.off()

# Inset
pdf("../../figures/raw/EDFigure10_inset.pdf") 
EDFigure10_inset <- ggplot() +
  geom_polygon(data=UGA_border, aes(x=long, y=lat, group = group), fill = "white") +
  geom_polygon(data=ug_water_bodies, aes(x=long, y=lat, group = group), fill = "lightcyan") +
  geom_polygon(data=KEN_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=TZA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=RWA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=SSD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=COD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_path(data=UGA_border, aes(y=lat, x=long, group = group), color = "grey50", size = .5) +
  geom_sf(data = roads_UG_2012, color = "grey30", size = .25) +
  geom_point(data = pca_output_ug_2016, aes(x = long, y = lat), size = .25) +
  coord_sf(xlim = c(29.55, 34.85), ylim = c(-1.3, 4.05)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))
EDFigure10_inset
dev.off()


