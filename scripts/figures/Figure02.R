# This code creates Figure 2 a-c. 

rm(list=ls())

library(ggplot2)
library(raster)
library(wesanderson)
library(rgdal)

setwd("data/figure_and_input_data/")

r_2006 <- raster("raster_2006_ug.tif")
pr_2006 <- projectRaster(r_2006,crs="+proj=longlat +datum=NAD27")
br_2006 <- unique(quantile(pr_2006,probs=seq(0,1,0.01)))
pal <- wes_palette("Zissou1", length(br_2006), type = "continuous")

r_2016 <- raster("raster_2016_ug.tif")
pr_2016 <- projectRaster(r_2016,crs="+proj=longlat +datum=NAD27")
br_2016 <- unique(quantile(pr_2016,probs=seq(0,1,0.01)))

UGA_border <- readOGR(dsn = "gadm36_UGA_0.shp", layer = "gadm36_UGA_0", stringsAsFactors = F)
UGA_border_1 <- readOGR(dsn = "gadm36_UGA_1.shp", layer = "gadm36_UGA_1", stringsAsFactors = F)
ug_water_bodies <- readOGR(dsn = "Ug_Waterbodies.shp", layer = "Ug_Waterbodies", stringsAsFactors = F)
COD_border <- readOGR(dsn = "gadm36_COD_0.shp", layer = "gadm36_COD_0", stringsAsFactors = F)
KEN_border <- readOGR(dsn = "gadm36_KEN_0.shp", layer = "gadm36_KEN_0", stringsAsFactors = F)
TZA_border <- readOGR(dsn = "gadm36_TZA_0.shp", layer = "gadm36_TZA_0", stringsAsFactors = F)
RWA_border <- readOGR(dsn = "gadm36_RWA_0.shp", layer = "gadm36_RWA_0", stringsAsFactors = F)
SSD_border <- readOGR(dsn = "gadm36_SSD_0.shp", layer = "gadm36_SSD_0", stringsAsFactors = F)

WI_community_change <- read.csv("WI_community_change.csv", check.names=FALSE)

# Figure 2a
pdf("../../figures/raw/WI_2006_raster_map.pdf") 
WI_2006_raster_map <- plot(pr_2006, col=pal, axes=F) 
plot(ug_water_bodies,add=TRUE,col="lightcyan",border="transparent")
plot(KEN_border,add=TRUE,col="white",border="white")
plot(TZA_border,add=TRUE,col="white",border="white")
plot(RWA_border,add=TRUE,col="white",border="white")
plot(SSD_border,add=TRUE,col="white",border="white")
plot(COD_border,add=TRUE,col="white",border="white")
WI_2006_raster_map
dev.off()

# Figure 2b
pdf("../../figures/raw/WI_2016_raster_map.pdf") 
WI_2016_raster_map <-plot(pr_2016,col=pal,axes=F) 
plot(ug_water_bodies,add=TRUE,col="lightcyan",border="transparent")
plot(KEN_border,add=TRUE,col="white",border="white")
plot(TZA_border,add=TRUE,col="white",border="white")
plot(RWA_border,add=TRUE,col="white",border="white")
plot(SSD_border,add=TRUE,col="white",border="white")
plot(COD_border,add=TRUE,col="white",border="white")
WI_2016_raster_map
dev.off()

# Figure 2c
pdf("../../figures/raw/raster_WI_change.pdf") 
raster_WI_change <- ggplot() +
  geom_polygon(data=UGA_border, aes(x=long, y=lat, group = group), fill = "grey70") +
  geom_polygon(data=ug_water_bodies, aes(x=long, y=lat, group = group), fill = "lightcyan") +
  geom_polygon(data=KEN_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=TZA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=RWA_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=SSD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_polygon(data=COD_border, aes(y=lat, x=long, group = group), fill = "white") +
  geom_path(data=UGA_border, aes(y=lat, x=long, group = group), size = .25) +
  geom_point(data = WI_community_change %>%
               arrange(change),
               aes(x = long, y = lat, color = -change), size = .65, show.legend = FALSE) +
  scale_color_gradientn(colours = rainbow(5)) +
  geom_path(data=UGA_border_1, aes(y=lat, x=long, group = group), size = .1, color = "black") +
  coord_sf(xlim = c(29.55, 34.85), ylim = c(-1.3, 4.05)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
raster_WI_change
dev.off()



 


