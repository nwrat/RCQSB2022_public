
setwd(dirname(getActiveDocumentContext()$path))   
gc(); rm(list = ls())
require(pacman)
p_load(ggplot2, tidyverse, rstudioapi, wesanderson, fixest)
​
​
df_treat <- read.csv("tc_2km_ug.csv")
df_dhs <- read.csv("WI_raw_survey_data.csv")
df_dhs <- df_dhs[grep("UG", df_dhs$DHSID), ]
df_using <- left_join(df_treat, df_dhs, by = "DHSID")
inc_years <- c(2006, 2009, 2011, 2014, 2016)
df_using <- df_using[df_using$year %in% inc_years, ]
df_using <- df_using[, !colnames(df_using) %in% c("lon", "lat.y", "adm1fips", "long",
                                                  "adm1dhs", "lat.x", "row_num", "hhid", 
                                                  "country_year", "cluster_id", "urban")]
​
treated_f <- function(x) { 
  if(x == 2006 | x == 2009) y <- "pre"
  if(x == 2011 | x == 2014 | x == 2016) y <- "post"
  return(y)
}
​
df_using$treat_period <- factor(sapply(df_using$year, treated_f), levels = c("pre", "post"))
​
treat_vars <- c("hv207", "hv208", "hv209", "hv210", "hv211", "hv212", "hv221", 
                "roof_qual", "floor_qual", "wall_qual", "toilet_qual", "water_qual")
treat_names <- c("radio", "tv", "refrig", "bicycle", "motorcycle", "car", "telephone", 
                 "roof", "floor", "walls", "toilet", "water")
treat_labels <- c("Radio", "TV", "Refrigerator", "Bicycle", "Motorcycle", "Car", 
                  "Telephone", "Roof", "Floor", "Walls", "Toilet", "Water Source")
​
## Take averages at the village level from DHS data
​
for(i in 1:length(treat_vars)){
  df_using[, treat_names[i]] <- scale(df_using[, treat_vars[i]])}
​
id_names <- c("DHSID", "year", "treated", "treat_period")
df_avg <- df_using %>% 
  group_by(DHSID, treated, treat_period) %>% 
  summarise_all(mean, na.rm = T)
​
​
results <- list()
df_results <- c()
for(i in 1:length(treat_names)){
  form <- as.formula(paste0(treat_names[i], " ~ treated*treat_period"))
  df_results <- rbind(df_results, summary(lm(form, data = df_using))$coefficients[4, ])
}
#df_results <- rbind(df_results, summary(feols(form, data = df_avg))$coeftable)
colnames(df_results) <- c("est", "std_error", "t", "p")
rownames(df_results) <- treat_labels
df_results <- as.data.frame(df_results)
df_results$ci_upper <- df_results$est + qnorm(0.975)*df_results$std_error
df_results$ci_lower <- df_results$est - qnorm(0.975)*df_results$std_error
df_results$varname <- rownames(df_results)
​
vars_rm <- c("Bicycle")
df_plot <- df_results[!df_results$varname %in% vars_rm, ]
treat_order <- rev(c("Floor", "Roof", "Walls", "Water Source", "Toilet", 
                     "Refrigerator", "Radio", "TV",
                     "Telephone", "Car", "Motorcycle"))
wes_colors <- wes_palette("Darjeeling1", 5)
df_plot$varname <- factor(df_plot$varname, levels = treat_order)
df_plot <- df_plot[rev(order(df_plot$varname)), ]
df_plot$treat_family <- factor(c(rep("Structure", 5), rep("Appliance",3), 
                                 rep("Object", 3)), 
                               levels = c("Structure", "Appliance", "Object"))
df_plot$treat_colors <- c(rep(wes_colors[1], 5), rep(wes_colors[2], 3), 
                          rep(wes_colors[4], 3))
​
ggplot(df_plot, aes(x = est, y = varname, color = treat_family)) + 
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, color = treat_family)) + 
  geom_point() + 
  geom_vline(xintercept = 0, col = "blue", lty = "dashed") + 
  xlab("Estimated Effect of Electrification 
       (Std Deviation Units)") + 
  ylab("Index Component") + 
  xlim(-0.45, 0.5) + theme_minimal() + 
  scale_color_manual(values = wes_colors[c(1,2,4)], 
                     name = "Item Type") 
​
ggsave("per_item_results.tiff", device = "tiff", dpi = 700, width = 7, height = 4)
