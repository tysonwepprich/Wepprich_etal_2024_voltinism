
library(sf)
library(tidyverse)

# Site Map
# Supplemental Figure 1
region <- st_read("data/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp")
oh <- region[region$name_en == "Ohio", ]
oh <- st_transform(oh, 4326)
# reg.points = fortify(oh, region="name_en")


allpops <- readRDS("data/modeling_data.rds")

# used sites
test <- allpops %>% 
  group_by(SiteID, lat, lon) %>% 
  summarise(obsYearbySite = length(unique(Year))) %>% 
  distinct()
# duplicates make figure weird
test <- test %>% 
  group_by(lat, lon) %>% 
  arrange(-obsYearbySite) %>% 
  mutate(keep = c("yes", rep("no", n() - 1))) %>% 
  arrange(lat, lon)


sites <- read.csv("data/OHsites2023update.txt") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0")) %>% 
  left_join(test) %>% 
  group_by(lat, lon) %>% 
  arrange(-obsYearbySite) %>% 
  mutate(keep = c("yes", rep("no", n() - 1))) %>%
  filter(keep == "yes") %>%
  filter(obsYearbySite >= 5) %>% 
  mutate(Years_monitored = case_when(is.na(obsYearbySite) == TRUE ~ "less than 5",
                                     obsYearbySite >= 5 & obsYearbySite <11 ~ "5-10",
                                     obsYearbySite >= 11 & obsYearbySite <20 ~ "11-19",
                                     obsYearbySite >= 20 ~ "20 or more"))

sites$Years_monitored = factor(sites$Years_monitored, levels = unique(sites$Years_monitored)[4:1])


cities <- data.frame(name = c("Columbus", "Cincinnati", "Cleveland", "Dayton", "Toledo", "Akron"),
                     lat = c(39.985, 39.14, 41.478, 39.777, 41.664, 41.08176),
                     lon = c(-82.985, -84.506, -81.679, -84.2, -83.582, -81.51145))


library(viridis)

mapplt <- ggplot(sites, aes(x = lon, y = lat, color = Years_monitored)) +
  geom_point(size = 5, alpha = .5) +
  scale_color_viridis(name = "Years monitored", begin = .8, end = 0, discrete = TRUE) +
  geom_point(data = cities, aes(x = lon, y = lat), color = "black", size = 3, inherit.aes = FALSE) +
  # geom_text(data = cities, aes(x = lon, y = lat, label=name),hjust=1, vjust=0, inherit.aes = FALSE) +
  theme_minimal(base_size = 24) +
  theme(legend.position = c(.85, .15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_sf(data = oh, fill = NA, color = "black", inherit.aes = FALSE, size = 1) +
  coord_sf(crs=st_crs(4326)) +
  # coord_fixed(1.3, xlim = c(-85.15, -80.2), ylim = c(38.35, 42.1), expand = FALSE, clip = "on") +
  xlab("Longitude") +
  ylab("Latitude")
mapplt

# Fig. S1 (cities added in Powerpoint) ----
ggsave(plot = mapplt, filename =  "map.png", device = "png", width = 10, height = 10, units = "in")
ggsave(filename = "figures/figS1.tiff", plot = mapplt, device = "tiff", dpi = "retina", 
       width = 10, height = 10, units = "in", compression = "lzw")

