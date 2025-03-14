library(tidyverse)
library(ggmap)
library(maps)
library(mapdata)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(viridis)

locations <- read_csv("analyses/data/LocationData.csv") %>% 
  filter(state == "CO")

counties <- map_data("county") %>% 
  filter(region == "colorado")

map1 <- ggplot() + 
  geom_polygon(data = counties, aes(x = long, y = lat, group = group), fill = "snow1",
               color = "grey10", linewidth=0.35) + 
  theme_few() +
  #  coord_fixed(xlim = c(-106, -104), ylim = c(37,41), ratio = 1.2) +
  geom_jitter(data = locations, aes(x = long, y = lat, fill = street), 
              shape = 21, stroke = 0.25, size = 6, alpha = 0.8, width = 0.15,
              height = 0.15) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                               "coral2", "coral4", "darkorange","darkorange2",
                               "darkorange4")) +
  labs(fill = "Street Name") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        legend.text=element_text(family="Helvetica", size=10),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20))
map1
ggsave("analyses/plots/llocation_map1.pdf", units = "in", width = 5, height = 3)
ggsave("analyses/plots/llocation_map1.jpg", units = "in", width = 5, height = 3)

counties <- map_data("county") %>% 
  filter(region == "colorado") %>% 
  filter(subregion %in% c("larimer", "jackson", "grand", "boulder", "weld"))

map2 <- ggplot() + 
  geom_polygon(data = counties, aes(x = long, y = lat, group = group), fill = "snow1",
               color = "grey10", linewidth=0.35) + 
  theme_few() +
#  coord_fixed(xlim = c(-106, -104), ylim = c(37,41), ratio = 1.2) +
  geom_jitter(data = locations, aes(x = long, y = lat, fill = street), 
              shape = 21, stroke = 0.25, size = 6, alpha = 0.8, width = 0.15,
              height = 0.15) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                                 "darkorange4")) +
  labs(fill = "Street Name") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        legend.text=element_text(family="Helvetica", size=10),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20))
map2
ggsave("analyses/plots/llocation_map2.pdf", units = "in", width = 5, height = 3)
ggsave("analyses/plots/llocation_map2.jpg", units = "in", width = 5, height = 3)

locations2 <- read_csv("analyses/data/LocationData.csv") %>% 

usa <- map_data("state")

map3 <- ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "snow1",
               color = "grey10", linewidth=0.35) + 
  theme_few() +
  #  coord_fixed(xlim = c(-106, -104), ylim = c(37,41), ratio = 1.2) +
  geom_jitter(data = locations2, aes(x = long, y = lat, fill = street), 
              shape = 21, stroke = 0.25, size = 6, alpha = 0.8, width = 0.4,
              height = 0.4) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "#00FFFF", "chocolate4", "coral",
                                 "coral2", "coral4", "darkorange","darkorange2",
                               "#006666", "darkorange4", "#00CCCC")) +
  labs(fill = "Street Name") +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        legend.text=element_text(family="Helvetica", size=10),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20))
map3
ggsave("analyses/plots/location_map3.pdf", units = "in", width = 6, height = 3)
ggsave("analyses/plots/location_map3.jpg", units = "in", width = 6, height = 3)
