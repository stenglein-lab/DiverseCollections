library(tidyverse)
library(ggmap)
library(ggspatial)

# Create a map of sample locations

# read in sample location information
all_locations <- read_csv("analyses/data/LocationData.csv") 
co_locations <- all_locations %>% filter(state == "CO")

# get stamen map data for Fort Collins inset map

# define lat/long borders for Fort Collins maps
foco_borders <- c(bottom  = 40.35,
                  top     = 40.65,
                  left    = -105.3,
                  right   = -104.9)

# download data from stadia maps
# note that these maps are permitted to be used for non-commerical academic purposes, including in papers
# see:
# https://stadiamaps.com/faqs/
# need to create an account with stadia and create an API key 
# and register the API key with register_stadiamaps()
foco_terrain_map <- get_stadiamap(foco_borders, zoom = 11, maptype = "stamen_terrain")

# make the Fort Collins inset map
foco_map <- 
  ggmap(foco_terrain_map) +
  geom_jitter(data = co_locations, aes(x = long, y = lat, fill = street), 
              shape = 21, stroke = 0.25, size = 4, alpha = 0.8, width = 0.005,
              height = 0.005) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "chocolate4", "coral",
                               "coral2", "coral4", "darkorange","darkorange2",
                               "darkorange4")) +
  labs(fill = "Street Name") +
  annotation_scale() +
  coord_sf(crs = 4326) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        legend.position = "none",
        legend.text=element_text(family="Helvetica", size=10),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20))

foco_map

# save as TIFF
# this tiff image will be inserted into Figure 1 using Affinity Designer
ggsave("analyses/plots/map_images/foco_terrain_map.tiff", width=2.75, height=2.75, units="in")

# now create the whole-USA map

# get the stadia/stamen map for whole of USA
usa_background <- c("left" = -125, "bottom" = 25, "right" = -67, "top" = 49.5) %>%
  get_stadiamap(zoom = 6, maptype = "stamen_terrain_background") 

# this defines USA state boundaries
usa_states <- map_data("state")

# make the USA map
usa_terrain_map <-
  ggmap(usa_background, darken = c(0.33, "grey95")) + 
  geom_polygon(data = usa_states, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey20", linewidth=0.2, alpha=1) +
  geom_point(data = all_locations, aes(x = long, y = lat, fill = street), 
              shape = 21, stroke = 0.25, size = 5, alpha = 0.8) +
  scale_fill_manual(values = c("chocolate", "chocolate1", "#00FFFF", "chocolate4", "coral",
                               "coral2", "coral4", "darkorange","darkorange2",
                               "#006666", "darkorange4", "#00CCCC")) +
  labs(fill = "Street Name") +
  coord_map() +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        legend.position = "none",
        legend.text=element_text(family="Helvetica", size=10),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 20))

usa_terrain_map

# save TIFF, which we will insert into Figure 1 using Affinity Designer
ggsave("analyses/plots/map_images/USA_terrain_map_few.tiff", width=5, height=4, units="in")


