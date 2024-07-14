# Script for exploring the dust metadata

# SETUP ####

## Packages ####
library(tidyverse)
library(phyloseq)
library(ggmap)
library(ggpp)
library(leaflet)
library(lmerTest)
library(htmltools)

## Extra functions ####
source("./R/googlemap_styling.R")

## Environment variables ####
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

## Load data ####
dat <- readRDS("./data/full_clean_metadata.RDS")

# PLAYING AROUND ####


## Canister position ####

# test effect of "position height"
mod1 <- lmer(data=dat,
             formula = its_reads ~ height + (1|site))
rsq::rsq(mod1)
summary(mod1)
mod1 <- lmer(data=dat,
             formula = ssu_reads ~ height + (1|site))
rsq::rsq(mod1)
summary(mod1)
mod1 <- lmer(data=dat,
             formula = ssu_asv_richness ~ height + (1|site))
rsq::rsq(mod1)
summary(mod1)

mod1 <- lmer(data=dat,
             formula = its_asv_richness ~ height + (1|site))
rsq::rsq(mod1)
summary(mod1)

## Site maps ####

# build map styling
mapstyle <- rjson::fromJSON(file = "./R/mapstyle3.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# Get googlemap
area <- 
  ggmap::get_googlemap(center = c(lon = dat$long_dd %>% mean(na.rm=TRUE), 
                                  lat = dat$lat_dd %>% mean(na.rm=TRUE)),
                       zoom = 3,
                       scale = 2,
                       style=mapstyle)

# simplify dat set to avoid redundancy
sites <- dat %>% 
  dplyr::select(long_dd,
                lat_dd,
                site,
                its_asv_richness,
                ssu_asv_richness,
                mean_canopy_wind_jun_nov_m_s,
                mean_temp_jun_nov_c,
                total_precip_jun_nov_mm) %>% 
  unique.data.frame() %>% 
  group_by(site) %>% 
  summarize(long_dd = mean(long_dd),
            lat_dd = mean(lat_dd),
            its_richness = mean(its_asv_richness,na.rm=TRUE),
            ssu_richness = mean(ssu_asv_richness,na.rm=TRUE),
            canopy_wind_jun_nov = mean(mean_canopy_wind_jun_nov_m_s,na.rm=TRUE),
            temp_jun_nov = mean(mean_temp_jun_nov_c,na.rm=TRUE),
            precip_jun_nov = mean(total_precip_jun_nov_mm,na.rm=TRUE))

# label nudge factors
sites$nudge_y <- -1
sites$nudge_y[which(sites$site == "CPER")] <- 0

### Static map ####
ggmap::ggmap(area) +
  geom_point(data=sites,aes(x=long_dd,y=lat_dd),color='white',size=1.5,alpha=.5) +
  geom_point(data=sites,aes(x=long_dd,y=lat_dd),color='darkorange',size=1) +
  geom_text_s(data=sites,color.target = 'all',
              aes(x=long_dd,y=lat_dd,label=site),
              nudge_x = 1,
              nudge_y = sites$nudge_y,
              color='white',size=2,
              point.padding = 0.4)


### Leaflet map ####

# custom color function
its_pal <- colorNumeric(palette = 'viridis',
                        domain = sites$its_richness)
pal <- its_pal(sites$its_richness)

lab <- paste0('<p>',sites$site,'<p></p>',
              '<p>',"canopy_wind: ",round(sites$canopy_wind_jun_nov,2),'<p></p>',
              '<p>',"temp: ",round(sites$temp_jun_nov,2),'<p></p>',
              '<p>',"precip: ",round(sites$precip_jun_nov,2),'<p></p>')

m = leaflet() %>% addTiles()
m  # a map with the default OSM tile layer

# map starting scope
m = m %>% setView(lng = dat$long_dd %>% mean(na.rm=TRUE), 
                  lat = dat$lat_dd %>% mean(na.rm=TRUE),
                  zoom=2)

# map with points added
m %>% 
  addCircleMarkers(lng=sites$long_dd,
                   lat=sites$lat_dd,
                   radius = 4,
                   color = pal,
                   fill = TRUE,
                   label = (as.list(lab) %>% map(HTML)),
                   opacity = 1) %>% 
  addLegend(position = "topright", pal = its_pal, values = sites$its_richness,
            title = "Mean ITS richness")




# POTENTIAL ANALYSES ####

# GDM for location/environment influence on community structure