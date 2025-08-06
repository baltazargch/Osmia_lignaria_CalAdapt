library(sf)
library(biomod2)
library(tidyverse)
library(terra)
library(ENMTools)
library(glue)
library(rgbif)
library(virtualspecies)
sf_use_s2(F)
source('R/udf_envSamp.R')

# Load occs for Osmia and and Target group
osli <- read_csv('inputs/records/cleaned_osmia_lignaria_records_august2025.csv')


tg_os <- read_csv('inputs/records/pres_abs_oslig_target_group.csv')

# Using NAD83 Equal Area Albers here. Later transform to 4269 for 
# CA_CC report
tg_os_sf <- st_as_sf(tg_os, coords=c('lon', 'lat'), crs=4326) %>% 
  st_transform(5070) %>% filter(class=='background')


osmia_occs_buff <- st_buffer(tg_os_sf, dist = 100000) %>% 
  st_union() 

study_area <- terra::rasterize(osmia_occs_buff %>% 
                                 st_transform(4326) %>% 
                                 vect(),
                               hist_rast[[1]],
                               field = 1) %>%
  terra::mask(hist_rast[[1]]) 

hist_rast  <- rast('outputs/bioclim_vars/biovars_2000_2014_historical_CNRM-ESM2-1.tif')
names(hist_rast)

names(hist_rast) <- c(
  "bio1"  = "Ann_Mean_Temp",
  "bio2"  = "Mean_Diurnal_Range",
  "bio3"  = "Isothermality",
  "bio4"  = "Temp_Seasonality",
  "bio5"  = "Max_Temp_Warm_Month",
  "bio6"  = "Min_Temp_Cold_Month",
  "bio7"  = "Ann_Temp_Range", 
  "bio8"  = "Mean_Temp_Wettest_Q", 
  "bio9"  = "Mean_Temp_Driest_Q",
  "bio10" = "Mean_Temp_Warm_Q",
  "bio11" = "Mean_Temp_Cold_Q",
  "bio12" = "Ann_Precip",
  "bio13" = "Precip_Wet_Month",
  "bio14" = "Precip_Dry_Month",
  "bio15" = "Precip_Season",
  "bio16" = "Precip_Wet_Q",
  "bio17" = "Precip_Dry_Q",
  "bio18" = "Precip_Warmest_Q",
  "bio19" = "Precip_Coldest_Q"
)

