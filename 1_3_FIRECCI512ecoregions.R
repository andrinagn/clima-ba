#####################################################################
# Script to calculate the average of grid points from an .nc file 
# and associate it with a region given by a shapefile.
#####################################################################

# --------------------------
# Clear Workspace & Clean-Up
# --------------------------

rm(list = ls())
graphics.off()
gc()

# -----------------------
# Load Required Libraries
# -----------------------

library(terra)       # Operations with raster and vector spatial data
library(sf)          # Simple Features: handling spatial vector data
library(tidyverse)   # Data manipulation
library(rmapshaper)  # Shapefile operations

# Output directory
dir_out = '~/Dropbox/model/out_ecoregions/'
# Directories for shapefiles and fire data
dir_shp = '/diskonfire/shapefiles/Ecoregions2017/'
dir_fire= '/diskonfire/FireCCI51'

# -------------------
# Load Shapefile Data
# -------------------

# ecoregion
file_shp <- file.path(dir_shp, "Ecoregions2017_repaired.shp")
eco <- st_read(file_shp) %>% st_make_valid()

# Load Raster Data
file_nc <- file.path(dir_fire, 'all_years_natural_fires_combined.nc')
r <- rast(file_nc)
# Extract BA for each ecoregion
obs_reg <- terra::extract(r, eco, fun = "sum")
obs_reg=obs_reg[,-1 ] #the first column is the ID, a progressive number
obs_reg <- t(obs_reg)
# Save computed data
save(obs_reg, file = file.path(dir_out, 'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020-natural-fires.Rdata'))


# Load Raster Data
file_nc <- file.path(dir_fire, 'all_years_burned_area_combined.nc')
r <- rast(file_nc)
# Extract BA for each ecoregion
obs_reg <- terra::extract(r, eco, fun = "sum")
obs_reg=obs_reg[,-1 ] #the first column is the ID, a progressive number
obs_reg <- t(obs_reg)
# Save computed data
save(obs_reg, file = file.path(dir_out, 'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020.Rdata'))
