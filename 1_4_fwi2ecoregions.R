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

# ---------------
# Fixed Parameters
# ---------------

where <- ifelse(Sys.info()['sysname'] == 'Darwin', 'mac', 'onfire') # Automatic machine detection

dir_base <- ifelse(where == 'mac', "~/Documents/", "/diskonfire/")

dir_out  <- '~/Dropbox/model/out_v2/'
dir_shp  <- file.path(dir_base, "virtualBox/temiav/Ecoregions2017/")
dir_fwi  <- file.path(dir_base, "dati/obs/ERA5/FWI/")
dir_spei <- file.path(dir_base, "dati/obs/ERA5/ERA5-sergio/")

# -------------------
# Load Shapefile Data
# -------------------

# ecoregion
file_shp <- file.path(dir_shp, "Ecoregions2017_repaired.shp")
eco <- st_read(file_shp) %>% st_make_valid()

# -----------------
# Load Raster Data
# -----------------

file_nc <- file.path(dir_fwi, 'FWI-2000-2021_05.nc')
r <- rast(file_nc)

# Extract mean BA for each ecoregion
fwi_reg <- terra::extract(r, eco, fun = "mean")
fwi_reg <- t(fwi_reg[-1, ])

save(fwi_reg, file = file.path(dir_out, 'fwi-2000-2021-ecoregions.Rdata'))

# ---------------------
# Compute STI Values
# ---------------------

STI <- function(tmean, sc) {
  nm <- length(tmean)
  
  # Get the accumulated data for the time scale sc
  A1 <- matrix(NA, (length(tmean) - sc + 1), sc)
  for (i in 1:sc) {
    A1[, i] <- tmean[i:(length(tmean) - sc + i)]
  }
  Y <- apply(A1, 1, mean, na.rm = TRUE)
  
  c(matrix(NA, sc - 1, 1), Y)
}

compute_fwi <- function(fwi_reg, sc) {
  result <- fwi_reg * NA
  for (i in 1:dim(fwi_reg)[2]) {
    if (!is.na(sum(fwi_reg[, i]))) {
      result[, i] <- STI(fwi_reg[, i], sc)
    }
  }
  result
}

fwi3 <- compute_fwi(fwi_reg, 3)
fwi6 <- compute_fwi(fwi_reg, 6)
fwi12 <- compute_fwi(fwi_reg, 12)

# Save computed data
save(fwi3, file = file.path(dir_out, "sfwi3-2000-2021-ecoregions.Rdata"))
save(fwi6, file = file.path(dir_out, "sfwi6-2000-2021-ecoregions.Rdata"))
save(fwi12, file = file.path(dir_out, "sfwi12-2000-2021-ecoregions.Rdata"))
