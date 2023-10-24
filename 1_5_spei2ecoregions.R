#####################################################################
## File to compute the mean of grid points from a .nc file
## and associate it with a region given by a shapefile.
#####################################################################

# Clear workspace
rm(list = ls())
graphics.off()
gc()

##########################
## LOAD LIBRARIES ##
##########################

library(ncdf4)
library(rgdal)
library(raster)
library(SPEI)

######################
## SET PARAMETERS ##
######################

where <- ifelse(Sys.info()['sysname'] == 'Darwin', 'mac', 'onfire') # Automatic machine detection

dir_settings <- list(
  mac = list(
    dir_out = '~/Dropbox/model/out_v2/',
    dir_shp = '~/Documents/virtualBox/temiav/Ecoregions2017/',
    dir_fwi= '~/Documents/dati/obs/ERA5/FWI/',
    dir_spei='~/Documents/dati/obs/ERA5/ERA5-sergio/',
    dir_prec='~/Documents/dati/obs/MSWEP/'
  ),
  onfire = list(
    dir_out = '~/Dropbox/model/out_v2/',
    dir_shp = '/diskonfire/shapefiles/Ecoregions2017/',
    dir_fwi= '/diskonfire/ERA5/FWI/',
    dir_spei='/diskonfire/ERA5-sergio/',
    dir_prec='/diskonfire/DROP/MSWEP/'
  )
)

dirs <- dir_settings[[where]]

#######################
## LOAD FILES ##
#######################

# Load shapefile
shp <- readOGR(file.path(dirs$dir_shp, "Ecoregions2017_repaired.shp"))
shp_ll <- spTransform(shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Function to process data
process_data <- function(filepath, shp) {
  data_nc = brick(filepath)
  data_shp = extract(data_nc, shp, na.rm = TRUE, fun = mean)
  return(t(data_shp))
}

# Process precipitation data
pre_reg <- process_data(file.path(dirs$dir_prec, 'PREC-MSWEP-1998-2021_05.nc'), shp_ll)
save(pre_reg, file = file.path(dirs$dir_out, 'pre-mswep-1998-2021-ecoregions.Rdata'))

# Process PET data
pet_reg <- process_data(file.path(dirs$dir_spei, 'eto_1998_2021.nc'), shp_ll)
save(pet_reg, file = file.path(dirs$dir_out, 'pet-era5-1998-2021-ecoregions.Rdata'))

# Print summaries
print(summary(as.vector(pre_reg)))
print(summary(as.vector(pet_reg)))

## Compute SPEI
compute_spei <- function(pre_data, pet_data, months) {
  result <- pre_data * NA
  for (i in 1:dim(pre_data)[2]) {
    if (!is.na(sum(pre_data[, i]))) {
      spei_data <- spei(pre_data[, i] - pet_data[, i], months, na.rm = TRUE)$fitted
      result[, i] = spei_data
    }
  }
  return(result)
}

spei3 <- compute_spei(pre_reg, pet_reg, 3)
spei6 <- compute_spei(pre_reg, pet_reg, 6)
spei12 <- compute_spei(pre_reg, pet_reg, 12)

# Adjust data and save
adjust_and_save <- function(data, filepath) {
  data[is.infinite(data)] <- NA
  data <- data[13:nrow(data),]
  save(data, file = filepath)
}

adjust_and_save(spei3, file.path(dirs$dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
adjust_and_save(spei6, file.path(dirs$dir_out, "spei6-1999-2021-ecoregions-mswep-era5.RData"))
adjust_and_save(spei12, file.path(dirs$dir_out, "spei12-1999-2021-ecoregions-mswep-era5.RData"))
