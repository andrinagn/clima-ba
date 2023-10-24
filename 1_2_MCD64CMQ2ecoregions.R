####################################################################################
## File to average grid points from a .nc file and associate it with a region 
## defined by a shapefile.
####################################################################################

# Clear the workspace
rm(list = ls())
graphics.off()
gc()

#######################
## PACKAGES & SOURCE ##
#######################

library(rgdal)        # For reading shapefiles
library(raster)       # For associating nc files with shapefiles


######################
## FIXED PARAMETERS ##
######################

# Define the directory paths based on the computer being used
where = 'mac'

if (where == 'mac') {
  dir_out = '~/Dropbox/model/out_v2/'
  dir_shp = '/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'
  dir_fire = '/Users/marco/Documents/dati/modis/MCD64CMQ/hdf'
} else if (where == 'onfire') {
  dir_out = '~/Dropbox/model/out_v2/'
  dir_shp = '/diskonfire/shapefiles/Ecoregions2017/'
  dir_fire = '/diskonfire/modis/MCD64CMQ'
}

years = 2001:2021
months = 1:12

###################
## LOAD FILES ##
###################

# Load shapefile data
file_shp = file.path(dir_shp, "Ecoregions2017_repaired.shp")
shp = readOGR(file_shp)
shp_ll = spTransform(shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


# Gather all .nc files from the fire directory
file_files = Sys.glob(file.path(dir_fire, "*.nc"))

# Matrix to store data
obs_reg = matrix(data = NA, nrow = length(years) * 12, ncol = length(shp_ll))

k = 0
for (i in 1:length(file_files)) {
  filename = file_files[i]
  print(filename)
  
  # Extract year and day of year from filename
  i1 = nchar(dir_fire) + 12
  year = substr(filename, i1, i1 + 3)
  day_of_year = substr(filename, i1 + 4, i1 + 6)
  
  # Convert day of year to month
  date = as.Date(paste0(year, day_of_year), format = "%Y%j")
  month = format(date, "%m")
  cat("Year:", year, "\n")
  cat("Month:", month, "\n")
  
  # Process files for specified year range
  if (year >= 2001 & year <= 2021) {
    k = k + 1
    f_pr = file.path(filename)
    pr_nc = brick(f_pr)
    pr_shp = extract(pr_nc, shp_ll, na.rm = T, fun = sum)
    obs_reg[k,] = pr_shp
  }
}

# Save results
save(obs_reg, file = file.path(dir_out, paste('MCD64CMQ-2001-2021.Rdata', "", sep = "")))
