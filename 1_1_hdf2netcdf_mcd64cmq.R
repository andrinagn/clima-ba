# Clear the workspace and loaded graphical plots
rm(list = ls())
graphics.off()
gc()

# Load necessary libraries
library(stars)
library(maptools)
library(ncdf4)

# Set directory and grid parameters
data(wrld_simpl)

where <- ifelse(Sys.info()['sysname'] == 'Darwin', 'mac', 'onfire') # Automatic machine detection
dir_settings <- list(
  mac = list(
    dir_fire = '/Users/marco/Documents/dati/modis/MCD64CMQ/hdf/'
  ),
  onfire = list(
    dir_out = '~/Dropbox/model/out_v2/',
    dir_shp = '/diskonfire/shapefiles/Ecoregions2017/',
    dir_fwi= '/diskonfire/ERA5/FWI/',
    dir_spei='/diskonfire/ERA5-sergio/',
    dir_fire = '/diskonfire/modis/MCD64CMQ/'
  )
)
dirs <- dir_settings[[where]]

XGRID = 0:1439
YGRID = 0:719
bin_size = 0.25
lat <- (90.0 - bin_size / 2) - bin_size * (YGRID)
lon <- (-180.0 + bin_size / 2) + bin_size * (XGRID)
lat <- rev(lat)

# List all HDF files
file_files <- Sys.glob(file.path(dirs$dir_fire, "*.hdf"))

# Iterate over each HDF file
for (filename in file_files) {
  print(filename)
  s# Extract year and month from the file name
  i1 = nchar(dirs$dir_fire) + 12
  year <- as.integer(substr(filename, i1, i1 + 3))
  day_of_year <- substr(filename, i1 + 4, i1 + 6)
  date <- as.Date(paste0(year, day_of_year), format = "%Y%j")
  month <- format(date, "%m")
  cat("Year:", year, "\nMonth:", month, "\n")
  
  # Read and process the HDF file
  r <- read_stars(gdal_subdatasets(filename)[[1]]) * 0.01
  cat("Total burned area is", format(sum(as.numeric(unlist(r))) / 1.0e6, nsmall = 2), "Mha\n")
  BA <- matrix(as.numeric(unlist(r)), ncol = length(lat), nrow = length(lon))
  
  # If the year is between 2001 and 2021, save as a NetCDF file
  if (year >= 2001 & year <= 2021) {
    londim <- ncdim_def("lon", "degrees_east", as.single(lon))
    latdim <- ncdim_def("lat", "degrees_north", as.single(lat))
    mean_def <- ncvar_def("BA", "ha", list(londim, latdim), 1e32, prec = "double")
    ncfname <- paste0(substr(filename, 1, nchar(filename)-3), "nc")
    
    ncout <- nc_create(ncfname, list(mean_def))
    ncvar_put(ncout, mean_def, BA)
    ncatt_put(ncout, "lon", "axis", "X")
    ncatt_put(ncout, "lat", "axis", "Y")
    nc_close(ncout)
    
    cat("NetCDF file created:", ncfname, "\n")
  }
}
