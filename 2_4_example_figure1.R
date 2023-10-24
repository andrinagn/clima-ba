rm(list = ls())
graphics.off()
gc()


library(RColorBrewer)
library(terra)
library(sf)
library(tidyverse)
library(dplyr)
library(rmapshaper)
library(rnaturalearth)
library(pals)
library(magrittr)
library(sf)
library(R.matlab)

### fix parameters ANDRINA VENTISCA
# dir_out = '/home/andrina/Desktop/climate-fire/out_v2/'
# dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
# dir_shp='/home/andrina/Desktop/climate-fire/datos/shp/'

### fix parameters A
dir_modis = '~/Dropbox/model/datos/'
dir_out = '~/Dropbox/model/out_ecoregions/'
dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
dir_shp = '/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'

years = 2001:2021

### version
transf = "log" #"log" "radq"

if (transf == "log") {
  list_version = c(
    'spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018',
    'spei_spei_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018'
  )
} else if (transf == "radq") {
  list_version = c(
    'spei_sfwi_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018',
    'spei_spei_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018'
  )
}
list_version_names = c('spei-sfwi', 'spei-spei')
list_models = c('spei-sfwi', 'spei-spei', 'ANTspei', 'CONsfwi', 'CONspei')

load(paste0(dir_out,
            "mask_ecoregions_season_ja_2018.RData"))

nreg = dim(mask)

all_r2 <- array(NA, dim = c(nreg, length(list_version)))
all_sfwi <- array(NA, dim = c(nreg, length(list_version)))
all_spei <- array(NA, dim = c(nreg, length(list_version)))
all_m_sfwi <- array(NA, dim = c(nreg, length(list_version)))
all_m_spei <- array(NA, dim = c(nreg, length(list_version)))
all_t_sfwi <- array(NA, dim = c(nreg, length(list_version)))
all_t_spei <- array(NA, dim = c(nreg, length(list_version)))

for (k in 1:length(list_version)) {
  version = list_version[k]
  load(paste0(dir_out, "corr_reconstruction_", version, ".RData"))
  load(paste0(dir_out, "sfwi_coef_reconstruction_", version, ".RData"))
  load(paste0(dir_out, "spei_coef_reconstruction_", version, ".RData"))
  
  load(paste0(dir_out, "best_m_SPI_fin_", version, ".RData")) #c(-14,-2)
  load(paste0(dir_out, "best_t_SPI_fin_", version, ".RData"))
  load(paste0(dir_out, "best_m_sfwi_fin_", version, ".RData")) #c(-12, 0)
  load(paste0(dir_out, "best_t_sfwi_fin_", version, ".RData"))
  # best_m_sfwi_fin
  # best_t_sfwi_fin
  # best_m_SPI_fin
  # best_t_SPI_fin
  
  all_r2[, k] = corr * corr
  all_spei[, k] = spei_coef
  all_sfwi[, k] = sfwi_coef
  all_t_spei[, k] = best_t_SPI_fin
  all_t_sfwi[, k] = best_t_sfwi_fin
  all_m_spei[, k] = best_m_SPI_fin
  all_m_sfwi[, k] = best_m_sfwi_fin
}


best_mod <- array(NA, dim = c(nreg))
best_r2 <- array(NA, dim = c(nreg))
best_m_spei <- array(NA, dim = c(nreg))
best_m_sfwi <- array(NA, dim = c(nreg))
best_t_spei <- array(NA, dim = c(nreg))
best_t_sfwi <- array(NA, dim = c(nreg))
best_spei_coef <- array(NA, dim = c(nreg))
best_sfwi_coef <- array(NA, dim = c(nreg))

for (ireg in 1:nreg) {
  if (length(which(is.na(all_r2[ireg, ]))) < length(list_version)) {
    # if (length(which(is.na(all_r2[i,j,])))==1) {cc}
    aux = which(all_r2[ireg, ] == max(all_r2[ireg, ], na.rm = T))
    if (length(aux) == 2) {
      aux = aux[1]
    }
    tmp = gregexpr(pattern = '-', list_version_names[aux])
    best_r2[ireg] = all_r2[ireg, aux]
    # list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')
    if (is.na(all_spei[ireg, aux]) & !is.na(all_sfwi[ireg, aux])) {
      aux2 = substr(list_version_names[aux],
                    as.numeric(tmp[[1]]) + 1,
                    nchar(list_version_names[aux]))
      aux2 = paste0('^CON', aux2, '$')
      best_m_sfwi[ireg] <- all_m_sfwi[ireg, aux]
      best_t_sfwi[ireg] <- all_t_sfwi[ireg, aux]
      best_sfwi_coef[ireg] <- all_sfwi[ireg, aux]
    }
    if (!is.na(all_spei[ireg, aux]) & is.na(all_sfwi[ireg, aux])) {
      aux2 = substr(list_version_names[aux], 1, as.numeric(tmp[[1]]) - 1)
      aux2 = paste0('^ANT', aux2, '$')
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_spei_coef[ireg] <- all_spei[ireg, aux]
    }
    if (!is.na(all_spei[ireg, aux]) &
        !is.na(all_sfwi[ireg, aux])) {
      aux2 = list_version_names[aux]
      best_m_sfwi[ireg] <- all_m_sfwi[ireg, aux]
      best_t_sfwi[ireg] <- all_t_sfwi[ireg, aux]
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_sfwi_coef[ireg] <- all_sfwi[ireg, aux]
      best_spei_coef[ireg] <- all_spei[ireg, aux]
    }
    
    best_mod[ireg] = grep(aux2, list_models)
    if (best_mod[ireg] == 2 ||
        best_mod[ireg] == 5) {
      best_sfwi_coef[ireg] = -best_sfwi_coef[ireg]
    }
  }
}



summary(as.vector(best_r2))

print(
  paste0(
    "Variance explained (spatial median) = ",
    round(median(best_r2, na.rm = TRUE), digits = 2),
    "\n",
    "Percentage valid models = ",
    round(100 * length(which(!is.na(
      best_r2
    ))) / length(which(mask == 1)), digits = 0),
    '%'
  )
)


### load data to analyse the model

## BA obs_reg
load(paste0(dir_modis,
            paste('MCD64CMQ-2001-2021.Rdata')))

nreg = dim(obs_reg)[2]
## fire season
# step 0: only cells with annual burned area series (BA>=2) that have at least 2 years were considered for further analysis.”).
BAy = array(0, dim = c(nreg, length(years)))
for (ireg in 1:nreg) {
  for (iyear in 1:length(years)) {
    i1 = (iyear - 1) * 12 + 1
    i2 = (iyear - 1) * 12 + 12
    BAy[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE) #*inout[i,j]
  }
}

mask0 = array(0, dim = c(nreg))
for (ireg in 1:nreg) {
  if (length(which(BAy[ireg,] > 0)) >= 2) {
    mask0[ireg] = 1
  }
}


firstmonth = readMat(file.path(dir_out, 'firstmonth_reg_modis.mat'))
endmonth = readMat(file.path(dir_out, 'endmonth_reg_modis.mat'))

firstmonth2 = firstmonth$firstmonth
endmonth2 = endmonth$endmonth

## aggregate BA over the fire season
BA_esa_y = array(NA, dim = c(nreg, length(years)))
for (ireg in 1:nreg) {
  if (mask0[ireg] == 1) {
    for (iyear in 1:length(years)) {
      if (firstmonth2[ireg] >= 1) {
        i1 = (iyear - 1) * 12 + firstmonth2[ireg]
        i2 = (iyear - 1) * 12 + endmonth2[ireg]
        BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
      } else
        if (iyear == 1) {
          next
        } else {
          i1 = (iyear - 1) * 12 + firstmonth2[ireg]
          i2 = (iyear - 1) * 12 + endmonth2[ireg]
          BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
        }
    }
  }
}

auxsum = apply(BA_esa_y, c(1), sum, na.rm = T)
mask = array(0, dim = c(nreg))
for (ireg in 1:nreg) {
  if (length(which(BA_esa_y[ireg, ] > 0)) >= 10 &&
      100 * auxsum[ireg] / sum(auxsum) > 0.001) {
    mask[ireg] = 1
  }
}


if (transf == "log") {
  BA = log(BA_esa_y + 1)
} else if (transf == "radq") {
  BA = sqrt(BA_esa_y)
}

## load climate
load(paste0(dir_out, "sfwi3-2000-2021-ecoregions.RData"))
load(paste0(dir_out, "sfwi6-2000-2021-ecoregions.RData"))
load(paste0(dir_out, "sfwi12-2000-2021-ecoregions.RData"))

load(paste0(dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
load(paste0(dir_out, "spei6-1999-2021-ecoregions-mswep-era.RData"))
load(paste0(dir_out, "spei12-1999-2021-ecoregions-mswep-era.RData"))
years_sfwi = 1999:2021 #this is different from the spei_sfwi model
iok_sfwi = match((years[1] - 1):years[length(years)], years_sfwi) #no need 1999
#sfwi coef is suppose to be >0, while the spi effect is supposet to be <0. to have the same considtions below, i changes the sign of the spi

spei3_cc = -spei3[(((iok_sfwi[1] - 1) * 12) + 1):(iok_sfwi[length(iok_sfwi)] *
                                                    12),]
spei6_cc = -spei6[(((iok_sfwi[1] - 1) * 12) + 1):(iok_sfwi[length(iok_sfwi)] *
                                                    12),]
spei12_cc = -spei12[(((iok_sfwi[1] - 1) * 12) + 1):(iok_sfwi[length(iok_sfwi)] *
                                                      12),]
## load SPEI AC
load(paste0(dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
load(paste0(dir_out, "spei6-1999-2021-ecoregions-mswep-era.RData"))
load(paste0(dir_out, "spei12-1999-2021-ecoregions-mswep-era.RData"))

sfwi = array(NA, dim = c(dim(fwi3)[1], nreg, 3))
sfwi[, , 1] = fwi3
sfwi[, , 2] = fwi6
sfwi[, , 3] = fwi12
rm(fwi3)
rm(fwi6)
rm(fwi12)

spei_cc = array(NA, dim = c(dim(sfwi)[1], nreg, 3))
spei_cc[, , 1] = spei3_cc
spei_cc[, , 2] = spei6_cc
spei_cc[, , 3] = spei12_cc
rm(spei3_cc)
rm(spei6_cc)
rm(spei12_cc)

spi = array(NA, dim = c(dim(spei3)[1], nreg, 3))
spi[, ,  1] = spei3
spi[, ,  2] = spei6
spi[, ,  3] = spei12
rm(spei3)
rm(spei6)
rm(spei12)

## the model

sfwi_coef = array(NA, dim = c(nreg))
spei_coef = array(NA, dim = c(nreg))
pval_sfwi_coef = array(NA, dim = c(nreg))
pval_spei_coef = array(NA, dim = c(nreg))


#check single regions

## climatologies
load(paste0(dir_out, "AI_ERA5_MSWEP_2001_2020_ecoregions.RData")) #AI
load(paste0(
  dir_out,
  "PREC_ERA5_MSWEP_2001_2020_ecoregions_annual.RData"
)) #AI
load(paste0(dir_out, "PET_ERA5_MSWEP_2001_2020_ecoregions_annual.RData")) #AI

plot(AI, best_spei_coef)
cor.test(AI, best_spei_coef)

file_shp = file.path(dir_shp, "Ecoregions2017_repaired.shp")
eco <- st_read(file_shp) %>% st_make_valid()

reg_test = which(eco$OBJECTID == 633)

eco_only <- filter(eco, mask==1)
reg_test_eco_only = which(eco_only$OBJECTID == 633)

wld <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") #https://github.com/ropensci/rnaturalearth
# Exclude Antarctica 
wld <- dplyr::filter(wld, name_long != "Antarctica")


g <- ggplot() +
  geom_sf(data = wld, fill = "lightgray", colour = "white", linewidth = 0.1) +  # Plotting the world (wld)
  geom_sf(data = eco_only, fill = "lightgray", colour = "white", linewidth = 0.1) +  # Plotting the entire eco dataset in light gray
  geom_sf(data = eco[reg_test, ], fill = "red", colour = "white", linewidth = 0.1) +  # Plotting only the reg_test region in red
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(legend.position = "none") +  # Remove legend since it's not relevant now
  guides(fill = T)  # Remove the fill guide

# Save the image
ggsave(paste0(dir_out, "ecoregions_4fig1.png"), g, height = 10, width = 15, unit = "in",
       bg = "white", dpi = 600)





best_mod[reg_test]

eco[reg_test,]

ireg=reg_test
if (length(which(is.na(all_r2[ireg, ]))) < length(list_version)) {
  # if (length(which(is.na(all_r2[i,j,])))==1) {cc}
  aux = which(all_r2[ireg, ] == max(all_r2[ireg, ], na.rm = T))
  if (length(aux) == 2) {
    aux = aux[1]
  }
  tmp = gregexpr(pattern = '-', list_version_names[aux])
  # best_r2[ireg] = all_r2[ireg, aux]
  # list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')
  if (is.na(all_spei[ireg, aux]) & !is.na(all_sfwi[ireg, aux])) {
    aux2 = substr(list_version_names[aux],
                  as.numeric(tmp[[1]]) + 1,
                  nchar(list_version_names[aux]))
    aux2 = paste0('^CON', aux2, '$')
  }
  if (!is.na(all_spei[ireg, aux]) & is.na(all_sfwi[ireg, aux])) {
    aux2 = substr(list_version_names[aux], 1, as.numeric(tmp[[1]]) - 1)
    aux2 = paste0('^ANT', aux2, '$')
  }
  if (!is.na(all_spei[ireg, aux]) &
      !is.na(all_sfwi[ireg, aux])) {
    aux2 = list_version_names[aux]
  }
  
  mydata_ba = data.frame("y" = BA[ireg,],
                         "x1" = (years))
  mydata_ba_det = data.frame("y" = scale(mydata_ba$y - predict(lm(y ~ x1 , data = mydata_ba), newdata =
                                                                 mydata_ba)))
  
  best_mod[ireg] = grep(aux2, list_models)
  # [1] "spei-sfwi" "spei-spei" "ANTspei"   "CONsfwi"   "CONspei"
  if (best_mod[ireg] == 1) {
    im_tmx = best_m_sfwi[ireg] + endmonth2[ireg,]
    if (im_tmx <= 0) {
      im_ok_tmx = 12 + im_tmx
      dum = sfwi[1:(dim(sfwi)[1] - 12), ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    } else {
      im_ok_tmx = im_tmx
      dum = sfwi[13:dim(sfwi)[1], ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    }
    im = best_m_spei[ireg] + firstmonth2[ireg,]
    if (im <= -12) {
      im_ok = 24 + im
      dum = spi[1:(dim(spi)[1] - 24), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else if (im > -12 & im <= 0) {
      im_ok = 12 + im
      dum = spi[13:(dim(spi)[1] - 12), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else {
      im_ok = im
      dum = spi[25:dim(spi)[1], ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    }
    
    mydata_clim = data.frame("y1" = spi_aux,
                             "y2" = sfwi_aux,
                             "x1" = (years))
    
    spi_aux = scale(mydata_clim$y1 - predict(lm(y1 ~ x1 , data = mydata_clim), newdata =
                                               mydata_clim))
    
    sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                mydata_clim))
    
    mydata_train_det = data.frame("y" = scale(mydata_ba_det$y),
                                  "x3" = sfwi_aux,
                                  "x4" = spi_aux)
    fit34 <- lm(y ~ x3 + x4, data = mydata_train_det)
    pre = fit34$coefficients[1] + fit34$coefficients[2] * mydata_train_det$x3 + fit34$coefficients[3] * mydata_train_det$x4
    sfwi_coef[ireg] = fit34$coefficients[2]
    spei_coef[ireg] = fit34$coefficients[3]
    aux = summary(fit34)
    VIF(fit34)
    cor.test(mydata_train_det$x3, mydata_train_det$x4)
    pcor(mydata_train_det)
    
    pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
    pval_spei_coef[ireg] = aux$coefficients[3, 4]
    
    plot.ts(mydata_train_det$y)
    lines(pre, col = "red")
    
    
  } else if (best_mod[ireg] == 2) {
    im_tmx = best_m_sfwi[ireg] + endmonth2[ireg,]
    if (im_tmx <= 0) {
      im_ok_tmx = 12 + im_tmx
      dum = spei_cc[1:(dim(sfwi)[1] - 12), ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    } else {
      im_ok_tmx = im_tmx
      dum = spei_cc[13:dim(sfwi)[1], ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    }
    im = best_m_spei[ireg] + firstmonth2[ireg,]
    if (im <= -12) {
      im_ok = 24 + im
      dum = spi[1:(dim(spi)[1] - 24), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else if (im > -12 & im <= 0) {
      im_ok = 12 + im
      dum = spi[13:(dim(spi)[1] - 12), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else {
      im_ok = im
      dum = spi[25:dim(spi)[1], ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    }
    
    mydata_clim = data.frame("y1" = spi_aux,
                             "y2" = sfwi_aux,
                             "x1" = (years))
    
    spi_aux = scale(mydata_clim$y1 - predict(lm(y1 ~ x1 , data = mydata_clim), newdata =
                                               mydata_clim))
    
    sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                mydata_clim))
    
    mydata_train_det = data.frame("y" = scale(mydata_ba_det$y),
                                  "x3" = sfwi_aux,
                                  "x4" = spi_aux)
    fit34 <- lm(y ~ x3 + x4, data = mydata_train_det)
    pre = fit34$coefficients[1] + fit34$coefficients[2] * mydata_train_det$x3 + fit34$coefficients[3] * mydata_train_det$x4
    sfwi_coef[ireg] = fit34$coefficients[2]
    spei_coef[ireg] = fit34$coefficients[3]
    aux = summary(fit34)
    
    pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
    pval_spei_coef[ireg] = aux$coefficients[3, 4]
    
    plot.ts(mydata_train_det$y)
    lines(pre, col = "red")
    
    plot.ts(mydata_train_det$x3)
    
    plot.ts(mydata_train_det$x4)
    
    pre = fit34$coefficients[1] + fit34$coefficients[2] * mydata_train_det$x3 + fit34$coefficients[3] * mydata_train_det$x4
    
    
    cor.test(mydata_train_det$x3, mydata_train_det$x4)
    
  } else if (best_mod[ireg] == 3) {
    im = best_m_spei[ireg] + firstmonth2[ireg,]
    if (im <= -12) {
      im_ok = 24 + im
      dum = spi[1:(dim(spi)[1] - 24), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else if (im > -12 & im <= 0) {
      im_ok = 12 + im
      dum = spi[13:(dim(spi)[1] - 12), ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    } else {
      im_ok = im
      dum = spi[25:dim(spi)[1], ireg, best_t_spei[ireg]]
      spi_aux = dum[seq(im_ok, length(dum), 12)]
    }
    
    mydata_clim = data.frame("y2" = spi_aux,
                             "x1" = (years))
    
    spi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                               mydata_clim))
    
    mydata_train_det = data.frame("y" = scale(mydata_ba_det$y),
                                  "x4" = spi_aux)
    fit4 <- lm(y ~ x4 , data = mydata_train_det)
    pre = fit4$coefficients[1] + fit4$coefficients[2] * mydata_train_det$x4
    aux = summary(fit4)
    spei_coef[ireg] = fit4$coefficients[2]
    pval_spei_coef[ireg] = aux$coefficients[2, 4]
    
    plot.ts(mydata_train_det$y)
    lines(pre, col = "red")
    
    
    
  } else if (best_mod[ireg] == 4) {
    im_tmx = best_m_sfwi[ireg] + endmonth2[ireg,]
    if (im_tmx <= 0) {
      im_ok_tmx = 12 + im_tmx
      dum = sfwi[1:(dim(sfwi)[1] - 12), ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    } else {
      im_ok_tmx = im_tmx
      dum = sfwi[13:dim(sfwi)[1], ireg,  best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    }
    
    mydata_clim = data.frame("y2" = sfwi_aux,
                             "x1" = (years))
    
    sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                mydata_clim))
    
    mydata_train_det = data.frame("y" = scale(mydata_ba_det$y),
                                  "x3" = sfwi_aux)
    fit3 <- lm(y ~ x3 , data = mydata_train_det)
    pre = fit3$coefficients[1] + fit3$coefficients[2] * mydata_train_det$x3
    aux = summary(fit3)
    sfwi_coef[ireg] = fit3$coefficients[2]
    pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
    
    plot.ts(mydata_train_det$y)
    lines(pre, col = "red")
    
  } else if (best_mod[ireg] == 5) {
    im_tmx = best_m_sfwi[ireg] + endmonth2[ireg,]
    if (im_tmx <= 0) {
      im_ok_tmx = 12 + im_tmx
      dum = spei_cc[1:(dim(sfwi)[1] - 12), ireg, best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    } else {
      im_ok_tmx = im_tmx
      dum = spei_cc[13:dim(sfwi)[1], ireg,  best_t_sfwi[ireg]]
      sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
    }
    
    mydata_clim = data.frame("y2" = sfwi_aux,
                             "x1" = (years))
    
    sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                mydata_clim))
    
    mydata_train_det = data.frame("y" = scale(mydata_ba_det$y),
                                  "x3" = sfwi_aux)
    fit3 <- lm(y ~ x3 , data = mydata_train_det)
    pre = fit3$coefficients[1] + fit3$coefficients[2] * mydata_train_det$x3
    aux = summary(fit3)
    sfwi_coef[ireg] = fit3$coefficients[2]
    pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
    
    plot.ts(mydata_train_det$y)
    lines(pre, col = "red")
    
    
  }
}


cc

plot.ts(mydata_train_det$y)
lines(pre, col = "red")



plot.ts(BA_esa_y[reg_test,])

plot.ts(log(BA_esa_y[reg_test,])+1)

BA=(log(BA_esa_y[reg_test,])+1)
mydata_ba = data.frame("y" = BA,
                       "x1" = (years))

mydata_ba_det = data.frame("y" = scale(mydata_ba$y - predict(lm(y ~ x1 , data = mydata_ba), newdata =
                                                               mydata_ba)))

plot.ts(mydata_ba_det)
plot.ts(BA_esa_y[reg_test,])

#  data
df <- data.frame(
  Year = years,
  BurnedArea = BA_esa_y[reg_test,],
  BurnedArea_std = mydata_ba_det,
  Model= pre
)

pdf(paste0(dir_out,"time_series_BA_",reg_test,".pdf"), width = 5, height = 3.5)  # Puoi personalizzare width e height come preferisci

# Create an empty plot with log scale
plot(df$Year, df$y, type = "n", 
     xlab = "Year", ylab = "Burned Area (log-transformed, detrended and standardized)", 
     ylim = c(min(df$y, na.rm = TRUE), max(df$y, na.rm = TRUE)),
     col = "black", cex.axis=1.2, cex.lab=1.2, cex.main=1.2, xaxt="n")

# Add custom x-axis to label the first year
axis(1, at = seq(2001, 2021, by = 5), labels = seq(2001, 2021, by = 5), cex.axis = 1.2)

# Add the time series for y
lines(df$Year, df$y, col = "black", lwd = 2)
points(df$Year, df$y, pch = 19, col = "black", cex=1.2)

# Add the time series for pre (Model) in red
lines(df$Year, df$Model, col = "red", lwd = 2)
points(df$Year, df$Model, pch = 19, col = "red", cex=1.2)

# Add legend
legend("topright", 
       legend = c("Observation", "Simulation"), 
       col = c("black", "red"), 
       lwd = 2, 
       pch = 19, 
       cex = 1.2, 
       box.lwd = 1.2)

dev.off()






########## annual cycle



## plot annual cycle
#a monthly climatology of burned area was produced from the entire time series.”
BAmm = array(NA, dim = c(12, nreg))
for (ireg in 1:nreg) {
  if (mask0[ireg] == 1) {
    BAm = array(NA, dim = c(12, dim(obs_reg)[1] / 12))
    for (im in (1:12)) {
      BAm[im, ] = obs_reg[seq(im, dim(obs_reg)[1], 12), ireg]
    }
    BAmm[,ireg] = apply(BAm, c(1), mean, na.rm = TRUE)
  }
}



barplot(BAmm[,reg_test])
plot.ts(BAmm[,reg_test])

burned_area <- BAmm[,reg_test]
endmonth2[reg_test]
firstmonth2[reg_test]
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# Create a vector of colors based on the month
colors <- ifelse(months %in% c( "Jul", "Aug", "Sep"), "red", "white")

pdf(paste0(dir_out,"annual_cycle_BA_",reg_test,".pdf"), width = 5, height = 3.5)  # Puoi personalizzare width e height come preferisci
bp <- barplot(burned_area, names.arg = months, col = colors, border = "black",
              main = "Monthly Burned Area", xlab = "Month", ylab = "Burned Area [ha]",
              las = 1, cex.names = 1.2,ylim = c(0, max(burned_area)*1.2))
# rect(xleft = (bp[firstmonth2[reg_test]-1] + bp[firstmonth2[reg_test]]) / 2, ybottom = 0, xright = (bp[endmonth2[reg_test]] + bp[endmonth2[reg_test]+1]) / 2, ytop = max(burned_area), col = "gray", border = NA)
barplot(burned_area, names.arg = months, col = colors, border = "black",
        main = "Monthly Burned Area", xlab = "Month", ylab = "Burned Area [ha]",
        las = 1, cex.names = 1.2, add = TRUE,ylim = c(0, max(burned_area)*1.2))
abline(v =(bp[firstmonth2[reg_test]-2] + bp[firstmonth2[reg_test]-1]) / 2, col = "black", lwd = 2)
abline(v = (bp[endmonth2[reg_test]] + bp[endmonth2[reg_test]+1])  / 2, col = "black", lwd = 2)
# Calcola la posizione verticale della freccia, mettiamola un po' sopra la barra più alta tra maggio e settembre
y_position <- max(burned_area[5:9]) * 1.05
# Linea con freccia tra maggio e settembre
arrows(x0 = (bp[firstmonth2[reg_test]-2] + bp[firstmonth2[reg_test]-1])  / 2, y0 = y_position, 
       x1 = (bp[endmonth2[reg_test]] + bp[endmonth2[reg_test]+1])  / 2, y1 = y_position, 
       angle = 30, code = 3, length = 0.1, lty = 2)
text(x = mean(c((bp[firstmonth2[reg_test]-2] + bp[firstmonth2[reg_test]-1])  / 2, (bp[endmonth2[reg_test]] + bp[endmonth2[reg_test]+1]) / 2)),  # Posizione media tra inizio e fine freccia
     y = y_position * 1.05,  # Un po' sopra la posizione della freccia
     labels = "CC", cex = 1.2)  # 'cex' controlla la dimensione del testo
arrows(x0 = bp[1], y0 = y_position, 
       x1 =  (bp[firstmonth2[reg_test]-2] + bp[firstmonth2[reg_test]-1])   / 2, y1 = y_position, 
       angle = 30, code = 2, length = 0.1, lty = 2) 
text(x = (bp[1] +  (bp[firstmonth2[reg_test]-2] + bp[firstmonth2[reg_test]-1])   / 2) / 2,  # Posizione media tra inizio e fine freccia
     y = y_position * 1.05,  # Un po' sopra la posizione della freccia
     labels = "AC", cex = 1.2)  # 'cex' controlla la dimensione del testo
dev.off()
