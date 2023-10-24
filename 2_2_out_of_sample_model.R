#===============================================================================
# Preliminary Setup
#===============================================================================

# Clearing the workspace and graphic device
rm(list = ls())
graphics.off()
gc()


################## TO CHOOSE HERE
product='firecci51-nat' #'modis or 'firecci51' or or 'firecci51-nat'

################## TO SELECT THE VERSION BELOW


# Directories
# Setting directories
# Adjust `local_environment` flag according to the environment you're working on
local_environment <- TRUE 

if (local_environment) {
  dir_data = '/Users/marco/Documents/dati/fire_climate_data/climate_fire/datos/'
  dir_fwi = '/Users/marco/Documents/dati/obs/ERA5/FWI/'
} else {
  dir_data = '/diskonfire/CLIMATE_FIRE/datos/'
  dir_fwi = '/diskonfire/ERA5/FWI/'
}
dir_modis = '~/Dropbox/model/datos/'
dir_out = '~/Dropbox/model/out_ecoregions/'



if (product=='firecci51') {
  years = 2001:2020
  load(paste0(dir_out, 'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020.Rdata'))
  # version = 'spei_sfwi_firecci51_mswep_era5_fireseason_ecoregions_log'
  # version = 'spei_spei_firecci51_mswep_era5_fireseason_ecoregions_log'
} else if (product == 'firecci51-nat') {
  years = 2001:2020
  load(paste0(dir_out, paste('ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020-natural-fires.Rdata')))
  version = 'spei_sfwi_firecci51-nat_mswep_era5_fireseason_ecoregions_log'
  # version = 'spei_spei_firecci51-nat_mswep_era5_fireseason_ecoregions_log'
} else if (product == 'modis') {
  years = 2001:2021
  load(paste0(dir_modis, paste('MCD64CMQ-2001-2021.Rdata')))
  # version = 'spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log'
  # version = 'spei_spei_modis_mswep_era5_fireseason_ecoregions_log'
  
  # version = 'spei_sfwi_modis_mswep_era5_fireseason_ecoregions'
  # version = 'spei_spei_modis_mswep_era5_fireseason_ecoregions'
}

#===============================================================================
# Data Loading and Processing
#===============================================================================

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
  if (length(which(BAy[ireg, ] > 0)) >= 2) {
    mask0[ireg] = 1
  }
}

load(paste0(dir_out, 'firstmonth_reg_',product,'.RData'))
load(paste0(dir_out, 'endmonth_reg_',product,'.RData'))

## aggregate BA over the fire season
BA_esa_y = array(NA, dim = c(nreg, length(years)))
for (ireg in 1:nreg) {
  if (mask0[ireg] == 1) {
    for (iyear in 1:length(years)) {
      if (firstmonth[ireg] >= 1) {
        i1 = (iyear - 1) * 12 + firstmonth[ireg]
        i2 = (iyear - 1) * 12 + endmonth[ireg]
        BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
      } else
        if (iyear == 1) {
          next
        } else {
          i1 = (iyear - 1) * 12 + firstmonth[ireg]
          i2 = (iyear - 1) * 12 + endmonth[ireg]
          BA_esa_y[ireg, iyear] = sum(obs_reg[i1:i2, ireg], na.rm = TRUE)
        }
    }
  }
}

auxsum = apply(BA_esa_y, c(1), sum, na.rm = T)
mask = array(0, dim = c(nreg))
for (ireg in 1:nreg) {
  if (length(which(BA_esa_y[ireg,] > 0)) >= 10 &&
      100 * auxsum[ireg] / sum(auxsum) > 0.001) {
    mask[ireg] = 1
  }
}



if (grepl("^spei_sfwi_", version))  {
  # fwi_reg
  load(paste0(dir_out, "sfwi3-2000-2021-ecoregions.RData"))
  load(paste0(dir_out, "sfwi6-2000-2021-ecoregions.RData"))
  load(paste0(dir_out, "sfwi12-2000-2021-ecoregions.RData"))
  
  #spei_reg
  load(paste0(dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei6-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei12-1999-2021-ecoregions-mswep-era5.RData"))
  
} else if (grepl("^spei_spei_", version))  {
  load(paste0(dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei6-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei12-1999-2021-ecoregions-mswep-era5.RData"))
  years_sfwi = 1999:2021 #this is different from the spei_sfwi model
  iok_sfwi = match((years[1] - 1):years[length(years)], years_sfwi) #no need 1999
  #sfwi coef is suppose to be >0, while the spi effect is supposet to be <0. to have the same considtions below, i changes the sign of the spi
  
  fwi3 = -spei3[(((iok_sfwi[1] - 1) * 12) + 1):nrow(spei3), ]
  fwi6 = -spei6[(((iok_sfwi[1] - 1) * 12) + 1):nrow(spei3), ]
  fwi12 = -spei12[(((iok_sfwi[1] - 1) * 12) + 1):nrow(spei3), ]
  
  ## load SPEI
  load(paste0(dir_out, "spei3-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei6-1999-2021-ecoregions-mswep-era5.RData"))
  load(paste0(dir_out, "spei12-1999-2021-ecoregions-mswep-era5.RData"))
}

sfwi = array(NA, dim = c(dim(fwi3)[1], nreg, 3))
sfwi[, , 1] = fwi3
sfwi[, , 2] = fwi6
sfwi[, , 3] = fwi12
rm(fwi3)
rm(fwi6)
rm(fwi12)


spi = array(NA, dim = c(dim(spei3)[1], nreg, 3))
spi[, ,  1] = spei3
spi[, ,  2] = spei6
spi[, ,  3] = spei12
rm(spei3)
rm(spei6)
rm(spei12)

# delete 2021 if firecci51 based on version
if (grepl("firecci51", product)) {
  spi <- spi[1:(nrow(spi)-12), , ]
  sfwi <- sfwi[1:(nrow(sfwi)-12), , ]
} 

# Transformation based on version
if (grepl("log", version)) {
  BA <- log(BA_esa_y + 1)
} else {
  BA <- sqrt(BA_esa_y)
}


## find the best predictor
best_corr = array(NA, dim = c(nreg))
# num_iter = 324
# num_iter = 1350
num_iter = 1755
best_sig = array(NA, dim = c(nreg))
best_m_SPI_fin = array(NA, dim = c(nreg))
best_t_SPI_fin = array(NA, dim = c(nreg))
best_m_sfwi_fin = array(NA, dim = c(nreg))
best_t_sfwi_fin = array(NA, dim = c(nreg))
best_x3_fin = array(NA, dim = c(nreg, length(years)))
best_x4_fin = array(NA, dim = c(nreg, length(years)))

for (ireg in 1:nreg) {
  # for (i in 121) {
  
  if (mask[ireg] == 1) {
    mydata_ba = data.frame("y" = BA[ireg,],
                           "x1" = (years))
    
    mydata_ba_det = data.frame("y" = scale(mydata_ba$y - predict(lm(y ~ x1 , data = mydata_ba), newdata =
                                                                   mydata_ba)))
    
    corr = array(NA, dim = c(num_iter))
    sig = array(NA, dim = c(num_iter))
    best_m_SPI = array(NA, dim = c(num_iter))
    best_t_SPI = array(NA, dim = c(num_iter))
    best_m_sfwi = array(NA, dim = c(num_iter))
    best_t_sfwi = array(NA, dim = c(num_iter))
    
    best_x3_3 = array(NA, dim = c(num_iter, length(years)))
    best_x3_34 = array(NA, dim = c(num_iter, length(years)))
    
    best_x4_4 = array(NA, dim = c(num_iter, length(years)))
    best_x4_34 = array(NA, dim = c(num_iter, length(years)))
    
    best_x3 = array(NA, dim = c(num_iter, length(years)))
    best_x4 = array(NA, dim = c(num_iter, length(years)))
    
    k = 0
    for (im in (firstmonth[ireg] - 14):(firstmonth[ireg] - 2)) {
      # for (im in ((1 - 14):1-2)) {
      
      for (isc in 1:3) {
        # for (im_tmx in ((1 - 1):12)) {
        for (im_tmx in ((firstmonth[ireg] - 1):endmonth[ireg])) {
          # for (im_tmx in ((12 - 5):12)) {
          for (isc_tmx in 1:3) {
            k = k + 1
            print(paste0('reg ',
                         ireg,
                         '/',
                         nreg,
                         '; step ',
                         k,
                         '/',
                         num_iter))
            
            spi_aux  = vector()
            sfwi_aux = vector()
            if (im <= -12) {
              im_ok = 24 + im
              dum = spi[1:(dim(spi)[1] - 24), ireg, isc]
              spi_aux = dum[seq(im_ok, length(dum), 12)]
            } else if (im > -12 & im <= 0) {
              im_ok = 12 + im
              dum = spi[13:(dim(spi)[1] - 12), ireg, isc]
              spi_aux = dum[seq(im_ok, length(dum), 12)]
            } else {
              im_ok = im
              dum = spi[25:dim(spi)[1], ireg, isc]
              spi_aux = dum[seq(im_ok, length(dum), 12)]
            }
            
            
            if (im_tmx <= 0) {
              im_ok_tmx = 12 + im_tmx
              dum = sfwi[1:(dim(sfwi)[1] - 12), ireg, isc_tmx]
              sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
            } else {
              im_ok_tmx = im_tmx
              dum = sfwi[13:dim(sfwi)[1], ireg, isc_tmx]
              sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
            }
            
            if (length(which(is.na(spi_aux))) >= round((length(spi_aux) -
                                                        1) * 0.9)) {
              next
            }
            
            if (length(which(is.na(spi_aux))) == length(spi_aux)) {
              next
            }
            
            if (length(which(is.na(sfwi_aux))) == length(sfwi_aux)) {
              next
            }
            
            
            mydata_clim = data.frame("y1" = spi_aux,
                                     "y2" = sfwi_aux,
                                     "x1" = (years))
            
            spi_aux = scale(mydata_clim$y1 - predict(lm(y1 ~ x1 , data = mydata_clim), newdata =
                                                       mydata_clim))
            
            sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                        mydata_clim))
            
            # spi_aux=scale(spi_aux)
            # sfwi_aux=scale(sfwi_aux)
            
            pre3 = vector()
            pre4 = vector()
            pre34 = vector()
            
            for (iy in 1:length(years)) {
              if (iy == 1 && firstmonth[ireg] <= 0) {
                next
              }
              BA_train = BA[ireg,-iy]
              spi_train = spi_aux[-iy]
              sfwi_train = sfwi_aux[-iy]
              years_train = years[-iy]
              
              spi_test = spi_aux[iy]
              sfwi_test = sfwi_aux[iy]
              years_test = years[iy]
              
              mydata_train = data.frame("y" = (BA_train),
                                        "x1" = (years_train))
              
              mydata_train_det = data.frame(
                "y" =
                  scale(
                    mydata_train$y - predict(lm(y ~ x1 , data = mydata_train), newdata =
                                               mydata_train)
                  ),
                # "x1" = (years_train),
                "x3" = sfwi_train,
                "x4" = spi_train
              )
              
              mydata_test = data.frame("x3" = (sfwi_test),
                                       "x4" = (spi_test))
              
              fit3 <- lm(y ~ x3 , data = mydata_train_det)
              fit4 <- lm(y ~ x4 , data = mydata_train_det)
              fit34 <- lm(y ~ x3 + x4, data = mydata_train_det)
              
              # pre3[iy] = fit3$coefficients[1] +fit3$coefficients[2]*mydata_test$x1 + fit3$coefficients[3] * mydata_test$x3
              # pre4[iy] = fit4$coefficients[1] +fit4$coefficients[2]*mydata_test$x1 + fit4$coefficients[3] * mydata_test$x4
              # pre34[iy] = fit34$coefficients[1] +fit34$coefficients[2]*mydata_test$x1 + fit34$coefficients[3] * mydata_test$x3 + fit34$coefficients[4] * mydata_test$x4
              
              pre3[iy] = fit3$coefficients[1] + fit3$coefficients[2] * mydata_test$x3
              pre4[iy] = fit4$coefficients[1] + fit4$coefficients[2] * mydata_test$x4
              pre34[iy] = fit34$coefficients[1] + fit34$coefficients[2] * mydata_test$x3 + fit34$coefficients[3] * mydata_test$x4
              
              
              # best_x3_3[k, iy] = fit3$coefficients[3]
              # best_x3_34[k, iy] = fit34$coefficients[3]
              #
              # best_x4_4[k, iy] = fit4$coefficients[3]
              # best_x4_34[k, iy] = fit34$coefficients[4]
              best_x3_3[k, iy] = fit3$coefficients[2]
              best_x3_34[k, iy] = fit34$coefficients[2]
              
              best_x4_4[k, iy] = fit4$coefficients[2]
              best_x4_34[k, iy] = fit34$coefficients[3]
              
            }
            
            
            # rho3 = cor.test(scale(BA[i,j,]),
            rho3 = cor.test(mydata_ba_det$y,
                            pre3,
                            use = "pairwise.complete.obs",
                            alternative = "greater")
            rho4 = cor.test(mydata_ba_det$y,
                            pre4,
                            use = "pairwise.complete.obs",
                            alternative = "greater")
            rho34 = cor.test(mydata_ba_det$y,
                             pre34,
                             use = "pairwise.complete.obs",
                             alternative = "greater")
            
            ## only x3 and x4 positive
            if (rho3$estimate > rho4$estimate &
                rho3$estimate > rho34$estimate &
                min(best_x3_3[k, ], na.rm = TRUE) > 0) {
              corr[k] = rho3$estimate
              sig[k] = rho3$p.value
              best_x3[k, ] = best_x3_3[k, ]
              best_t_sfwi[k] = isc_tmx
              best_m_sfwi[k] = im_tmx - endmonth[ireg]
            }
            if (rho4$estimate > rho3$estimate &
                rho4$estimate > rho34$estimate &
                min(best_x4_4[k, ], na.rm = TRUE) > 0) {
              corr[k] = rho4$estimate
              sig[k] = rho4$p.value
              best_x4[k, ] = best_x4_4[k, ]
              best_t_SPI[k] = isc
              best_m_SPI[k] = im - (firstmonth[ireg])
            }
            if (rho34$estimate > rho3$estimate &
                rho34$estimate > rho4$estimate &
                min(best_x3_34[k, ], na.rm = TRUE) > 0 &
                min(best_x4_34[k, ], na.rm = TRUE) > 0) {
              corr[k] = rho34$estimate
              sig[k] = rho34$p.value
              best_x3[k, ] = best_x3_34[k, ]
              best_x4[k, ] = best_x4_34[k, ]
              best_t_sfwi[k] = isc_tmx
              best_m_sfwi[k] = im_tmx - endmonth[ireg]
              best_t_SPI[k] = isc
              best_m_SPI[k] = im - (firstmonth[ireg])
            }
          }
        }
      }
    }
    
    
    sig = p.adjust(sig, method = "fdr")
    corr[sig > 0.05] = NA
    
    if (length(which(is.na(corr))) == length(corr)) {
      next
    }
    
    dum = max((corr), na.rm = TRUE)
    idx = which(abs(corr) == dum, arr.ind = TRUE)
    idx = idx[1]
    best_corr[ireg] = corr[idx]
    best_sig[ireg] = sig[idx]
    best_m_SPI_fin[ireg] = best_m_SPI[idx]
    best_t_SPI_fin[ireg] = best_t_SPI[idx]
    best_m_sfwi_fin[ireg] = best_m_sfwi[idx]
    best_t_sfwi_fin[ireg] = best_t_sfwi[idx]
    #best_sig = sig[idx]
    
    
    best_x3_fin[ireg,] = best_x3[idx,]
    best_x4_fin[ireg,] = best_x4[idx,]
    
    
  }
}



length(which(is.na(best_corr))) / length(best_corr)



save(best_corr,
     file = paste0(dir_out, "corr_", version, "_fire_season_ja_2018.RData"))
save(best_sig,
     file = paste0(dir_out, "sig_", version, "_fire_season_ja_2018.RData"))
save(
  best_m_SPI_fin,
  file = paste0(
    dir_out,
    "best_m_SPI_fin_",
    version,
    "_fire_season_ja_2018.RData"
  )
)
save(
  best_t_SPI_fin,
  file = paste0(
    dir_out,
    "best_t_SPI_fin_",
    version,
    "_fire_season_ja_2018.RData"
  )
)
save(
  best_m_sfwi_fin,
  file = paste0(
    dir_out,
    "best_m_sfwi_fin_",
    version,
    "_fire_season_ja_2018.RData"
  )
)
save(
  best_t_sfwi_fin,
  file = paste0(
    dir_out,
    "best_t_sfwi_fin_",
    version,
    "_fire_season_ja_2018.RData"
  )
)

save(best_x3_fin,
     file = paste0(dir_out,
                   "best_x3_",
                   version,
                   "_fire_season_ja_2018.RData"))
save(best_x4_fin,
     file = paste0(dir_out,
                   "best_x4_",
                   version,
                   "_fire_season_ja_2018.RData"))

