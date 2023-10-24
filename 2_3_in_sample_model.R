rm(list = ls())
graphics.off()
gc()

library(regclass)
library(lmtest)
library(car)


################## TO CHOOSE HERE
product='firecci51-nat' #'modis or 'firecci51' or or 'firecci51-nat'
### version
transf="log" #"log" "radq"
### local_environment
local_environment <- "onfire" #"onfire" "mac" "andrina"

if (local_environment=="onfire") {
  dir_out = '~/Dropbox/model/out_ecoregions/'
  dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
  dir_shp='/diskonfire/shapefiles/Ecoregions2017/'
} else if (local_environment=="mac") {
  dir_modis = '~/Dropbox/model/datos/'
  dir_out = '~/Dropbox/model/out_ecoregions/'
  dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
} else if (local_environment=="andrina") {
  dir_modis = '/home/andrina/Dropbox/model/datos/'
  dir_out = '/home/andrina/Dropbox/model/out_ecoregions/'
  dir_data = '/home/andrina/Desktop/climate-fire/datos/'
}


if (transf=="log") {
  if (product == "modis") {
    years = 2001:2021
    ## BA obs_reg
    load(paste0(dir_modis,
                paste('MCD64CMQ-2001-2021.Rdata')))
    # version = 'spei_spei_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018_fire_season_ja_2018'
    # version = 'spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018_fire_season_ja_2018'
  } else if (product == "firecci51") {
    years = 2001:2020
    load(paste0(
      dir_out,
      'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020.Rdata'
    ))
    version = 'spei_sfwi_firecci51_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018'
    # version = 'spei_spei_firecci51_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018'
  } else if (product == "firecci51-nat") {
    years = 2001:2020
    load(
      paste0(
        dir_out,
        'ESACCI-L4_FIRE-BA-MODIS-fv5.1-2001-2020-natural-fires.Rdata'
      )
    )
    # version = 'spei_sfwi_firecci51-nat_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018'
    version = 'spei_spei_firecci51-nat_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018'
  }
} else if (transf == "radq") {
  years = 2001:2021
  # version = 'spei_sfwi_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018'
  # version = 'spei_spei_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018'
}


# Function definition
check_model_residuals <- function(model) {
  resid_values <- as.numeric(model$residuals)
  ks_test <- ks.test(resid_values, "pnorm", mean(resid_values), sd(resid_values))$p.value
  # 
  # Shapiro-Wilk normality test
  # shapiro_test <- shapiro.test(residuals(model))$p.value
  
  # Breusch-Pagan test for heteroscedasticity
  bp_test <- bptest(model)$p.value
  
  # Durbin-Watson test for autocorrelation
  dw_test <- durbinWatsonTest(model)$p
  
  # Return 1 if all tests pass, 0 otherwise
  if (ks_test > 0.05 & bp_test > 0.05 & dw_test > 0.05) {
    return(1)
  } else {
    return(0)
  }
}


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
  if (length(which(BA_esa_y[ireg, ] > 0)) >= 10 &&
      100 * auxsum[ireg] / sum(auxsum) > 0.001) {
    mask[ireg] = 1
  }
}


if (transf=="log") {
  BA = log(BA_esa_y+1)
} else if (transf=="radq") {
  BA = sqrt(BA_esa_y)
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


### 
corr = array(NA, dim = c(nreg))
sig = array(NA, dim = c(nreg))
sfwi_coef = array(NA, dim = c(nreg))
spei_coef = array(NA, dim = c(nreg))
pval_sfwi_coef = array(NA, dim = c(nreg))
pval_spei_coef = array(NA, dim = c(nreg))


load(paste0(dir_out, "sig_", version, ".RData"))

load(paste0(dir_out, "best_m_SPI_fin_", version, ".RData")) #c(-14,-2)
load(paste0(dir_out, "best_t_SPI_fin_", version, ".RData"))
load(paste0(dir_out, "best_m_sfwi_fin_", version, ".RData")) #c(-12, 0)
load(paste0(dir_out, "best_t_sfwi_fin_", version, ".RData"))


pval_sfwi_coef_num = 0
pval_spei_coef_num = 0
pval_sfwi_coef_num_2 = 0
pval_spei_coef_num_2 = 0


kk = 0
  for (ireg in 1:nreg) {
  print(paste0('reg ', ireg, '/', nreg))
  
  if (!is.na(best_sig[ireg])) {
    mydata_ba = data.frame("y" = BA[ireg,],
                           "x1" = (years))
    
    mydata_ba_det = data.frame("y" = scale(mydata_ba$y - predict(lm(y ~ x1 , data = mydata_ba), newdata =
                                                                   mydata_ba)))
    
    
    if (!is.na(best_m_sfwi_fin[ireg]) &
        is.na(best_m_SPI_fin[ireg])) {
      #c(-12, 0)
      # for (im_tmx in ((first_month[i,j] - 1):last_month[i,j])) {
      # best_m_sfwi[k] = im_tmx-last_month[i,j]
      # --> im_tmx=best_m_sfwi[k]+last_month[i,j]
      im_tmx = best_m_sfwi_fin[ireg] + endmonth[ireg]
      if (im_tmx <= 0) {
        im_ok_tmx = 12 + im_tmx
        dum = sfwi[1:(dim(sfwi)[1] - 12),ireg,best_t_sfwi_fin[ireg]]
        sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
      } else {
        im_ok_tmx = im_tmx
        dum = sfwi[13:dim(sfwi)[1],ireg,  best_t_sfwi_fin[ireg]]
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
      
      # resid_values <- as.numeric(fit3$residuals)
      # shapiro.test(resid_values)
      # hist(resid_values)
      # shapiro.test(residuals(fit3))
      # 
      # ks.test(resid_values, "pnorm", mean(resid_values), sd(resid_values))
      # 
      # 
      # bptest(fit3)$p.value
      # 
      # # Durbin-Watson test for autocorrelation
      # dw_test <- durbinWatsonTest(fit3)$p
      # bgtest(fit3, order = 1)
      
      result <- check_model_residuals(fit3)
      
      
    }
    
    if (!is.na(best_m_SPI_fin[ireg]) &
        is.na(best_m_sfwi_fin[ireg])) {
      im = best_m_SPI_fin[ireg] + firstmonth[ireg]
      
      if (im <= -12) {
        im_ok = 24 + im
        dum = spi[ 1:(dim(spi)[1] - 24),ireg, best_t_SPI_fin[ireg]]
        spi_aux = dum[seq(im_ok, length(dum), 12)]
      } else if (im > -12 & im <= 0) {
        im_ok = 12 + im
        dum = spi[13:(dim(spi)[1] - 12), ireg, best_t_SPI_fin[ireg]]
        spi_aux = dum[seq(im_ok, length(dum), 12)]
      } else {
        im_ok = im
        dum = spi[25:dim(spi)[1],ireg, best_t_SPI_fin[ireg]]
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
      
      result <- check_model_residuals(fit4)
      
    }
    
    if (!is.na(best_m_SPI_fin[ireg]) &
        !is.na(best_m_sfwi_fin[ireg])) {
      im_tmx = best_m_sfwi_fin[ireg] + endmonth[ireg]
      if (im_tmx <= 0) {
        im_ok_tmx = 12 + im_tmx
        dum = sfwi[1:(dim(sfwi)[1] - 12), ireg, best_t_sfwi_fin[ireg]]
        sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
      } else {
        im_ok_tmx = im_tmx
        dum = sfwi[13:dim(sfwi)[1],ireg, best_t_sfwi_fin[ireg]]
        sfwi_aux = dum[seq(im_ok_tmx, length(dum), 12)]
      }
      im = best_m_SPI_fin[ireg] + firstmonth[ireg]
      if (im <= -12) {
        im_ok = 24 + im
        dum = spi[ 1:(dim(spi)[1] - 24), ireg,best_t_SPI_fin[ireg]]
        spi_aux = dum[seq(im_ok, length(dum), 12)]
      } else if (im > -12 & im <= 0) {
        im_ok = 12 + im
        dum = spi[13:(dim(spi)[1] - 12), ireg,best_t_SPI_fin[ireg]]
        spi_aux = dum[seq(im_ok, length(dum), 12)]
      } else {
        im_ok = im
        dum = spi[ 25:dim(spi)[1], ireg, best_t_SPI_fin[ireg]]
        spi_aux = dum[seq(im_ok, length(dum), 12)]
      }
      
      
      mydata_clim = data.frame("y1" = spi_aux,
                               "y2" = sfwi_aux,
                               "x1" = (years))
      
      spi_aux = scale(mydata_clim$y1 - predict(lm(y1 ~ x1 , data = mydata_clim), newdata =
                                                 mydata_clim))
      
      sfwi_aux = scale(mydata_clim$y2 - predict(lm(y2 ~ x1 , data = mydata_clim), newdata =
                                                  mydata_clim))
      
      mydata_train_det = data.frame(
        "y" = scale(mydata_ba_det$y),
        "x3" = sfwi_aux,
        "x4" = spi_aux
      )
      fit34 <- lm(y ~ x3 + x4, data = mydata_train_det)
      
      result <- check_model_residuals(fit34)
      if (result ==1) {
        
      
      pre = fit34$coefficients[1] + fit34$coefficients[2] * mydata_train_det$x3 + fit34$coefficients[3] * mydata_train_det$x4
      sfwi_coef[ireg] = fit34$coefficients[2]
      spei_coef[ireg] = fit34$coefficients[3]
      aux = summary(fit34)
      pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
      pval_spei_coef[ireg] = aux$coefficients[3, 4]
      
      if (pval_sfwi_coef[ireg] > 0.05 &
          pval_spei_coef[ireg] > 0.05) {
        sfwi_coef[ireg] = NA
        spei_coef[ireg] = NA
      } else if (pval_sfwi_coef[ireg] > 0.05 &
                 pval_spei_coef[ireg] <= 0.05) {
        sfwi_coef[ireg] = NA
        fit34_1 <- lm(y ~  x4, data = mydata_train_det)
        pre = fit34_1$coefficients[1]  + fit34_1$coefficients[2] * mydata_train_det$x4
        spei_coef[ireg] = fit34_1$coefficients[2]
        aux = summary(fit34_1)
        pval_spei_coef[ireg] = aux$coefficients[2, 4]
        if (pval_spei_coef[ireg] > 0.05) {
          sfwi_coef[ireg] = NA
          spei_coef[ireg] = NA
        }
      } else if (pval_spei_coef[ireg] > 0.05 &
                 pval_sfwi_coef[ireg] <= 0.05) {
        pval_spei_coef_num_2 = pval_spei_coef_num_2 + 1
        fit34_2 <- lm(y ~  x3, data = mydata_train_det)
        spei_coef[ireg] = NA
        pre = fit34_2$coefficients[1]  + fit34_2$coefficients[2] * mydata_train_det$x3
        sfwi_coef[ireg] = fit34_2$coefficients[2]
        aux = summary(fit34_2)
        pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
        if (pval_sfwi_coef[ireg] > 0.05) {
          sfwi_coef[ireg] = NA
          spei_coef[ireg] = NA
        }
      }
      
      # clean_data <- na.omit(mydata_train_det)
      # 
      # par_cor=pcor(clean_data)
      # if (par_cor$p.value[1,3] > 0.05) {
      #   spei_coef[ireg] = NA
      #   fit34_2 <- lm(y ~  x3, data = mydata_train_det)
      #   pre = fit34_2$coefficients[1]  + fit34_2$coefficients[2] * mydata_train_det$x3
      #   sfwi_coef[ireg] = fit34_2$coefficients[2]
      #   aux = summary(fit34_2)
      #   pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
      #   if (pval_sfwi_coef[ireg] > 0.05) {
      #     sfwi_coef[ireg] = NA
      #     spei_coef[ireg] = NA
      #   }
      # }

      vif34=VIF(fit34)
      vif_thres=2
      if (vif34[1] > vif_thres &
          vif34[2] > vif_thres) {
        spei_coef[ireg] = NA
        fit34_2 <- lm(y ~  x3, data = mydata_train_det)
        pre = fit34_2$coefficients[1]  + fit34_2$coefficients[2] * mydata_train_det$x3
        sfwi_coef[ireg] = fit34_2$coefficients[2]
        aux = summary(fit34_2)
        pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
        if (pval_sfwi_coef[ireg] > 0.05) {
          sfwi_coef[ireg] = NA
          spei_coef[ireg] = NA
        }
      } else if (vif34[1] > vif_thres &
                 vif34[2] <= vif_thres) {
        sfwi_coef[ireg] = NA
        pval_sfwi_coef_num_2 = pval_sfwi_coef_num_2 + 1
        fit34_1 <- lm(y ~  x4, data = mydata_train_det)
        pre = fit34_1$coefficients[1]  + fit34_1$coefficients[2] * mydata_train_det$x4
        spei_coef[ireg] = fit34_1$coefficients[2]
        aux = summary(fit34_1)
        pval_spei_coef[ireg] = aux$coefficients[2, 4]
        if (pval_spei_coef[ireg] > 0.05) {
          sfwi_coef[ireg] = NA
          spei_coef[ireg] = NA
        }
      } else if (vif34[1] <= vif_thres &
                 vif34[2] > vif_thres) {
        pval_spei_coef_num_2 = pval_spei_coef_num_2 + 1
        fit34_2 <- lm(y ~  x3, data = mydata_train_det)
        spei_coef[ireg] = NA
        pre = fit34_2$coefficients[1]  + fit34_2$coefficients[2] * mydata_train_det$x3
        sfwi_coef[ireg] = fit34_2$coefficients[2]
        aux = summary(fit34_2)
        pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
        if (pval_sfwi_coef[ireg] > 0.05) {
          sfwi_coef[ireg] = NA
          spei_coef[ireg] = NA
        }
      }
      
      } else {
        
        fit34_2 <- lm(y ~  x3, data = mydata_train_det)
        
        pre = fit34_2$coefficients[1]  + fit34_2$coefficients[2] * mydata_train_det$x3
        sfwi_coef[ireg] = fit34_2$coefficients[2]
        spei_coef[ireg] = NA
        aux = summary(fit34_2)
        pval_sfwi_coef[ireg] = aux$coefficients[2, 4]
        
        if (pval_sfwi_coef[ireg] > 0.05 | check_model_residuals(fit34_2)==0) {
          sfwi_coef[ireg] = NA
          
          fit34_3 <- lm(y ~  x4, data = mydata_train_det)
          aux_fit34_3 = summary(fit34_3)
          pval_spei_coef[ireg] = aux_fit34_3$coefficients[2, 4]
          if (pval_spei_coef[ireg] < 0.05 & check_model_residuals(fit34_3)==1) {
            pre = fit34_3$coefficients[1]  + fit34_3$coefficients[2] * mydata_train_det$x4
            spei_coef[ireg] = fit34_3$coefficients[2]
            if(spei_coef[ireg] <0) {oooo}
          } 
            
        }
        
        
        
      }
    }
    
    
    if (length(which(!is.na(c(
      sfwi_coef[ireg], spei_coef[ireg]
    )))) >= 1) {
      kk = kk + 1
      rho = cor.test(mydata_ba_det$y,
                     pre,
                     use = "pairwise.complete.obs",
                     alternative = "greater")
      
      
      corr[ireg] = rho$estimate
      sig[ireg] = rho$p.value
      
    }
    
  }
}

# cc2


length(which(!is.na(corr))) / length(which(mask==1))
# length(which(!is.na(sfwi_coef))) / length(sfwi_coef)
length(which(!is.na(spei_coef))) 

summary(as.vector(corr))
summary(as.vector(corr * corr))

sfwi_coef[19]
spei_coef[19]

save(corr,
     file = paste0(dir_out, "corr_reconstruction_", version, ".RData"))
save(sfwi_coef,
     file = paste0(dir_out, "sfwi_coef_reconstruction_", version, ".RData"))
save(spei_coef,
     file = paste0(dir_out, "spei_coef_reconstruction_", version, ".RData"))
save(
  pval_sfwi_coef,
  file = paste0(dir_out, "pval_sfwi_coef_reconstruction_", version, ".RData")
)
save(
  pval_spei_coef,
  file = paste0(dir_out, "pval_spei_coef_reconstruction_", version, ".RData")
)

