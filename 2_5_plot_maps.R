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

################## TO CHOOSE HERE
product='modis' #'modis or 'firecci51' or or 'firecci51-nat'
### version
transf="log" #"log" "radq"
### local_environment
local_environment <- "onfire" #"onfire" "mac" "andrina"

if (local_environment=="onfire") {
  dir_out = '~/Dropbox/model/out_ecoregions/'
  dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
  dir_shp='/diskonfire/shapefiles/Ecoregions2017/'
} else if (local_environment=="mac") {
  dir_out = '~/Dropbox/model/out_ecoregions/'
  dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
  dir_shp='/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'
} else if (local_environment=="andrina") {
  dir_out = '/home/andrina/Dropbox/model/out_ecoregions/'
  dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
  dir_shp='/home/andrina/Desktop/climate-fire/datos/shp/'
}


if (transf=="log") {
  if (product=="modis") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018') 
  } else if (product=="firecci51") {
    list_version=c('spei_sfwi_firecci51_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018','spei_spei_firecci51_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018') 
  } else if (product=="firecci51-nat") {
    list_version=c('spei_sfwi_firecci51-nat_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018','spei_spei_firecci51-nat_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018') 
  }
} else if (transf=="radq") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018')
}

list_version_names=c('spei-sfwi','spei-spei')
list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')

load(paste0(dir_out,
            "mask_ecoregions_season_ja_2018_",product,".RData"))

nreg=dim(mask)

all_r2 <- array(NA, dim = c(nreg,length(list_version)))
all_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_spei <- array(NA, dim = c(nreg,length(list_version)))
all_m_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_m_spei <- array(NA, dim = c(nreg,length(list_version)))
all_t_sfwi <- array(NA, dim = c(nreg,length(list_version)))
all_t_spei <- array(NA, dim = c(nreg,length(list_version)))

for (k in 1:length(list_version)) {
  version=list_version[k]
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
 
  all_r2[,k]=corr*corr
  all_spei[,k]=spei_coef
  all_sfwi[,k]=sfwi_coef
  all_t_spei[,k]=best_t_SPI_fin
  all_t_sfwi[,k]=best_t_sfwi_fin
  all_m_spei[,k]=best_m_SPI_fin
  all_m_sfwi[,k]=best_m_sfwi_fin
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
      best_sfwi_coef[ireg] <- all_sfwi[ireg,aux]
    }
    if (!is.na(all_spei[ireg, aux]) & is.na(all_sfwi[ireg, aux])) {
      aux2 = substr(list_version_names[aux], 1, as.numeric(tmp[[1]]) - 1)
      aux2 = paste0('^ANT', aux2, '$')
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_spei_coef[ireg] <- all_spei[ireg,aux]
    }
    if (!is.na(all_spei[ireg, aux]) &
        !is.na(all_sfwi[ireg, aux])) {
      aux2 = list_version_names[aux]
      best_m_sfwi[ireg] <- all_m_sfwi[ireg, aux]
      best_t_sfwi[ireg] <- all_t_sfwi[ireg, aux]
      best_m_spei[ireg] <- all_m_spei[ireg, aux]
      best_t_spei[ireg] <- all_t_spei[ireg, aux]
      best_sfwi_coef[ireg] <- all_sfwi[ireg,aux]
      best_spei_coef[ireg] <- all_spei[ireg,aux]
    }
    
    best_mod[ireg] = grep(aux2, list_models)
    if (best_mod[ireg] == 2 ||
        best_mod[ireg] == 5) {
      best_sfwi_coef[ireg] = -best_sfwi_coef[ireg]
    }
  }
}

summary(as.vector(best_r2))

print(paste0(
  "Variance explained (spatial median) = ",
  round(median(best_r2, na.rm = TRUE), digits = 2),
  "\n",
  "Percentage valid models = ",
  round(100 * length(which(!is.na(
    best_r2
  ))) / length(which(mask == 1)), digits = 0),
  '%'))




########################PLOT 
file_shp = file.path(dir_shp, "Ecoregions2017_repaired.shp")
eco <- st_read(file_shp) %>% st_make_valid()

area_eco <- st_area(eco)
sum(area_eco[!is.na(best_r2)])/sum(area_eco[mask == 1])

wld <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") #https://github.com/ropensci/rnaturalearth
# Exclude Antarctica 
wld <- dplyr::filter(wld, name_long != "Antarctica")



# unir a los datos vectoriales
eco <- mutate(eco, best_r2)
#
eco_only <- filter(eco, !is.na(best_r2))

colors <- rev(inferno(5))
colors <- c( "darkgray", colors)

# Modificaciones en la parte de ggplot
f1 <- ggplot(eco_only) +
  geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
  geom_sf(aes(fill = cut(best_r2, breaks = c(0, 0.25, 0.50 , 0.75, 1), include.lowest = TRUE)), 
          colour = "white", linewidth = 0.01) +
  geom_sf(data = subset(eco, is.na(best_r2) & mask == 1),
          aes(fill = "No effect"), colour = "white", linewidth = 0.1) +
  scale_fill_manual(name = "Variance explained",
                    values = c("darkgray", inferno(5)[5], inferno(5)[4], inferno(5)[3], inferno(5)[2], inferno(5)[1]), 
                    breaks = c("No effect", "[0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]"),
                    labels = c("No effect", "0-0.25", "0.25-0.5", "0.5-0.75","0.75-1")) +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),  # Tamaño aún más reducido del cuadro de color
        legend.key.width = unit(1.2, "cm"), #nuevo
        legend.key.height = unit(0.5, "cm"),#nuevo
        legend.text = element_text(size = 14),  # Tamaño de texto reducido
        legend.title = element_text(size = 13, face = "bold")) +  # Tamaño de título reducido
  
  
  guides(fill = guide_legend(nrow=1,
                             title.position = "top", 
                             label.position = "bottom",
                             keyheight = unit(0.4, "cm"),  # Altura de la leyenda aún más reducida
                             keywidth = unit(4, "cm")))  # Anchura de la leyenda aún más reducida

# Guarda la imagen
ggsave(paste0(dir_out, "best_r2_",product,"_",transf,".png"), f1, height = 10, width = 15, unit = "in",
       bg = "white", dpi = 600)


####################
# BEST_MODEL ECORREGIONS
#######################

# best_model ecoregions

# Unir a los datos vectoriales
eco <- mutate(eco, best_mod = ifelse(is.na(best_mod) & mask == 1, "NoEffect", as.character(best_mod)))
# Aquí, la categoría "NoEffect" se añade a `best_mod`
eco$best_mod <- factor(eco$best_mod, levels = c("NoEffect", "1", "2", "3", "4", "5"))
# Filtrar para sólo las observaciones que tienen un 'best_mod'
eco_only <- filter(eco, !is.na(best_mod))

# Generate the color palette
#colors = c('darkgray','#420075', '#8F00FF','#2edaec','#ff0000','#FFDB15')
# colors = c('darkgray', '#D1AFFF', '#6A00A7', '#008080', '#A32D2D', '#6A0000')
colors = c('darkgray','#420075', '#8F00FF','#2edaec','#ff0000','#FFDB15')
           
f2 <- ggplot(data = eco_only, aes(fill = best_mod)) +
  geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
  geom_sf(colour = "white", linewidth = 0.01) +
  scale_fill_manual(values = colors,
                    breaks = c("NoEffect", "1", "2", "3", "4", "5"),
                    labels = c("No Effect", "SPEIa-SFWIc", "SPEIa-SPEIc", "SPEIa", "SFWIc", "SPEIc")) +  # Utiliza los nombres de los modelos aquí
  labs(fill = "Best Model") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),  # Cambia esto para ajustar el tamaño total de las claves
        legend.key.width = unit(1.2, "cm"), # Aumentar el ancho de las claves
        legend.key.height = unit(0.5, "cm"), # Aumentar la altura de las claves
        legend.text = element_text(size = 14),  # Aumentar el tamaño del texto
        legend.title = element_text(size = 13, face = "bold"))  # Aumentar el tamaño y estilo del título de la leyenda

# Guardar el mapa como PNG

ggsave(paste0(dir_out, "best_model_",product,"_mswep_",transf,".png"), f2, height = 10, width = 15, unit = "in",
       bg = "white", dpi = 600)
# ggsave(paste0(dir_out, "best_model_modis_mswep.pdf"), plot = g, device = "pdf", width = 10, height = 7)






############################################# 
    #SPEI COEF
#############################################

# Define the colors and labels
col_prob_hex <- c(
  "darkgray",  # For no effect
  "#001F3F",   # New color for < -1
  "#053061",
  "#2166ac",
  "#4393c3",
  "#92c5de",
  "#d1e5f0",
  "#f7f7f7",
  "#fddbc7",
  "#f4a582",
  "#d6604d",
  "#b2182b",
  "#67001f"   # Color for > 1
)

legend_labels <- c(
  "No effect",
  "< -1",
  "[-1, -0.8)",
  "[-0.8, -0.6)",
  "[-0.6, -0.4)",
  "[-0.4, -0.2)",
  "[-0.2, 0)",
  "[0, 0.2)",
  "[0.2, 0.4)",
  "[0.4, 0.6)",
  "[0.6, 0.8)",
  "[0.8, 1]",
  "> 1"
)

summary(best_spei_coef)
length(which(!is.na(best_spei_coef)))/length(which(mask==1))
sum(area_eco[!is.na(best_spei_coef)])/sum(area_eco[mask == 1])



# Classify the best_spei_coef values
eco <- mutate(eco, best_spei_coef)
eco_only <- eco %>%
  filter(!is.na(best_spei_coef)) %>%
  mutate(
    discrete_spei_coef = case_when(
      best_spei_coef < -1            ~ legend_labels[2],
      between(best_spei_coef, -1, -0.8)  ~ legend_labels[3],
      between(best_spei_coef, -0.8, -0.6) ~ legend_labels[4],
      between(best_spei_coef, -0.6, -0.4) ~ legend_labels[5],
      between(best_spei_coef, -0.4, -0.2) ~ legend_labels[6],
      between(best_spei_coef, -0.2, 0)    ~ legend_labels[7],
      between(best_spei_coef, 0, 0.2)     ~ legend_labels[8],
      between(best_spei_coef, 0.2, 0.4)   ~ legend_labels[9],
      between(best_spei_coef, 0.4, 0.6)   ~ legend_labels[10],
      between(best_spei_coef, 0.6, 0.8)   ~ legend_labels[11],
      best_spei_coef >= 0.8              ~ legend_labels[12],
      TRUE                               ~ legend_labels[1]   # NA or other cases
    )
  )

# Creating an empty data frame with all potential classes
all_classes_df <- data.frame(discrete_spei_coef = legend_labels)

# Plotting
f3 <- ggplot(eco_only) +
  # This layer ensures all levels are present but doesn't actually plot anything
  geom_blank(data = all_classes_df, aes(fill = discrete_spei_coef)) +
  
  # The rest of the plotting code remains the same
  geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
  geom_sf(aes(fill = discrete_spei_coef), colour = "white", linewidth = 0.01) +
  geom_sf(data = subset(eco, is.na(best_spei_coef) & mask == 1), fill = "darkgray", colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = col_prob_hex, breaks = legend_labels, labels = legend_labels, drop = FALSE) +
  labs(fill = "Antecedent coeficients") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 13, face = "bold")) +
  guides(fill = guide_legend(nrow = 2))

# Save the plot
ggsave(paste0(dir_out, "spei_coef_best_model_",product,"_mswep_", transf, ".png"), f3, height = 10, width = 15, unit = "in", bg = "white", dpi = 600)





############################################# 
    #SFWI    
#############################################

summary(best_sfwi_coef[best_sfwi_coef>0])
length(which(!is.na(best_sfwi_coef[best_sfwi_coef>0])))/length(which(mask==1))
sum(area_eco[!is.na(best_sfwi_coef[best_sfwi_coef>0])])/sum(area_eco[mask == 1])

summary(best_sfwi_coef[best_sfwi_coef<0])
length(which(!is.na(best_sfwi_coef[best_sfwi_coef<0])))/length(which(mask==1))
sum(area_eco[!is.na(best_sfwi_coef[best_sfwi_coef<0])])/sum(area_eco[mask == 1])


length(which(!is.na(best_sfwi_coef)))/length(which(mask==1))
sum(area_eco[!is.na(best_sfwi_coef)])/sum(area_eco[mask == 1])

# Classify the best_sfwi_coef values
eco <- mutate(eco, best_sfwi_coef)
eco_only <- eco %>%
  filter(!is.na(best_sfwi_coef)) %>%
  mutate(
    discrete_sfwi_coef = case_when(
      best_sfwi_coef < -1            ~ legend_labels[2],
      between(best_sfwi_coef, -1, -0.8)  ~ legend_labels[3],
      between(best_sfwi_coef, -0.8, -0.6) ~ legend_labels[4],
      between(best_sfwi_coef, -0.6, -0.4) ~ legend_labels[5],
      between(best_sfwi_coef, -0.4, -0.2) ~ legend_labels[6],
      between(best_sfwi_coef, -0.2, 0)    ~ legend_labels[7],
      between(best_sfwi_coef, 0, 0.2)     ~ legend_labels[8],
      between(best_sfwi_coef, 0.2, 0.4)   ~ legend_labels[9],
      between(best_sfwi_coef, 0.4, 0.6)   ~ legend_labels[10],
      between(best_sfwi_coef, 0.6, 0.8)   ~ legend_labels[11],
      best_sfwi_coef >= 0.8              ~ legend_labels[12],
      TRUE                               ~ legend_labels[1]   # NA or other cases
    )
  )

# Creating an empty data frame with all potential classes
all_classes_df <- data.frame(discrete_sfwi_coef = legend_labels)

# Plotting
f4 <- ggplot(eco_only) +
  # This layer ensures all levels are present but doesn't actually plot anything
  geom_blank(data = all_classes_df, aes(fill = discrete_sfwi_coef)) +
  
  # The rest of the plotting code remains the same
  geom_sf(data = wld, fill = "grey90", colour = "white", linewidth = 0.1) +
  geom_sf(aes(fill = discrete_sfwi_coef), colour = "white", linewidth = 0.01) +
  geom_sf(data = subset(eco, is.na(best_sfwi_coef) & mask == 1), fill = "darkgray", colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = col_prob_hex, breaks = legend_labels, labels = legend_labels, drop = FALSE) +
  labs(fill = "Antecedent coeficients") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 13, face = "bold")) +
  guides(fill = guide_legend(nrow = 2))

# Save the plot
ggsave(paste0(dir_out, "sfwi_coef_best_",product,"_modis_mswep_", transf, ".png"), f4, height = 10, width = 15, unit = "in", bg = "white", dpi = 600)


############################################# 
#SPEI m
#####################################
# Definir los valores para los breaks en orden inverso
brk_prob_reversed <- rev(seq(-14,-2))

# Paleta de colores similar a "Spectral", pero en hexadecimal
spectral_13_hex <-
  c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F28E2B",
    "#F46D43",
    "#D94801",
    "#D53E4F",
    "#9E0142")

# Combinar los colores Spectral con el color 'No effect'
combined_colors_reversed <-
  c("No effect" = "darkgray",
    setNames(spectral_13_hex, as.character(brk_prob_reversed)))

# Procesamiento previo de los datos
best_m_spei[is.na(best_spei_coef)] = NA
eco <- mutate(eco, best_m_spei )
eco_only <- filter(eco,!is.na(best_m_spei))
eco_only$best_m_spei <- as.factor(eco_only$best_m_spei)

# Generar el gráfico
f5 <- ggplot(eco_only) +
  geom_sf(
    data = wld,
    fill = "grey90",
    colour = "white",
    linewidth = 0.1
  ) +
  geom_sf(aes(fill = best_m_spei),
          colour = "white",
          linewidth = 0.01) +
  geom_sf(
    data = subset(eco, is.na(best_m_spei) & mask == 1),
    aes(fill = "No effect"),
    colour = "white",
    linewidth = 0.1
  ) +
  scale_fill_manual(
    values = combined_colors_reversed,
    breaks = c("No effect", brk_prob_reversed),
    labels = c("No effect", brk_prob_reversed)
  ) +
  labs(fill = "Best month") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 2))

# Guardar el gráfico
ggsave(
  paste0(dir_out, "best_m_spei_",product,"_mswep_",transf,".png"),
  f5,
  height = 10,
  width = 15,
  unit = "in",
  bg = "white",
  dpi = 600
)
    
#########################
#FWI m
#####################################
    
# Definir los valores para los breaks, empezando con "No effect", 0, -1, ..., -12
brk_prob_extended <- c("No effect", 0, seq(-1, -12))
# Crear un único punto y replicarlo para todos los niveles que quieres en tu leyenda
single_point <- st_point(c(0, 0))
fake_point <- st_sf(geometry = st_sfc(single_point), best_m_sfwi = factor("No effect", levels = brk_prob_extended))
fake_data <- fake_point[rep(1, length(brk_prob_extended)),]
fake_data$best_m_sfwi <- factor(brk_prob_extended, levels = brk_prob_extended)
# Asegúrate de que el CRS del objeto fake_data coincida con el de tus otros datos
st_crs(fake_data) <- st_crs(eco_only)
    
# Paleta de colores
spectral_14_hex <-
  c(
    "#4D3E82",
    "#7E2FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FFFFBF",
    "#FEE08B",
    "#FDAE61",
    "#F28E2B",
    "#F46D43",
    "#D94801",
    "#D53E4F",
    "#9E0142"
  )
combined_colors_extended <-
  setNames(c("darkgray", spectral_14_hex),
           as.character(brk_prob_extended))
# Procesamiento previo de los datos
best_m_sfwi[is.na(best_sfwi_coef)] = NA
eco <- mutate(eco, best_m_sfwi)
eco_only <- filter(eco,!is.na(best_m_sfwi))
eco_only$best_m_sfwi <- as.factor(eco_only$best_m_sfwi)

# Generar el gráfico
f6 <- ggplot() +
  geom_sf(
    data = wld,
    fill = "grey90",
    colour = "white",
    linewidth = 0.1
  ) +
  geom_sf(data = fake_data,
          aes(fill = best_m_sfwi),
          alpha = 0) +
  geom_sf(
    data = eco_only,
    aes(fill = best_m_sfwi),
    colour = "white",
    linewidth = 0.01
  ) +
  geom_sf(
    data = subset(eco, is.na(best_m_sfwi) & mask == 1),
    aes(fill = "No effect"),
    colour = "white",
    linewidth = 0.1
  ) +
  scale_fill_manual(values = combined_colors_extended,
                    breaks = brk_prob_extended,
                    labels = brk_prob_extended) +
  labs(fill = "Best month") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 2))

# Guardar el gráfico
ggsave(
  paste0(dir_out, "best_m_sfwi_",product,"_mswep_",transf,".png"),
  f6,
  height = 10,
  width = 15,
  unit = "in",
  bg = "white",
  dpi = 600
)

############################################ 
#SPEI t
###########################################

# Elegir 3 colores en formato hexadecimal de la paleta Dark2
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Nuevas categorías
new_brk_prob <- c("3", "6", "12")

# Combinar el color 'darkgray' con los colores seleccionados de la paleta Dark2
new_combined_colors <-
  c("No effect" = "darkgray",
    setNames(selected_Dark2_hex, as.character(new_brk_prob)))

# Procesamiento previo de los datos
best_t_spei[is.na(best_spei_coef)] <- NA
eco <- mutate(eco, best_t_spei)

# Cambiar los valores 1, 2, 3 a 3, 6, 12 respectivamente
eco <- eco %>%
  mutate(best_t_spei = case_when(
    best_t_spei == 1 ~ 3,
    best_t_spei == 2 ~ 6,
    best_t_spei == 3 ~ 12,
    TRUE ~ best_t_spei
  ))


# Filtrar los datos para eliminar los NA
eco_only <- filter(eco,!is.na(best_t_spei))
eco_only$best_t_spei <- as.factor(eco_only$best_t_spei)

# Generar el gráfico
f7 <- ggplot(eco_only) +
  geom_sf(
    data = wld,
    fill = "grey90",
    colour = "white",
    linewidth = 0.1
  ) +
  geom_sf(aes(fill = best_t_spei),
          colour = "white",
          linewidth = 0.01) +
  geom_sf(
    data = subset(eco, is.na(best_t_spei) & mask == 1),
    aes(fill = "No effect"),
    # Agregar el mapeo estético para la leyenda
    colour = "white",
    linewidth = 0.1
  ) +
  scale_fill_manual(
    values = new_combined_colors,
    breaks = c("No effect", new_brk_prob),
    labels = c("No effect", new_brk_prob)
  ) +
  labs(fill = "Best time scale (month)") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold")
  )

# Guardar el gráfico
ggsave(
  paste0(dir_out, "best_t_spei_",product,"_mswep_",transf,".png"),
  f7,
  height = 10,
  width = 15,
  unit = "in",
  bg = "white",
  dpi = 600
)

    
########################################### 
#SFWI t
###########################################   
    
# Procesamiento previo de los datos
best_t_sfwi[is.na(best_sfwi_coef)] <- NA
eco <- mutate(eco, best_t_sfwi)

# Cambiar los valores 1, 2, 3 a 3, 6, 12 respectivamente
eco <- eco %>%
  mutate(best_t_sfwi = case_when(
    best_t_sfwi == 1 ~ 3,
    best_t_sfwi == 2 ~ 6,
    best_t_sfwi == 3 ~ 12,
    TRUE ~ best_t_sfwi
  ))

# Filtrar los datos para eliminar los NA
eco_only <- filter(eco,!is.na(best_t_sfwi))
eco_only$best_t_sfwi <- as.factor(eco_only$best_t_sfwi)

# Generar el gráfico
f8 <- ggplot(eco_only) +
  geom_sf(
    data = wld,
    fill = "grey90",
    colour = "white",
    linewidth = 0.1
  ) +
  geom_sf(aes(fill = best_t_sfwi),
          colour = "white",
          linewidth = 0.01) +
  geom_sf(
    data = subset(eco, is.na(best_t_sfwi) & mask == 1),
    aes(fill = "No effect"),
    # Agregar el mapeo estético para la leyenda
    colour = "white",
    linewidth = 0.1
  ) +
  scale_fill_manual(
    values = new_combined_colors,
    breaks = c("No effect", new_brk_prob),
    labels = c("No effect", new_brk_prob)
  ) +
  labs(fill = "Best time scale (month)") +
  coord_sf(crs = "ESRI:54030") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold")
  )

# Guardar el gráfico
ggsave(
  paste0(dir_out, "best_t_sfwi_",product,"_mswep_",transf,".png"),
  f8,
  height = 10,
  width = 15,
  unit = "in",
  bg = "white",
  dpi = 600
)



#################
### SAVE FIGURES 
###############
library(cowplot)
library(patchwork)

###FIGURE 1
# Combinar las gráficas con plot_grid
plot_grid <- plot_grid(f1, f2, ncol = 1)
# Guardar la grilla de gráficas en un archivo PDF
ggsave(paste0(dir_out,"Figure1-",product,".pdf"), plot_grid, width = 15, height = 20, units = "in")
# Guardar la grilla de gráficas en un archivo PNG
ggsave(paste0(dir_out,"Figure1-",product,".png"), plot_grid, width = 15, height = 20, units = "in")


###FIGURE 2 CC 
# Combinar las gráficas con plot_grid
plot_grid2 <- plot_grid(f4, f6,f8,  ncol = 1)
# Guardar la grilla de gráficas en un archivo PDF
ggsave(paste0(dir_out,"Figure2-",product,".pdf"), plot_grid2, width = 15, height = 20, units = "in")
# Guardar la grilla de gráficas en un archivo PNG
ggsave(paste0(dir_out,"Figure2-",product,".png"), plot_grid2, width = 15, height = 20, units = "in")

###FIGURE 3 AC
# Combinar las gráficas con plot_grid
plot_grid3 <- plot_grid(f3, f5,f7,  ncol = 1)
# Guardar la grilla de gráficas en un archivo PDF
ggsave(paste0(dir_out,"Figure3-",product,".pdf"), plot_grid3, width = 15, height = 20, units = "in")
# Guardar la grilla de gráficas en un archivo PNG
ggsave(paste0(dir_out,"Figure3-",product,".png"), plot_grid3, width = 15, height = 20, units = "in")


