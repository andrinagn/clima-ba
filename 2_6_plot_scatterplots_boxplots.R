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
library(broom)
library(ggplot2)
library(gridExtra)
library(mgcv)


### fix parameters ANDRINA VENTISCA
# dir_out = '/home/andrina/Dropbox/model/out_ecoregions/'
# dir_data = '/home/andrina/Desktop/climate-fire/datos/'
# dir_shp='/home/andrina/Desktop/climate-fire/datos/shp/'

### fix parameters A
dir_out = '~/Dropbox/model/out_ecoregions/'
dir_data = '~/Documents/dati/fire_climate_data/climate_fire/datos/'
dir_shp='/Users/marco/Documents/virtualBox/temiav/Ecoregions2017/'

### version
transf="log" #"log" "radq"

if (transf=="log") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_log_fire_season_ja_2018')  
} else if (transf=="radq") {
  list_version=c('spei_sfwi_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018','spei_spei_modis_mswep_era5_fireseason_ecoregions_fire_season_ja_2018')
}
list_version_names=c('spei-sfwi','spei-spei')
list_models=c('spei-sfwi','spei-spei','ANTspei','CONsfwi','CONspei')

load(paste0(dir_out,
                   "mask_ecoregions_season_ja_2018.RData"))

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

## climatologies
load(paste0(dir_out, "AI_ERA5_MSWEP_2001_2021_ecoregions.RData")) #AI
load(paste0(dir_out, "PREC_ERA5_MSWEP_2001_2021_ecoregions_annual.RData")) #AI
load(paste0(dir_out, "PET_ERA5_MSWEP_2001_2021_ecoregions_annual.RData")) #AI
# load(paste0(dir_out, "TAS_ERA5_2001_2021_ecoregions.RData")) #AI


# Assuming your data frame is named df
df <- data.frame(
  AI = AI, 
  PREC = PREC,
  PET = PET,
  # TAS = TAS,
  best_r2 = best_r2,
  best_sfwi_coef = abs(best_sfwi_coef),
  best_m_sfwi = best_m_sfwi,
  best_t_sfwi = best_t_sfwi,
  best_spei_coef = best_spei_coef,
  # best_sfwi_coef = best_sfwi_coef,
  best_m_spei = best_m_spei,
  best_t_spei = best_t_spei
)


# A function to extract and format the correlation coefficient based on p-value
get_cor_p <- function(x, y) {
  test <- cor.test(x, y, method = "spearman")
  test <- cor.test(x, y, method = "pearson")
  cor_val <- round(test$estimate, 2) # rounding to 2 decimal places
  if (test$p.value < 0.01) {
    return(paste0(cor_val, "**"))
  } else if (test$p.value < 0.05) {
    return(paste0(cor_val, "*"))
  } else {
    return(as.character(cor_val))
  }
}

# Variables of interest
best_cols <- grep("^best_", names(df), value = TRUE)
fixed_cols <- c("AI", "PREC", "PET")

# Compute the correlations
results <- sapply(fixed_cols, function(fixed_col) {
  sapply(best_cols, function(best_col) {
    get_cor_p(df[[fixed_col]], df[[best_col]])
  })
}, simplify = "matrix")

# Convert to a data frame for better readability
results_df <- as.data.frame(results)
rownames(results_df) <- best_cols

results_df







#######

df_sfwi <- df[, c("AI", "PREC", "PET", "best_sfwi_coef", "best_m_sfwi", "best_t_sfwi")]
df_sfwi <- df_sfwi[complete.cases(df_sfwi), ]

df_spei <- df[, c("AI", "PREC", "PET", "best_spei_coef", "best_m_spei", "best_t_spei")]
df_spei <- df_spei[complete.cases(df_spei), ]

###### plot



# Helper function to format R^2 and p-value
format_labels <- function(r2, p_value){
  asterisks <- ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))
  return(sprintf("%.1f%% %s", r2*100, asterisks))
}


# For p2: best_sfwi_coef vs PET
model_p2 <- lm(best_sfwi_coef ~ PET, data = df_sfwi)
glance_p2 <- glance(model_p2)
p2 <- ggplot(df_sfwi, aes(x = PET, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PET [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), label = "a)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels(glance_p2$r.squared, glance_p2$p.value), hjust = 1.1, vjust = 1.5)

# For p3: best_sfwi_coef vs PREC
model_p3 <- lm(best_sfwi_coef ~ PREC, data = df_sfwi)
glance_p3 <- glance(model_p3)
p3 <- ggplot(df_sfwi, aes(x = PREC, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PREC [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), label = "b)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels(glance_p3$r.squared, glance_p3$p.value), hjust = 1.1, vjust = 1.5)


# For p5: best_spei_coef vs PET
model_p5 <- lm(best_spei_coef ~ PET, data = df_spei)
glance_p5 <- glance(model_p5)
p5 <- ggplot(df_spei, aes(x = PET, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PET [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PET), y = max(df_spei$best_spei_coef), label = "c)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PET), y = max(df_spei$best_spei_coef), 
           label = format_labels(glance_p5$r.squared, glance_p5$p.value), hjust = 1.1, vjust = 1.5)

# For p6: best_spei_coef vs PREC
model_p6 <- lm(best_spei_coef ~ PREC, data = df_spei)
glance_p6 <- glance(model_p6)
p6 <- ggplot(df_spei, aes(x = PREC, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  labs(x = "PREC [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PREC), y = max(df_spei$best_spei_coef), label = "d)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PREC), y = max(df_spei$best_spei_coef), 
           label = format_labels(glance_p6$r.squared, glance_p6$p.value), hjust = 1.1, vjust = 1.5)

# Arrange plots
# plot_grid <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
plot_grid <- grid.arrange( p2, p3, p5, p6, ncol = 2)


ggsave(paste0(dir_out,"scatterplots_lm.pdf"), plot_grid, width = 11, height = 8.5, units = "in")
ggsave(paste0(dir_out,"scatterplots_lm.png"), plot_grid, width = 11, height = 8.5, units = "in")




####### GAM




# Helper function to format R^2 and p-value for GAM models
format_labels_gam <- function(model){
  r2 <- summary(model)$r.sq
  p_value <- coef(summary(model))[,4][1]  # Assuming you're interested in the p-value for the smooth term
  asterisks <- ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", ""))
  return(sprintf("%.1f%% %s", r2*100, asterisks))
}




df_r2 <- df[, c("AI", "PREC", "PET", "best_r2")]
df_r2 <- df_r2[complete.cases(df_r2), ]


model_r2 <- gam(best_r2 ~ s(PET), data = df_r2)
p1 <- ggplot(df_r2, aes(x = PET, y = best_r2)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PET [mm/y]", y = "Variance explained") +
  theme_minimal() 
 p1
summary(model_r2)

model_r2 <- gam(best_r2 ~ s(PREC), data = df_r2)
p1 <- ggplot(df_r2, aes(x = PREC, y = best_r2)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PREC [mm/y]", y = "Variance explained") +
  theme_minimal() 

p1
summary(model_r2)

# For p2: best_sfwi_coef vs PET
model_p2 <- gam(best_sfwi_coef ~ s(PET), data = df_sfwi)
summary(model_p2)
p2 <- ggplot(df_sfwi, aes(x = PET, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PET [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), label = "a)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PET), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels_gam(model_p2), hjust = 1.1, vjust = 1.5)

# For p3: best_sfwi_coef vs PREC
model_p3 <- gam(best_sfwi_coef ~ s(PREC), data = df_sfwi)
summary(model_p3)
p3  <- ggplot(df_sfwi, aes(x = PREC, y = best_sfwi_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PREC [mm/y]", y = "Concurrent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), label = "b)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_sfwi$PREC), y = max(df_sfwi$best_sfwi_coef), 
           label = format_labels_gam(model_p3), hjust = 1.1, vjust = 1.5)

# For p5: best_spei_coef vs PET
model_p5 <- gam(best_spei_coef ~ s(PET), data = df_spei)
summary(model_p5)
p5 <- ggplot(df_spei, aes(x = PET, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PET [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PET), y = max(df_spei$best_spei_coef), label = "c)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PET), y = max(df_spei$best_spei_coef), 
           label = format_labels_gam(model_p5), hjust = 1.1, vjust = 1.5)

# For p6: best_spei_coef vs PREC
model_p6 <- gam(best_spei_coef ~ s(PREC), data = df_spei)
summary(model_p6)
p6 <- ggplot(df_spei, aes(x = PREC, y = best_spei_coef)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE,color="darkred") +
  labs(x = "PREC [mm/y]", y = "Antecedent Climate coefficients") +
  theme_minimal() +
  annotate("text", x = min(df_spei$PREC), y = max(df_spei$best_spei_coef), label = "d)", hjust = -0.1, vjust = 1.5) +
  annotate("text", x = max(df_spei$PREC), y = max(df_spei$best_spei_coef), 
           label = format_labels_gam(model_p6), hjust = 1.1, vjust = 1.5)

# Arrange plots
plot_grid <- grid.arrange(p2, p3, p5, p6, ncol = 2)
ggsave(paste0(dir_out,"scatterplots_gam.pdf"), plot_grid, width = 11, height = 8.5, units = "in")
ggsave(paste0(dir_out,"scatterplots_gam.png"), plot_grid, width = 11, height = 8.5, units = "in")


## boxplot t

########################
##ANTECEDENT CLIMATE
#######################

### BEST M PET
# Create color palette
spectral_13_hex <- c(
  
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
  "#9E0142"
)

# Especificamos el mapeo de los valores a los colores
color_mapping <- setNames(spectral_13_hex, c( -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14))

# Calculate counts, percentages, and median
counts <- df_spei %>%
  group_by(best_m_spei) %>%
  summarise(n = n(),
            median_PET = median(PET, na.rm = TRUE)) %>%
  mutate(n = round(100 * n / length(which(mask == 1))),
         label = paste0(n, "%"))

# Create boxplot
f1 <- ggplot(df_spei, aes(x = best_m_spei, y = PET, group = best_m_spei, fill = as.factor(best_m_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 610, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  scale_x_continuous(breaks = unique(df_spei$best_m_spei),
                     labels = as.character(unique(df_spei$best_m_spei))) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = FALSE) +
  labs(x = "Best month for AC", y = "PET [mm/y]") + 
  coord_cartesian(ylim = c(0, 620)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

# Print the plot
print(f1)

### BEST M PPREC



# Create color palette
spectral_13_hex <- c(
  
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
  "#9E0142"
)

# Especificamos el mapeo de los valores a los colores
color_mapping <- setNames(spectral_13_hex, c( -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14))

# Calculate counts, percentages, and median
counts <- df_spei %>%
  group_by(best_m_spei) %>%
  summarise(n = n(),
            median_PREC = median(PREC, na.rm = TRUE)) %>%
  mutate(n = round(100 * n / length(which(mask == 1))),
         label = paste0(n, "%"))

# Create boxplot
f2 <- ggplot(df_spei, aes(x = best_m_spei, y = PREC, group = best_m_spei, fill = as.factor(best_m_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 3500, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  scale_x_continuous(breaks = unique(df_spei$best_m_spei),
                     labels = as.character(unique(df_spei$best_m_spei))) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = FALSE) +
  labs(x = "Best month for AC", y = "PREC [mm/y]") + 
  coord_cartesian(ylim = c(0, 3500)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

# Print the plot
print(f2)



#BEST t PET 
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_spei$best_t_spei))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_spei", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PET for positioning labels
label_position <- round(quantile(df_spei$PET, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f3 <- ggplot(df_spei, aes(x=factor(best_t_spei), y=PET, fill=factor(best_t_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(x = factor(best_t_spei), y = 700, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  
  coord_cartesian(ylim = c(0, 700)) +
  scale_fill_manual(values = selected_Dark2_hex) +
  scale_x_discrete(labels = c("3", "6", "12")) +
  
  guides(fill=FALSE) +
  
  labs(x = "Best time scale for AC", y = "PET [mm/y]") + 
  theme_minimal()
f3 <- f3 + theme(
  axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
  axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20)
)

print(f3)



#BEST t PREC
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_spei$best_t_spei))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_spei", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PET for positioning labels
label_position <- round(quantile(df_spei$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f4 <- ggplot(df_spei, aes(x=factor(best_t_spei), y=PREC, fill=factor(best_t_spei))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(x = factor(best_t_spei), y = 3500, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  
  coord_cartesian(ylim = c(0, 3500)) +
  scale_fill_manual(values = selected_Dark2_hex) +
  scale_x_discrete(labels = c("3", "6", "12")) +
  
  guides(fill=FALSE) +
  
  labs(x = "Best time scale for AC", y = "PREC [mm/y]") + 
  theme_minimal()
f4 <- f4 + theme(
  axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
  axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20)
)

print(f4)


########################
## CONCURRENT CLIMATE 
#######################

# Create color palette
spectral_13_hex <- c(
  "#4D3E82",
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
  "#9E0142"
)

# Especificamos el mapeo de los valores a los colores
color_mapping <- setNames(spectral_13_hex, c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13))

# Calculate counts, percentages, and median
counts <- df_sfwi %>%
  group_by(best_m_sfwi) %>%
  summarise(n = n(),
            median_PET = median(PET, na.rm = TRUE)) %>%
  mutate(n = round(100 * n / length(which(mask == 1))),
         label = paste0(n, "%"))

# Create boxplot
f5 <- ggplot(df_sfwi, aes(x = best_m_sfwi, y = PET, group = best_m_sfwi, fill = as.factor(best_m_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 600, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  scale_x_continuous(breaks = unique(df_sfwi$best_m_sfwi),
                     labels = as.character(unique(df_sfwi$best_m_sfwi))) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = FALSE) +
  labs(x = "Best month for CC", y = "PET [mm/y]") + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

# Print the plot
print(f5)


### BEST M PPREC
# Calculate the number of samples for each value of best_m_sfwi

# Create color palette
spectral_13_hex <- c(
  "#4D3E82",
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
  "#9E0142"
)

# Especificamos el mapeo de los valores a los colores
color_mapping <- setNames(spectral_13_hex, c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13))

# Calculate counts, percentages, and median
counts <- df_sfwi %>%
  group_by(best_m_sfwi) %>%
  summarise(n = n(),
            median_PET = median(PREC, na.rm = TRUE)) %>%
  mutate(n = round(100 * n / length(which(mask == 1))),
         label = paste0(n, "%"))

# Create boxplot
f6 <- ggplot(df_sfwi, aes(x = best_m_sfwi, y = PREC, group = best_m_sfwi, fill = as.factor(best_m_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 4000, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  scale_x_continuous(breaks = unique(df_sfwi$best_m_sfwi),
                     labels = as.character(unique(df_sfwi$best_m_sfwi))) +
  scale_fill_manual(values = color_mapping) +
  guides(fill = FALSE) +
  labs(x = "Best month for CC", y = "PREC [mm/y]") + 
  coord_cartesian(ylim = c(0, 4200)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

# Print the plot
print(f6)


#BEST t PET 
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_sfwi$best_t_sfwi))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_sfwi", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Create the boxplot without outliers
f7 <- ggplot(df_sfwi, aes(x=factor(best_t_sfwi), y=PET, fill=factor(best_t_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(x = factor(best_t_sfwi), y = 600, label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  
  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  
  # Change x-axis labels
  scale_x_discrete(labels = c("3", "6", "12")) +
  
  # Remove the legend
  guides(fill=FALSE) +
  
  # Labels for the axis
  labs(x = "Best time scale for CC", y = "PET [mm/y]") + 
  
  # Define y-axis limits
  coord_cartesian(ylim = c(0, 600)) +
  
  # Customizing axis and text
  theme_minimal() +
 theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

print(f7)







#BEST t PREC
# Paleta de colores similar a "Dark2", pero en hexadecimal
selected_Dark2_hex <- c("#1B9E77", "#D95F02", "#7570B3")

# Calculate the number of samples for each class
counts <- as.data.frame(table(df_sfwi$best_t_sfwi))
counts$Freq <- round(100*counts$Freq/length(which(mask==1)))
colnames(counts) <- c("best_t_sfwi", "n")

# Convert counts to percentage and create the label
counts$label <- paste0(counts$n, "%")

# Calculate the 97.5th percentile of PET for positioning labels
label_position <- round(quantile(df_sfwi$PREC, 0.98, na.rm=TRUE))

# Create the boxplot without outliers
f8 <- ggplot(df_sfwi, aes(x=factor(best_t_sfwi), y=PREC, fill=factor(best_t_sfwi))) + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_label(data = counts, aes(y = 4000,  label = label), 
             position = position_dodge(width = 0.75),
             vjust = 0.5, fontface="bold", fill="grey99", label.size=0.5, size=6) +
  

  # Setting color palette for the boxes
  scale_fill_manual(values = selected_Dark2_hex) +
  # Change x-axis labels
  scale_x_discrete(labels = c("3", "6", "12")) +
  # Remove the legend
  guides(fill=FALSE) +
  labs(x = "Best time scale for CC", y = "PREC [mm/y]") + 
  coord_cartesian(ylim = c(0, 4000)) +
  # Customizing axis and text
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.title.y = element_text(face = "bold", size = 20),  # I increased the size for emphasis
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

print(f8)

########GUARDAR FIGURAS- combina ejes automáticamente
#install.packages("patchwork")
#install.packages("cowplot")
library(patchwork)
library(cowplot)
# Arrange plots
#guardarlo de forma simple
#plot_grid <- f5 + f6 + f7 + f8 +f1+f2+f3+f4 + plot_layout(ncol = 2)

library(cowplot)


# Creamos una grilla con espacio entre las figuras
plot_grid <- plot_grid(
  f5, f6, f7, f8, f1, f2, f3, f4, 
  ncol = 2,
  rel_widths = c(25, 25),  # Puedes ajustar estos valores para cambiar el ancho relativo de las columnas
  rel_heights = c(25, 25, 25, 25)  # Puedes ajustar estos valores para cambiar el alto relativo de las filas
)


# Añadir márgenes a cada gráfica
margen <- margin(1.5, 1.5, 1.5, 1.5, "cm")

f1 <- f1 + theme(plot.margin = margen)
f2 <- f2 + theme(plot.margin = margen)
f3 <- f3 + theme(plot.margin = margen)
f4 <- f4 + theme(plot.margin = margen)
f5 <- f5 + theme(plot.margin = margen)
f6 <- f6 + theme(plot.margin = margen)
f7 <- f7 + theme(plot.margin = margen)
f8 <- f8 + theme(plot.margin = margen)

# Combinar las gráficas con plot_grid
plot_grid <- plot_grid(f5, f6, f7, f8, f1, f2, f3, f4, ncol = 2)


ggsave(paste0(dir_out,"boxplot_cc_ac.pdf"), plot_grid, width = 22, height = 32, units = "in")
# Guardar la grilla de gráficas en un archivo PNG
ggsave(paste0(dir_out,"boxplot_cc_ac.png"), plot_grid, width = 22, height = 32, units = "in")

