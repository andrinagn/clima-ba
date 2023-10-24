#################################################################################################
#		
# Title:	The interannual variability of global burned area is mostly explained by climatic 
#               drivers
#
# 				
# Authors: 	Andrina Gincheva, University of Murcia (andrina@um.es)
#               Marco Turco, University of Murcia (marco.turco@um.es) 
#               
#
#################################################################################################

#################################################################################################
# A. General instructions 
#################################################################################################

This project is designed to be executed with shell scripts and R codes. 
Execute script files in the order they are listed.

Data sources:

- Burned area data from MCD64CMQ:
sftp://fuoco.geog.umd.edu/data/MODIS/C6/MCD64CMQ

- Burned area data from FIRECCI51:
https://catalogue.ceda.ac.uk/uuid/58f00d8814064b79a0c49662ad3af537

- ERA5 variables:
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels

- Fire Weather index:
https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-fire-historical-v1

- Precipitation MSWEP:
https://www.gloh2o.org/mswep/

- Ecoregions:
https://ecoregions.appspot.com/

Notes regarding reproducibility:

Script files starting with "1_" in their name are for data preprocessing. 
Most of these script files will NOT run because we do not include the 
raw data because files are simply too large to conveniently share. 
We suggest you run scripts starting with "2_" which directly reproduces the 
results in the paper. 

If you have any questions or wish to express any comment to the authors, please 
contact Andrina Gincheva and/or Dr. Marco Turco at the emails indicated above.


#################################################################################################
# B. Description of script files
#################################################################################################

Scripts for data preparation

- 1_1_hdf2netcdf_mcd64cmq.R
Read hdf MCD64CMQ and transform to netcdf.

- 1_2_MCD64CMQ2ecoregions.R
Prepare ecoregion-level MCD64CMQ dataset.

- 1_3_FIRECCI512ecoregions.R
Prepare ecoregion-level FIRECCI51 dataset.

- 1_4_fwi2ecoregions.R
Prepare ecoregion-level SWFI dataset.

- 1_5_spei2ecoregions.R
Prepare ecoregion-level SPEI dataset.

- 1_6_fire_season.R
Estimate the fire season as in Abatzoglou et al. (2018).


Scripts to reproduce the results in the paper.

- 2_1_plot_fireseason.R
Plot figures relative to the fire season (Supp. Material)

- 2_2_out_of_sample_model.R
Determine the potential predictors for the climate-fire model through out-of-sample calibration. The process may require several hours.

- 2_3_in_sample_model.R
Determine the best models checking for the significance of the predictors, residual assumptions, and multicollinearity. The process is very fast.

- 2_4_example_figure1.R
Compute fire season and climate-fire model for a specific ecoregion (Sierra Nevada, California). 
Plot figures to make Figure 1.

- 2_5_plot_maps.R
Plot maps for Figures 2-4.

- 2_6_plot_scatterplots_boxplots.R
Plot scatterplots for Figure 5 and box plots for Figure 6.