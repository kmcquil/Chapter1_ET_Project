#################################################################################################
## Process the HAND .tif from GEE to match the MODIS data and also to match the landsat data 
##################################################################################################
library(raster)
library(rgdal)
library(gdalUtils)

# base directory
home <- "G:/My Drive/Chapter1_ET_Project/Data"
source(paste0(home, "/Processing_scripts/helpers.R"))

# elevation input 
fin <- paste0(home, "/Topography/HeightAboveNearestDrainage/hand30_100.tif")
roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp")
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))

# MODIS
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif"))
fout <- paste0(home,"/Topography/HeightAboveNearestDrainage/hand_modis.tif" )
warpMn(fin, modis_template, fout, 'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)

# Landsat 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
fout <- paste0(home,"/Topography/HeightAboveNearestDrainage/hand_landsat.tif" )
warpMn(fin, landsat_template, fout,'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)
