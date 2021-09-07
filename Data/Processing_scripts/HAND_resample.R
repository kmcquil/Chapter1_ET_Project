#################################################################################################
## Process the HAND rasters from GEE to match MODIS and Landsat
##################################################################################################
library(raster)
library(rgdal)
library(gdalUtils)

# base directory
home <- "G:/My Drive/Chapter1_ET_Project/Data"
source(paste0(home, "/Processing_scripts/helpers.R"))

# HAND input
fin <- paste0(home, "/Topography/HeightAboveNearestDrainage/hand30_100.tif")
roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp")
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))

# resample to match MODIS and mask to the ROI
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif"))
fout <- paste0(home,"/Topography/HeightAboveNearestDrainage/hand_modis.tif" )
warpMn(fin, modis_template, fout, 'bilinear') 
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)

# resample to match Landsat and mask to the ROI 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
fout <- paste0(home,"/Topography/HeightAboveNearestDrainage/hand_landsat.tif" )
warpMn(fin, landsat_template, fout,'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)
