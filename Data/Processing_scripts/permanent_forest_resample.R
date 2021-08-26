############################################################################################
## resample forest cover masks to match landsat and modis 
############################################################################################
library(raster)
library(rgdal)
library(gdalUtils)

home <- "G:/My Drive/Chapter1_ET_Project/Data" # base directory
source(paste0(home, "/Processing_scripts/helpers.R")) # source a helper function 

roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp") # bring in the roi path 
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))

# first do MODIS 
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) #MODIS template to resample data to
forest_mask_modis <- paste0(home, "/landcover/MODIS_FOREST/modis_permanent_forest.tif")
fout <- paste0(home, "/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif")
warpMn(forest_mask_modis, modis_template, fout, 'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)



# next to landsat 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
forest_mask_landsat <- paste0(home, "/landcover/LANDSAT_FOREST/landsat_permanent_forest.tif")
fout <- paste0(home,  "/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif" )
warpMn(forest_mask_landsat, landsat_template, fout,'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)


