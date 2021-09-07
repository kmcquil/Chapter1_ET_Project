#################################################################################################
## Resample the elevation rasters from GEE to MODIS and Landsat resolution/extent/crs 
##################################################################################################
library(raster)
library(rgdal)
library(gdalUtils)
home <- "G:/My Drive/Chapter1_ET_Project/Data"  # home directory
source(paste0(home, "/Processing_scripts/helpers.R"))  # script has functions to help resample

# elevation input 
fin <- paste0(home, "/Topography/usgsNED_elevation/elevation10m.tif")
roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp")
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))

# resmaple to match MODIS 
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) # template with parameters for resampling
fout <- paste0(home,"/Topography/usgsNED_elevation/elevation500m.tif" ) # file name for output
warpMn(fin, modis_template, fout, 'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask) # mask the new raster to the ROI 
writeRaster(r, fout, overwrite = T) # resave 

# resample to match Landsat 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
fout <- paste0(home,"/Topography/usgsNED_elevation/elevation30m.tif" )
warpMn(fin, landsat_template, fout, 'bilinear')
r <- raster(fout)
r <- mask(r, roi_mask)
writeRaster(r, fout, overwrite = T)

