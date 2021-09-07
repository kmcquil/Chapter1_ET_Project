#################################################################################
# resample SPI to MODIS and landsat resolutions 
#################################################################################
library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(doParallel)

# initialize a cluster 
UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data" # set a home directory 
source(paste0(home, "/Processing_scripts/helpers.R")) # source a helper function 

roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp") # the file path of the ROI 
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))

spi <- list.files(paste0(home, "/SPI/SPI90_OG"), full.names = T, pattern = ".tif$") # list spi files 
spi_short <- list.files(paste0(home, "/SPI/SPI90_OG"), full.names = F, pattern = ".tif$")

# first resample to MODIS using a template of final MODIS data  
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) #MODIS template to resample data to
foreach(i = 1:length(spi))%dopar%{
  library(raster)
  library(rgdal)
  library(gdalUtils)
  source(paste0(home, "/Processing_scripts/helpers.R"))
  fout <- paste0(home, "/SPI/SPI90_MODIS/", spi_short[i])
  warpMn(spi[i], modis_template, fout, 'bilinear')  # resample 
  r <- raster(fout)
  r <- mask(r, roi_mask) # mask to the ROI 
  writeRaster(r, fout, overwrite = T)
}


# resample to Landsat template  
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
foreach(i = 1:length(spi))%dopar%{
  library(raster)
  library(rgdal)
  library(gdalUtils)
  source(paste0(home, "/Processing_scripts/helpers.R"))
  fout <- paste0(home, "/SPI/SPI90_Landsat/", spi_short[i], ".tif")
  warpMn(spi[i], landsat_template, fout, 'bilinear') # resample 
  r <- raster(fout)
  r <- mask(r, roi_mask) # mask to the ROI 
  writeRaster(r, fout, overwrite = T)
}




