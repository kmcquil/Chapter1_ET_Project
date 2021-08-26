#################################################################################
# resample SPI to MODIS and landsat resolutions 
#################################################################################
library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(doParallel)

UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)


home <- "G:/My Drive/Chapter1_ET_Project/Data" # base directory
#home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
source(paste0(home, "/Processing_scripts/helpers.R")) # source a helper function 

roi <- paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp")
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))


spi <- list.files(paste0(home, "/SPI/SPI1Y_OG"), full.names = T, pattern = ".tif$") # spi files 
spi_short <- list.files(paste0(home, "/SPI/SPI1Y_OG"), full.names = F, pattern = ".tif$")


# first do MODIS 
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) #MODIS template to resample data to

foreach(i = 1:length(spi))%dopar%{
  library(raster)
  library(rgdal)
  library(gdalUtils)
  source(paste0(home, "/Processing_scripts/helpers.R"))
  
  fout <- paste0(home, "/SPI/SPI1Y_MODIS/", spi_short[i])
  warpMn(spi[i], modis_template, fout, 'bilinear')
  r <- raster(fout)
  r <- mask(r, roi_mask)
  writeRaster(r, fout, overwrite = T)
  
  
}





