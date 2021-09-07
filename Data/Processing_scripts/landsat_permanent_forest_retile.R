########################################################################################################
## Retile the Landsat permanent forest raster to match the retiled Landsat ET 
## This is necessary to reduce memory requirements and speed up computation 
#########################################################################################################
library(gdalUtils)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)

home <- "G:/My Drive/Chapter1_ET_Project/Data" # set a home directory 

# break the raster into 42 tiles 
# do x first (ncols) = 15096 so we divide 15096 by 6
xwin <- seq(0, 15096, 2516)
x_tile_size <-rep(2516, 6)

# next do y (nrows) = 11501 so we divide 11501 by 7 
ywin <- seq(0, 11501, 1643)
y_tile_size <- rep(1643, 7)

# use gdal to break the raster into tiles 
fin <- paste0(home, "/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif") # input landsat forest file
for(k in 1:length(x_tile_size)){
  for(j in 1:length(y_tile_size)){
      fout <- paste0(home, "/landcover/LANDSAT_FOREST/retile/", "x", sprintf("%05d",xwin[k]), "y", sprintf("%05d", ywin[j]), ".tif")
      gdal_translate(fin, fout, srcwin=c(xwin[k], ywin[j], x_tile_size[k], y_tile_size[j]))
  }
}
  





