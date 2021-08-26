library(gdalUtils)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)


home <- "G:/My Drive/Chapter1_ET_Project/Data"

# break each raster into 9 tiles 
# do x first (ncols) = 15096
# xwin <- c(0, 5032, 10064)
# x_tile_size <- c(5032, 5032, 5032)
xwin <- seq(0, 15096, 2516)
x_tile_size <-rep(2516, 6)


# next do y (nrows) = 11501
# ywin <- c(0, 3833, 7666)
# y_tile_size <- c(3833, 3833, 3835)
ywin <- seq(0, 11501, 1643)
y_tile_size <- rep(1643, 7)


fin <- paste0(home, "/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
  
for(k in 1:length(x_tile_size)){
  for(j in 1:length(y_tile_size)){
      fout <- paste0(home, "/landcover/LANDSAT_FOREST/retile/", "x", sprintf("%05d",xwin[k]), "y", sprintf("%05d", ywin[j]), ".tif")
      gdal_translate(fin, fout, srcwin=c(xwin[k], ywin[j], x_tile_size[k], y_tile_size[j]))
  }
}
  





