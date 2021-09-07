#################################################################################
## Retile the spi to match the landsat retile 
#################################################################################
library(gdalUtils)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)

# initialize cluster 
UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)
print(cl)


home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
efile <- list.files(paste0(home, "/SPI/SPI90_Landsat"), full.names = T) # get the filepaths for spi rasters
efile_short <- list.files(paste0(home, "/SPI/SPI90_Landsat"), full.names = F) # get the filepaths for spi rasters

print(length(efile))

# break each raster into 42 tiles 
xwin <- seq(0, 15096, 2516)
x_tile_size <-rep(2516, 6)


# next do y (nrows) = 11501
ywin <- seq(0, 11501, 1643)
y_tile_size <- rep(1643, 7)

# for every spi raster, retile and save with tile + date as ID 
foreach(i = 1:length(efile))%dopar%{
  library(gdalUtils)
  library(raster)
  library(rgdal)
  fin <- efile[i]
  finshort <- substr(efile_short[i], 1, 10)
  
  for(k in 1:length(x_tile_size)){
    for(j in 1:length(y_tile_size)){
      fout <- paste0(home, "/SPI/retile/", "x", sprintf("%05d",xwin[k]), "y", sprintf("%05d", ywin[j]), "_", finshort)
      
      gdal_translate(fin, fout, srcwin=c(xwin[k], ywin[j], x_tile_size[k], y_tile_size[j]))
    }
  }
  
}





