#################################################################################
## Landsat resample files are apparently too big so i need to retile now fuck me 
#################################################################################

library(gdalUtils)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)

UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)
print(cl)


home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
efile <- list.files(paste0(home, "/Landsat_ET/tifs_resampled"), full.names = T) # get the filepaths for et rasters
efile_short <- list.files(paste0(home, "/Landsat_ET/tifs_resampled"), full.names = F) # get the filepaths for et rasters

print(length(efile))

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

foreach(i = 1:length(efile))%dopar%{
  library(gdalUtils)
  library(raster)
  library(rgdal)
  fin <- efile[i]
  finshort <- efile_short[i]
  
  for(k in 1:length(x_tile_size)){
    for(j in 1:length(y_tile_size)){
      fout <- paste0(home, "/Landsat_ET/retile/", "x", sprintf("%05d",xwin[k]), "y", sprintf("%05d", ywin[j]), "_", finshort)
      
      gdal_translate(fin, fout, srcwin=c(xwin[k], ywin[j], x_tile_size[k], y_tile_size[j]))
    }
  }
  
}




