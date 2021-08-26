# retile the long term mean and standard deviation rasters so that I can check those out before I start rerunning everything 
library(raster)
library(gdalUtils)
library(data.table)
library(rgdal)
library(lubridate)
library(doParallel)
library(foreach)

UseCores <- 5
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
ltm <- list.files(paste0(home, "/Landsat_ET/ltm_month"), full.names = T) # bring in ltm and ltsd to index for anom calcs
ltsd <- list.files(paste0(home, "/Landsat_ET/ltsd_month"), full.names = T)
template <- paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif")


idx <- c("_01", "_02", "_03", "_04", "_05", "_06", "_07", "_08", "_09", "_10", "_11", "_12")
foreach(i = 1:length(idx))%dopar%{
  library(gdalUtils)
  library(raster)
  
  m <- ltm[grep(idx[i], ltm)]
  sdd <- ltsd[grep(idx[i], ltsd)]
  
  out_m <- paste0("/share/klmarti3/kmcquil/Chapter1_ET_Project/ltm", idx[i], ".tif")
  out_sd <- paste0("/share/klmarti3/kmcquil/Chapter1_ET_Project/ltsd", idx[i], ".tif")
  
  align_rasters(m, template, dstfile = out_m, overwrite = T)
  align_rasters(sdd, template, dstfile = out_sd, overwrite = T)
  
}



rasts <- list.files("/share/klmarti3/kmcquil/Chapter1_ET_Project", full.names = T, pattern = ".tif$")

plot(rasts[1])







