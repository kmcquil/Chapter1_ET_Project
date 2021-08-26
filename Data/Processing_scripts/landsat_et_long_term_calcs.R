####################################################################################################
## Calculate monthly landsat ET long term mean and standard deviation 
####################################################################################################

# gotta do it for each tile
library(raster)
library(gdalUtils)
library(data.table)
library(rgdal)
library(lubridate)
library(doParallel)
library(foreach)

#home <- "G:/My Drive/Chapter1_ET_Project/Data"
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"

## create data.table to easily subset all of the files by date and tile
efile <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = T) # get the filepaths for et rasters
efile_short <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = F) # get the raster names 
file_dt <- data.table(file = efile, 
                      file_short = efile_short, 
                      date = as.Date(substr(efile_short, 17, 24), "%Y%m%d")) # get the date of files in a data.table

file_dt$month <- month(file_dt$date) # find month and year to help with subsetting 
file_dt$year <- year(file_dt$date)
file_dt$YM <- paste0(file_dt$year, sprintf("%02d", file_dt$month))
## Add in column for tile 
file_dt$tile <- substr(efile_short, 1, 12)

unique_months <- unique(file_dt$month) # for the long term calcs, we will calculate for each month 

unique_tiles <- unique(file_dt$tile)

beginCluster(32)

LTs <- function(TD){
  tile_dt <- file_dt[tile == TD]
  
  
  for(i in 1:length(unique_months)){
    sub <- tile_dt[month == unique_months[i]]$file
    stk <- do.call('stack', lapply(sub, raster))
    print("finished stacking")
    print(nlayers(stk))
    avg_stk <- clusterR(stk, fun = calc, args =list(fun = mean, na.rm = T))
    sd_stk <- clusterR(stk, fun = calc, args = list(fun=sd, na.rm=T))
    
    writeRaster(sd_stk, paste0(home, "/Landsat_ET/ltsd_month/", TD, "_",sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    
    writeRaster(avg_stk, paste0(home, "/Landsat_ET/ltm_month/et_",TD, "_", sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    print(i)
  }
  
  
}


LTs(unique_tiles[1])
#lapply(unique_tiles, LTs)

endCluster()