library(raster)
library(gdalUtils)
library(data.table)
library(rgdal)
library(lubridate)
library(doParallel)
library(foreach)


#home <- "G:/My Drive/Chapter1_ET_Project/Data"
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"


## calculate long term mean and standard deviation on a pixel basis 
efile <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = T) # get the filepaths for et rasters
efile_short <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = F) # get the raster names 
file_dt <- data.table(file = efile, 
                      file_short = efile_short, 
                      date = as.Date(substr(efile_short, 17, 24), "%Y%m%d")) # get the date of files in a data.table

file_dt$month <- sprintf("%02d", month(file_dt$date)) # find month and year to help with subsetting 
file_dt$year <- year(file_dt$date)
file_dt$YM <- paste0(file_dt$year, file_dt$month)
file_dt$tile <- substr(efile_short, 1,12)

unique_YM <- unique(file_dt$YM) # get a list of unique year-month combos

ltm <- list.files(paste0(home, "/Landsat_ET/ltm_month"), full.names = T) # bring in ltm and ltsd to index for anom calcs
ltsd <- list.files(paste0(home, "/Landsat_ET/ltsd_month"), full.names = T)


UseCores <- 30
cl <- makeCluster(UseCores)
registerDoParallel(cl)

for(i in 1:length(unique_YM)){
  
  # subset to just the files corresponding to that date and grab the unique tiles - should be all but just in case
  sub <- file_dt[YM == unique_YM[i]]
  sub_unique_tiles <- unique(sub$tile)
  
  # loop through each tile for each date and calculate anomalies 
  ptm <- proc.time()
  foreach(j =  1:length(sub_unique_tiles))%dopar%{
   library(data.table)
    library(raster)
    subsub <- sub[tile == sub_unique_tiles[j]]$file
    
    stk <- do.call('stack', lapply(subsub, raster))
    if(nlayers(stk) == 1){
      mn_stk <- stk
    }else{mn_stk <- mean(stk, na.rm=T)}
    
    
    # bring in the long term mean and standard deviation 
    mo <- paste0(substr(unique_YM[i], 5, 6), ".tif") # grab the month we are looking at 
    
    ltm_mo <- ltm[grep(mo, ltm)]
    ltm_mo_tile <- ltm_mo[grep(sub_unique_tiles[j], ltm_mo)]
    
    ltsd_mo <- ltsd[grep(mo, ltsd)]
    ltsd_mo_tile <- ltsd_mo[grep(sub_unique_tiles[j], ltsd_mo)]
    
    ltm_month <- raster(ltm_mo_tile)
    ltsd_month <- raster(ltsd_mo_tile)
    
    
     anom <- (mn_stk - ltm_month)/ltsd_month
    
    writeRaster(mn_stk, paste0(home, "/Landsat_ET/monthlyMean/",sub_unique_tiles[j], "_" ,unique_YM[i], ".tif"), overwrite = T, format = "GTiff")
    writeRaster(anom, paste0(home, "/Landsat_ET/monthlyAnom/", sub_unique_tiles[j], "_" ,unique_YM[i], ".tif"), overwrite = T, format = "GTiff")

  }
  print(i)
  print(proc.time() - ptm)
  
}


