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

UseCores <- 15
cl <- makeCluster(UseCores)
registerDoParallel(cl)

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

LTs <- function(TD){
  tile_dt <- file_dt[tile == TD]
  
  for(i in 1:length(unique_months)){
    sub <- tile_dt[month == unique_months[i]]$file
    print(length(sub))
    
    # create a data.frame from the rasters 
    sub_dt <- foreach(i = 1:length(sub), .combine = cbind)%dopar%{
      library(raster)
      as.data.frame(raster(sub[i]))
    }
    paste0('finished creating dataframe for month ', unique_months[i])
    
    # attach cell number so things don't get out of order 
    cellnum <- data.table(cellnum = seq(1, nrow(sub_dt)))
    sub_dt <- cbind(cellnum, sub_dt)
    print(dim(sub_dt))
    
    # drop rows that are totally empty
    NAcount <- data.table(NAcount = apply(sub_dt[,2:ncol(sub_dt)], 1, function(x){sum(is.na(x))}))
    Ncols <- ncol(sub_dt[,2:ncol(sub_dt)])
    sub_dt <- sub_dt[-NAcount$NAcount == Ncols]
    print(dim(sub_dt))
    

    if(dim(sub_dt)[1] == 0){next}
    
    # calculate the mean and sd for each cell 
    sub_dt$meanET <- as.numeric(rep(NA, nrow(sub_dt)))
    sub_dt$sdET <- as.numeric(rep(NA, nrow(sub_dt)))
    
    sel_cols <- sub_dt[, names(sub_dt)[seq(2, ncol(sub_dt)-2)]]
    
    subsub <- sub_dt[, ..sel_cols]
    c_meanET <- foreach(i = 1:nrow(subsub), .combine = c)%dopar%{
      library(data.table)
      y <- mean(as.numeric(subsub[i]), na.rm = T)
    }
    
    c_sdET <- foreach(i = 1:nrow(subsub), .combine = c)%dopar%{
      library(data.table)
      y <- sd(as.numeric(subsub[i]), na.rm = T)
    }
    
    sub_dt[,meanET := c_meanET]
    sub_dt[, sdET := c_sdET]
    
    
    # put back in raster format 
    # start with standard deviation
    cin <- c("sdET", "cellnum")
    sd_fin <- merge(cellnum, sub_dt[, ..cin], by = "cellnum", all.x=T)
    output_matrix <- matrix(sd_fin$sdET, nrow=nrow(ExampleRaster),ncol=ncol(ExampleRaster),byrow=T)
    new_output_raster<-raster(output_matrix,xmn=xmin(ExampleRaster),ymn=ymin(ExampleRaster),xmx=xmax(ExampleRaster),
                              ymx=ymax(ExampleRaster), crs=proj)
    writeRaster(new_output_raster, paste0(home, "/Landsat_ET/ltsd_month/", TD, "_",sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    
    # then do the mean
    cin <- c("meanET", "cellnum")
    mean_fin <- merge(cellnum, sub_dt[, ..cin], by = "cellnum", all.x=T)
    output_matrix <- matrix(mean_fin$meanET, nrow=nrow(ExampleRaster),ncol=ncol(ExampleRaster),byrow=T)
    new_output_raster<-raster(output_matrix,xmn=xmin(ExampleRaster),ymn=ymin(ExampleRaster),xmx=xmax(ExampleRaster),
                              ymx=ymax(ExampleRaster), crs=proj)
    writeRaster(new_output_raster, paste0(home, "/Landsat_ET/ltm_month/et_",TD, "_", sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    
    print(i)
    
  }
  
  
}




# apply the function to each 
lapply(unique_tiles, LTs)





TD <- unique_tiles[5]
tile_dt <- file_dt[tile == TD]
sub <- tile_dt[month == unique_months[4]]$file
rr <- lapply(sub, function(x){minValue(raster(x))})
unique(rr)


i = 1 : NA
i = 2 : NA
i = 3 : NA
i = 4 : NA
i = 5 : NA
i = 6 : NA
i = 7: 
  
  
  
  
clust <- makeCluster(15)
clusterExport(clust, c("sub"), envir=environment())
isNA <- function(x){library(raster); minValue(raster(x))}
a <- parLapply(clust, sub, isNA)
stopCluster(clust)







