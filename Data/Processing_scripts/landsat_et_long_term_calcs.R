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

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
source(paste0(home, "/Processing_scripts/helpers.R"))

## create data.table to easily subset all of the files by date and tile
efile <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = T) # get the filepaths for et rasters
efile_short <- list.files(paste0(home, "/Landsat_ET/retile"), full.names = F) # get the raster names 
file_dt <- data.table(file = efile, 
                      file_short = efile_short, 
                      date = as.Date(substr(efile_short, 17, 24), "%Y%m%d")) # get the date of files in a data.table
file_dt$month <- month(file_dt$date) # find month and year to help with subsetting 
file_dt$year <- year(file_dt$date)
file_dt$YM <- paste0(file_dt$year, sprintf("%02d", file_dt$month))
file_dt$tile <- substr(efile_short, 1, 12)

# list the unique months and tiles 
unique_months <- unique(file_dt$month) 
unique_tiles <- unique(file_dt$tile)

# start a cluster 
beginCluster(32)
# calculate the long term mean andn long term standard deviation of each month for each tile 
apply(unique_tiles, LTs)
endCluster()

