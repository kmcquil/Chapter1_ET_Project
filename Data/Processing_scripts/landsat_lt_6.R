####################################################################################################
## Calculate monthly landsat ET long term mean and standard deviation 
####################################################################################################
tileID <- seq(26, 30)
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
## Add in column for tile 
file_dt$tile <- substr(efile_short, 1, 12)

unique_months <- unique(file_dt$month) # for the long term calcs, we will calculate for each month 
unique_tiles <- unique(file_dt$tile)

beginCluster(30)

untls <- unique_tiles[tileID]
apply(untls , LTs)

endCluster()