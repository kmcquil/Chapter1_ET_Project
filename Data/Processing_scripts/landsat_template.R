################################################################################
## Create landsat ET template to use to mosaic data onto a common grid 
## cropped to SBR ecoregion and crs = 32617
################################################################################
library(raster)
library(data.table)
library(rgdal)
library(gdalUtils)
library(stringr)
library(foreach)
library(doParallel)
library(lubridate)

UseCores <- detectCores() -1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "G:/My Drive/Chapter1_ET_Project/Data" # set home directory 

# first create a template that will be used to resample based on all landsat tiles used in the analysis 
et_long <- list.files(paste0(home, "/Landsat_ET/original_data"), full.names = T)
et_short <- list.files(paste0(home, "/Landsat_ET/original_data"), full.names = F)

# there are 12 tiles included 
tiles <- substr(et_short, 5, 10)
un_tiles <- unique(tiles)

files <- data.table(file = et_long, 
                    short_file = et_short, 
                    date = as.Date(substr(et_short, 11, 18), "%Y%m%d"), 
                    tile = substr(et_short, 5, 10))

# get the first record for each tile
record <- files[, .SD[1], tile]

# untar all of these to a template folder 
for(i in 1:nrow(record)){
  untar(record$file[i], 
        exdir = paste0(home, "/Landsat_ET/landsat_template"))
  
}

# grab all of the eta files from the template folder, mosaic them, and then mask and crop to the ROI 
temp_files <- list.files(paste0(home, "/Landsat_ET/landsat_template"), full.names = T)
temp_files <- temp_files[grep("eta", temp_files)]
crss <- lapply(temp_files, function(x) { crs(raster(x)) })

mosaic_rasters(gdalfile=temp_files, 
               dst_dataset=paste0(home, "/Landsat_ET/landsat_template/landsat_template.tif"), 
               t_srs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs", 
               of="GTiff", 
               overwrite = T)

gdalwarp(paste0(home, "/Landsat_ET/landsat_template/landsat_template.tif"), 
         paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"), 
         cutline = paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"), 
         crop_to_cutline = T, 
         overwrite = T)
