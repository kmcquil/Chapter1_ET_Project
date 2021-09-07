################################################################################
################################################################################
## Process MODIS ET so it is ready to use for analysis 
### .hdf to .tif 
## 
################################################################################
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

home <- "G:/My Drive/Chapter1_ET_Project/Data" # set a home directory

# I originally downloaded all of the .hdf files for apps from maine to alabama. 
# Now I only want to .hdf files converted to .tif files that intersect the SBR ecoregion 
# So only keep v5h11 and v5h10 
roi <- readOGR(paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp")) # bring in the ROI 
hdf <- list.files(paste0(home, "/MODIS_ET/HDF"), full.names = T, pattern = ".hdf$") # get list of all.hdf files 
hdf_sub <- hdf[grep("h11v05|h10v05", hdf)]

hdf_short <- list.files(paste0(home, "/MODIS_ET/HDF"), full.names = F, pattern = ".hdf$") # get list of all.hdf files 
hdf_short_sub <- hdf_short[grep("h11v05|h10v05", hdf_short)]

gtiffs <- gsub("hdf", "tif", hdf_sub) # replace the .hdf with .tif for saving new files 
gtiffshort <- gsub("hdf", "tif", hdf_short_sub)

################################################################################################################
# Convert .hdf to .tif
# convert all no-data values to NA
# looking at the data, the scale has already been done --- i need to find my google doc that looked at where i downloaded this from
lyr <- 1 # this corresponds to the ET 
foreach(i = 1:length(hdf_sub))%dopar%{
  library(raster)
  library(gdalUtils)
  library(rgdal)
  
  # convert the hdf to .tif files 
  gdal_translate(hdf_sub[i],
                 paste0(gtiffs[i]),
                 sd_index=lyr, overwrite = TRUE) # extract the specified HDF layer and save it as a Geotiff
  
  # read in the raster and mask out NoData values 
  r <- raster(gtiffs[i])
  r[r > 3276] <- NA
  writeRaster(r, gtiffs[i], overwrite = T, format="GTiff")
  
}

####################################################################################################
# mosaic by date
file_dt <- data.table(file = gtiffs)
file_dt$date <- as.Date(substr(file_dt$file, 62, 68), "%Y%j")
dates <- unique(file_dt$date)

foreach(i = 1:length(dates))%dopar%{
  library(data.table)
  library(gdalUtils)
  library(rgdal)
  library(raster)
  library(lubridate)
  files <- file_dt[date == dates[i]]$file
  mosaic_rasters(gdalfile=files,
                 dst_dataset = paste0(home, "/MODIS_ET/mosaiced/", year(dates[i]), sprintf("%02d",month(dates[i])), sprintf("%02d", day(dates[i])), '.tif'),
                 of="GTiff",
                 overwrite = TRUE)
}

######################################################################################################
# reporject to 32617, crop, and mask to ROI
frmSRS = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" # original HDF SRS
toSRS = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" # desired GeoTIFF SRS
mos_files <- list.files(paste0(home, "/MODIS_ET/mosaiced"), full.names=T)
mos_files_short <- list.files(paste0(home, "/MODIS_ET/mosaiced"), full.names=F)
foreach(i = 1:length(mos_files))%dopar%{
  library(raster)
  library(gdalUtils)
  library(rgdal)
  gdalwarp(mos_files[i], 
           paste0(home, "/MODIS_ET/clean/et_", mos_files_short[i]), 
           s_srs=frmSRS, 
           t_srs=toSRS, 
           overwrite = T, 
           cutline = paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"), 
           crop_to_cutline = T
           )
}



##############################################################################################################
## there was only one file for "2016-03-29" instead of two
## Do that file separate by resampling it onto the grid of a date that was mosaiced already 
xx <- list.files("G:/My Drive/Chapter1_ET_Project/Data/MODIS_ET/HDF", pattern = ".hdf$", full.names = T)
dtt <- data.table(hdf_file = hdf_sub, date = as.Date(substr(hdf_short_sub, 12, 18), "%Y%j"))
dtt[date == "2016-03-29"]

fin <- "G:/My Drive/Chapter1_ET_Project/Data/MODIS_ET/HDF/MOD16A2GF.A2016089.h11v05.006.2020021143244.tif"
fout <- "G:/My Drive/Chapter1_ET_Project/Data/MODIS_ET/mosaiced/20160329.tif"
template <- raster("G:/My Drive/Chapter1_ET_Project/Data/MODIS_ET/mosaiced/20000101.tif")

res <- res(template)
t1 <- c(xmin(template), ymin(template), 
        xmax(template), ymax(template))
res_out <- crs(template)
re <- 'near'
gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r=re, overwrite = T)

frmSRS = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" # original HDF SRS
toSRS = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs" # desired GeoTIFF SRS
gdalwarp(fout, 
         paste0(home, "/MODIS_ET/clean/et_", "20160329.tif"), 
         s_srs=frmSRS, 
         t_srs=toSRS, 
         overwrite = T, 
         cutline = paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"), 
         crop_to_cutline = T
)