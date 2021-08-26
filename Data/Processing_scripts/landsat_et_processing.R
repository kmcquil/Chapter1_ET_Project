################################################################################
################################################################################
## Process landsat ET 
library(raster)
library(data.table)
library(rgdal)
library(gdalUtils)
library(stringr)
library(foreach)
library(doParallel)
library(lubridate)

# UseCores <- detectCores() -1
# cl <- makeCluster(UseCores)
# registerDoParallel(cl)
# print(cl)

# base directory
home <- "G:/My Drive/Chapter1_ET_Project/Data"
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"

## All data should be cropped to the SBR ecoregion and projected to EPSG:32617
################################################################################
################################################################################
roi_mask <- readOGR( paste0(home, "/NA_CEC_Eco_Level3/blue_ridge.shp"))


# make a list of the tar.fz files 
# files <- list.files(paste0(home, "/Landsat_ET/original_data"), full.names = T, pattern = ".tar.gz")
# files_short <- list.files(paste0(home, "/Landsat_ET/original_data"), full.names = F, pattern = ".tar.gz")
# 
# # untar the files and output .tif files for ETa and ETqa
# foreach(i = 1:length(files))%dopar%{
#   library(raster)
#   library(rgdal)
#   library(gdalUtils)
#   library(lubridate)
#   library(data.table)
# 
#   # get the names of files want from the tar.gz
#   ls <- untar(files[i], list = T)
#   ls <- ls[grep("eta|etqa", ls)]
#   # untar the file
#   untar(tarfile = files[i],
#         files = ls,
#         exdir = paste0(home, "/Landsat_ET/original_data/untar"))
# }
# 
# print('finished untarring')
# 
# # scale and QC each ETa file 
# # list the eta files 
# eta_files <- list.files(paste0(home, "/Landsat_ET/original_data/untar"), full.names = T, pattern = ".tif$")
# eta_files <- eta_files[grep("eta", eta_files)]
# eta_files_short <- list.files(paste0(home, "/Landsat_ET/original_data/untar"), full.names = F, pattern = ".tif$")
# eta_files_short <- eta_files_short[grep("eta", eta_files_short)]
# 
# # list the qc files 
# qc_files <- list.files(paste0(home, "/Landsat_ET/original_data/untar"), full.names = T, pattern = ".tif$")
# qc_files <- qc_files[grep("etqa", qc_files)]
# 
# checking <- data.table(filename = as.character(rep(NA, 6961)))
# 
# foreach(i = 1:length(eta_files))%dopar%{
#   library(raster)
#   library(rgdal)
#   library(gdalUtils)
#   library(lubridate)
#   library(data.table)
#   
#   eta <- raster(eta_files[i])
#   eta <- eta*0.001
#   
#   # grab the corresponding qc file based on tile id and date
#   id <- substr(eta_files_short[i], 11, 25)
#   qcf <- raster(qc_files[grep(id, qc_files)])
#   qcf[qcf >1] <- NA
#   
#   eta <- mask(eta, qcf)
#   writeRaster(eta, paste0(home, "/Landsat_ET/tifs_screened/", eta_files_short[i]), overwrite = T, format="GTiff")
# }
# 
# 
# stopCluster(cl)
# 
# print('finished scale and qc')
# 
# start a new cluster that uses less cores 
UseCores <- 20
cl <- makeCluster(UseCores)
registerDoParallel(cl)


# mosaic tiles from the same day and resample to the common grid 
et_files <- list.files(paste0(home, "/Landsat_ET/tifs_screened"), full.names=T)
et_files_short <- list.files(paste0(home, "/Landsat_ET/tifs_screened"), full.names=F)
files <- data.table(file = et_files, 
                    file_short = et_files_short, 
                    date = as.Date(substr(et_files_short, 18, 25), "%Y%m%d"))

dates <- unique(files$date)

template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
res <- res(template)
res_out <- crs(template)
t1 <- c(xmin(template), ymin(template), 
        xmax(template), ymax(template))

roi_mask <- spTransform(roi_mask, CRS(" +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84
+ +towgs84=0,0,0"))
                        

foreach(i = 1:length(dates))%dopar%{
  library(raster)
  library(rgdal)
  library(gdalUtils)
  library(lubridate)
  library(data.table)
  
  sub <- files[date == dates[i],]
  fout <- paste0(home, "/Landsat_ET/tifs_mosaiced/et_", year(dates[i]), sprintf("%02d", month(dates[i])), sprintf("%02d", day(dates[i])), ".tif")
  mosaic_rasters(gdalfile = sub$file, 
                 dst_dataset = fout, 
                 of="GTiff",
                 overwrite = TRUE)
  
  foutout <- paste0(home, "/Landsat_ET/tifs_resampled/et_", year(dates[i]), sprintf("%02d", month(dates[i])), sprintf("%02d", day(dates[i])), ".tif")
  gdalwarp(fout, foutout, t_srs=res_out, tr=res, te = t1, r='bilinear', overwrite = T, dstnodata = NA)
  r <- raster(foutout)
  r <- mask(r, roi_mask)
  writeRaster(r, foutout, overwrite = T, type = "GTiff")
  
  
}


