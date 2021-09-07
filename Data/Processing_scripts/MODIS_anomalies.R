################################################################################################
## Calculate MODIS ET long term mean and SD of average monthly ET rate
## average monthly ET rate and z-score seaasonally standardized anomalies  
################################################################################################
library(raster)
library(data.table)
library(rgdal)
library(gdalUtils)
library(stringr)
library(foreach)
library(doParallel)
library(lubridate)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"

# create a DT of files with the date to easily subset files for anomaly calculations 
files <- list.files(paste0(home, "/MODIS_ET/clean"), full.names = T, pattern = ".tif$")
files_short <- list.files(paste0(home, "/MODIS_ET/clean"), full.names = F, pattern = ".tif$")

file_dt <- data.table(file = files, 
                      file_short = files_short, 
                      date = as.Date(substr(files_short, 4, 11), "%Y%m%d")) # get the date of files in a data.table

file_dt$month <- month(file_dt$date) # find month and year to help with subsetting 
file_dt$year <- year(file_dt$date)
file_dt$YM <- paste0(file_dt$year, sprintf("%02d", file_dt$month))

#################################################################################################
# Calculate the long term mean and standard deviation of ET for every month 
unique_months <- unique(file_dt$month) 
beginCluster(n = 20) # start a cluster 
for(i in 1:length(unique_months)){
  sub <- file_dt[month == unique_months[i]]$file
  stk <- do.call('stack', lapply(sub, raster))
  avg_stk <- clusterR(stk, fun = calc, args =list(fun = mean, na.rm = T))
  sd_stk <- clusterR(stk, fun = calc, args = list(fun=sd, na.rm=T))
  
  writeRaster(avg_stk, paste0(home, "/MODIS_ET/ltm_month/et_", sprintf("%02d", unique_months[i]), ".tif"), 
              overwrite = T, format = "GTiff")
  writeRaster(sd_stk, paste0(home, "/MODIS_ET/ltsd_month/et_", sprintf("%02d", unique_months[i]), ".tif"), 
              overwrite = T, format = "GTiff")
  print(i)
}
print('longn term calculations are finished')

##############################################################################################
# then calculate average ET and anomaly for every month-year 
unique_YM <- unique(file_dt$YM) # get a list of unique year-month combos

# bring in ltm and ltsd to index for anom calcs
ltm <- list.files(paste0(home, "/MODIS_ET/ltm_month"), full.names = T) 
ltsd <- list.files(paste0(home, "/MODIS_ET/ltsd_month"), full.names = T)

for(i in 1:length(unique_YM)){
 
  # calc the mean of that year and month 
  sub <- file_dt[YM == unique_YM[i]]$file
  stk <- do.call('stack', lapply(sub, raster))
  mn_stk <- clusterR(stk, fun = calc, args = list(fun = mean, na.rm=T))
  
  # bring in the long term mean and standard deviation 
  mo <- substr(unique_YM[i], 5, 6) # grab the month we are looking at 
  ltm_month <- raster(ltm[grep(mo, ltm)])
  ltsd_month <- raster(ltsd[grep(mo, ltsd)])
  
  # calc z-score anoms 
  suub_fun <- function(x) { yy = ((x[1] - x[2])/x[3]); return(yy)}
  a_stk <- stack(mn_stk, ltm_month, ltsd_month)
  anom <- clusterR(a_stk, fun = calc, args = list(fun = suub_fun))
  
  writeRaster(mn_stk, paste0(home, "/MODIS_ET/monthlyMean/et_", unique_YM[i], ".tif"), overwrite = T, format = "GTiff")
  writeRaster(anom, paste0(home, "/MODIS_ET/monthlyAnom/et_", unique_YM[i], ".tif"), overwrite = T, format = "GTiff")
  print(i)
}
endCluster()



