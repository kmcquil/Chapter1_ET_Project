#####################################################################################################
## MODIS calculate the 5 year rolling coupling between growing season monthly SPI and ET anomalies 
## Calculate the sen's slope of coupling on a pixel basis 
######################################################################################################
library(rgdal)
library(raster)
library(lubridate)
library(trend)
library(data.table)
library(parallel)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
#home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

UseCores <- detectCores()-1
cl <- makeCluster(UseCores)


# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data/SPI/SPI1Y_MODIS"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data/SPI/SPI1Y_MODIS"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 

# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,1,4), substr(spi_short,5,6), "01"), "%Y%m%d"))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 4, 7), substr(et_short, 8, 9), "01"), "%Y%m%d"))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = 'date')
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

fmask <- paste0(home, "/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif")

rollingCoupling(files_dt_gs$spi_file, files_dt_gs$et_file, files_dt_gs$date, fmask, cl, "MODIS_1Y", tl = NA)

stopCluster(cl)







