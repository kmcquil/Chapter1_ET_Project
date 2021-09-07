###################################################################################################### 
## MODIS: Get the average ET anomaly and ET residual during non-drought periods to use as a baseline 
######################################################################################################
library(rgdal)
library(raster)
library(lubridate)
library(data.table)
library(parallel)
library(Rmpi)
library(parallel)
library(snow)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

# start an rmpi cluster
workers <- mpi.universe.size() -1
cl <- makeMPIcluster(workers, type='MPI')

# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data/SPI/SPI90_MODIS"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data/SPI/SPI90_MODIS"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 

# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,1,4), substr(spi_short,5,6), "01"), "%Y%m%d"))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 4, 7), substr(et_short, 8, 9), "01"), "%Y%m%d"))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = 'date')
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

betas_file <- paste0(home, "/Analysis/outputs/MODIS/model_beta_coef.csv")
fmask <- paste0(home, "/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif")

avgNA_tl(files_dt_gs$spi_file, files_dt_gs$et_file, betas_file, fmask, cl, "MODIS", tl = NA)

stopCluster(cl)
mpi.quit()




