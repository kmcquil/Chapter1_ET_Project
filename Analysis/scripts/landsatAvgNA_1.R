###########################################################################################
## Find the average ET anomaly and ET anomaly residual during normal wetness conditions 
## Normal means SPI < 0.5 and SPI > -0.5 
###########################################################################################
library(rgdal)
library(raster)
library(lubridate)
library(data.table)
library(Rmpi)
library(parallel)
library(snow)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

########################################################################################
# Identify droughts and retrieve ET anomalies at drought peak 
########################################################################################
# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data/SPI/retile"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data/SPI/retile"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data/Landsat_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data/Landsat_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 

# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,14,17), substr(spi_short,18,19), "01"), "%Y%m%d"), 
                  tile = substr(spi_short, 1, 12))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 14, 17), substr(et_short, 18, 19), "01"), "%Y%m%d"), 
                 tile = substr(et_short, 1, 12))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = c('date', 'tile'))

# get a list of the years and tiles 
years <- unique(year(files_dt$date))
tiles <- unique(files_dt$tile)

# subset the spi and et files to just growing season 
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

# grab the forest mask file path 
fmask <- list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), full.names = T, pattern = ".tif")

# grab the beta coefficient data table file paths which we will use to calculate residuals 
beta_file <- list.files(paste0(home, "/Analysis/outputs/Landsat"), full.names = T, pattern = ".csv$")
beta_file <- beta_file[grep("model_beta_coef", beta_file)]

# start an rmpi cluster  
workers <- mpi.universe.size() -1
cl <- makeMPIcluster(workers, type='MPI')

# for the first 14 tiles, calcualte the average no drought ET anomaly and residual by looping through each tile
# save the results in a .csv 
for(i in 1:14){
  print("start")
  ptm <- proc.time()
  
  fmask_sub <- fmask[grep(tiles[i], fmask)]
  beta_sub <- beta_file[grep(tiles[i], beta_file)]
  if(length(beta_sub) == 0){next} # if there is no beta file then the tile did not have any data and we can skip to the next one
  print(beta_sub)
  files_dt_gs_tl <- files_dt_gs[tile == tiles[i]]
  avgNA_tl(files_dt_gs_tl$spi_file, files_dt_gs_tl$et_file, beta_sub, fmask_sub, cl, "Landsat", tl = tiles[i])
  
  print(proc.time() - ptm)
  print(i)
}

stopCluster(cl)
mpi.quit()

