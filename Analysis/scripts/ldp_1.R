########################################################################################
## Drought analysis for Landsat data 
########################################################################################

library(rgdal)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n = 30 # how many cores do i want 

########################################################################################
# Identify droughts and retrieve ET anomalies at drought peak 
########################################################################################

# make sure to jsut delete spi files that are outside of the Landsat scope 
# modis data is from 2000-01 : 2019-12

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

# get a list of the years 
years <- unique(year(files_dt$date))
tiles <- unique(files_dt$tile)
paste0("The number of years ", length(years))

# only look in the growing season, so subset further from April (4) - September (9)
ptm <- proc.time()
for(i in 1:length(years)){
  y <- years[i]
  spi_sub <- files_dt[year(date) == y & month(date) >3 & month(date) < 10]
  et_sub <- files_dt[year(date) == y & month(date) >3 & month(date) < 10]
  print(y)
  for(j in 1:length(tiles)){
    tl <- tiles[j]
    spi_sub1 <- spi_sub[tile == tiles[j]]$spi_file
    et_sub1 <- et_sub[tile == tiles[j]]$et_file
    droughtID_tl(spi_sub1, et_sub1, n, home, y, tl,"Landsat")
    print(j)
  }
  
}

print("The drought peak rasters have been created!")
proc.time() - ptm
