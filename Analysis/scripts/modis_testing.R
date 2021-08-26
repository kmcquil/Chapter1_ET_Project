########################################################################################
## Drought analysis test on a small dataset 
########################################################################################

library(rgdal)
library(raster)
library(lubridate)
library(data.table)

home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n = 6 # how many cores do i want 

########################################################################################
# Identify droughts and retrieve ET anomalies at drought peak 
########################################################################################

# make sure to jsut delete spi files that are outside of the modis scope 
# modis data is from 2000-01 : 2019-12

# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data_Testing/SPI/SPI90_MODIS"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data_Testing/SPI/SPI90_MODIS"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data_Testing/MODIS_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data_Testing/MODIS_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 


# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,1,4), substr(spi_short,5,6), "01"), "%Y%m%d"))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 4, 7), substr(et_short, 8, 9), "01"), "%Y%m%d"))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = 'date')

# get a list of the years 
years <- unique(year(files_dt$date))

# only look in the growing season, so subset further from April (4) - September (9)
for(i in 1:length(years)){
  y <- years[i]
  spi_sub <- files_dt[year(date) == y & month(date) >3 & month(date) < 10]$spi_file
  et_sub <- files_dt[year(date) == y & month(date) >3 & month(date) < 10]$et_file
  
  droughtID(spi_sub, et_sub, n, home, y, "MODIS")
}



########################################################################################
# Fit the linear model between ETanom ~ SPI for each grid cell
########################################################################################

# subset the spi and et files to just growing season 
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

# grab the forest mask file path for the modis data 
fmask <- paste0(home, "/Data_Testing/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif")

# calculate the linear models and store B0 and B1 coefficients in DT for future use 
modFit(files_dt$spi_file, files_dt$et_file, fmask, n, "MODIS")

#rollingCoupling(files_dt_gs$spi_file, files_dt_gs$et_file, files_dt_gs$date, fmask, cl, "MODIS", tl = NA)


#########################################################################################
# Calculate residuals for each drought peak using the linear model fit for each grid cell
# and then get the average residual for all years for each pixel 
#########################################################################################

# for each year, grab the ET anom at peak and SPI at peak 
# merge the files into one dt to make sure date of spi and et peak files match
peak_files <- list.files(paste0(home, "/Analysis/outputs/MODIS/droughtPeak"), full.names = T, pattern = ".tif$")
etPeak_files <- data.table(et_file = peak_files[grep("droughtPeakETanom", peak_files)])
etPeak_files$date <- substr(etPeak_files$et_file, nchar(etPeak_files$et_file) - 7, nchar(etPeak_files$et_file) - 4)
spiPeak_files <- data.table(spi_file = peak_files[grep("droughtPeakSPI", peak_files)])
spiPeak_files$date <- substr(spiPeak_files$spi_file, nchar(spiPeak_files$spi_file)-7, nchar(spiPeak_files$spi_file) -4)

peak_dt <- merge(etPeak_files, spiPeak_files, by = "date")

# find the years that did not have any drought and drop them from data table 
minDrought <- unlist(lapply(peak_dt$spi_file, function(x) { minValue(raster(x))}))
peak_dt <- peak_dt[!is.na(minDrought) == T]

# beta coefficeint file 
beta_file <- paste0(home, "/Analysis/outputs/MODIS/model_beta_coef.csv")

UseCores <- n
cl <- makeCluster(UseCores)
registerDoParallel(cl)

foreach(i = 1:nrow(peak_dt))%dopar%{
  library(data.table)
  library(raster)
  spi_file <- peak_dt$spi_file[i]
  et_file <- peak_dt$et_file[i]
  y <- peak_dt$date[i]
  calcResiduals(spi_file, et_file, beta_file, y, "MODIS")
}

stopCluster(cl)

resid_files <- list.files(paste0(home, "/Analysis/outputs/MODIS/residuals"), full.names = T, pattern = ".csv$")
avgResid(resid_files, n, "MODIS")




