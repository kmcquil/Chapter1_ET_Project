##############################################################################################
## Landsat drought peak calcs part 3 
##############################################################################################
library(rgdal)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n <- 32

#########################################################################################
# Calculate residuals for each drought peak using the linear model fit for each grid cell
# and then get the average residual for all years for each pixel 
# also get the corresponding ET anomalies and spi values 
#########################################################################################

# for each year, grab the ET anom at peak and SPI at peak 
# merge the files into one dt to make sure date of spi and et peak files match
peak_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/droughtPeak"), full.names = T, pattern = ".tif$")
etPeak_files <- data.table(et_file = peak_files[grep("droughtPeakETanom", peak_files)])
etPeak_files$date <- substr(etPeak_files$et_file, nchar(etPeak_files$et_file) - 20, nchar(etPeak_files$et_file) - 17)
etPeak_files$tile <- substr(etPeak_files$et_file, nchar(etPeak_files$et_file) - 15, nchar(etPeak_files$et_file) - 4)
spiPeak_files <- data.table(spi_file = peak_files[grep("droughtPeakSPI", peak_files)])
spiPeak_files$date <- substr(spiPeak_files$spi_file, nchar(spiPeak_files$spi_file)-20, nchar(spiPeak_files$spi_file) -17)
spiPeak_files$tile <- substr(spiPeak_files$spi_file, nchar(spiPeak_files$spi_file)-15, nchar(spiPeak_files$spi_file) -4)

peak_dt <- merge(etPeak_files, spiPeak_files, by = c("date", "tile"), all.x=T, all.y=T)

# find the years/tiles that did not have any drought and drop them from data table 
#minDrought <- unlist(lapply(peak_dt$spi_file, function(x) { minValue(raster(x))}))
minDrought <- unlist(lapply(peak_dt$spi_file, function(x) {cellStats(raster(x), min)}))
peak_dt <- peak_dt[!minDrought == Inf]

# beta coefficeint file 
beta_file <- list.files(paste0(home, "/Analysis/outputs/Landsat"), full.names = T, pattern = ".csv$")
beta_file <- beta_file[grep("model_beta_coef", beta_file)]
beta_file <- data.table(beta_file = beta_file, tile = substr(beta_file, nchar(beta_file)-15, nchar(beta_file)-4))

peak_dt <- merge(peak_dt, beta_file, by ="tile", all.x=T)
peak_dt <- peak_dt[complete.cases(peak_dt),]
dim(peak_dt)


UseCores <- n
cl <- makeCluster(UseCores)
print(cl)
registerDoParallel(cl)

foreach(i = 1:nrow(peak_dt))%dopar%{
  library(data.table)
  library(raster)
  spi_file <- peak_dt$spi_file[i]
  et_file <- peak_dt$et_file[i]
  beta_file1 <- peak_dt$beta_file[i]
  y <- peak_dt$date[i]
  tl <- peak_dt$tile[i]
  calcResiduals_tl(spi_file, et_file, beta_file1, y, tl, "Landsat")
}

stopCluster(cl)



