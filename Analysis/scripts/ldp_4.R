##############################################################################################
## Landsat drought peak calcs part 4
##############################################################################################
library(rgdal)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n = 10

## create data tables of each variable (spi, et, and et residuals) 
# residual files 
resid_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/residuals/full"), full.names = T, pattern = ".csv$")
avgResid(resid_files, n, "Landsat")

# now do ET 
yrET_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/ETanom/full"), full.names = T, pattern = ".csv$")
avgET(yrET_files, n, "Landsat")

# now SPI 
yrSPI_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/spi/full"), full.names = T, pattern = ".csv$")
avgSPI(yrSPI_files, n, "Landsat")


print("The average drought residuals have been calcualted!")


