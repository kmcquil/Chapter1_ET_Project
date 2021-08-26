################################################################################################
## Convert the results of ldp_3.R (residduals, ETanoms, and spi) for each year to rasters and then to DF
## This way each year is just one file instead of 42 
################################################################################################
library(rgdal)
library(gdalUtils)
library(raster)
library(lubridate)
library(data.table)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
#home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

template <- paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif")
tiles <- list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), pattern = ".tif$", full.names=T)

# for the spi files 
spi_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/spi"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/spi/full/")
retiling(spi_files, tiles, template, output_dir, 10)

# for the residual files 
res_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/residuals"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/residuals/full/")
retiling(res_files, tiles, template, output_dir, 10)


# for the ET anomaly files 
et_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/ETanom"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/ETanom/full/")
retiling(et_files, tiles, template, output_dir, 10)

