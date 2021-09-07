################################################################################################
## Convert the results of ldp_3.R (residuals, ETanoms, and spi) for each year to rasters and then to DF
## This way each year is just one file instead of 42 
################################################################################################
library(rgdal)
library(gdalUtils)
library(raster)
library(lubridate)
library(data.table)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

# landsat template to retile all of the separate tiles onto 
template <- paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif")
# list of the all tiles 
tiles <- list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), pattern = ".tif$", full.names=T)

# retile to one raster and convert to a DF for the SPI files at drought peak 
spi_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/spi"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/spi/full/")
retiling(spi_files, tiles, template, output_dir, 10)

# do the same for residuals 
res_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/residuals"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/residuals/full/")
retiling(res_files, tiles, template, output_dir, 10)


# do the same for ET anomalies 
et_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/ETanom"), full.names=T, pattern = ".csv$")
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/ETanom/full/")
retiling(et_files, tiles, template, output_dir, 10)

