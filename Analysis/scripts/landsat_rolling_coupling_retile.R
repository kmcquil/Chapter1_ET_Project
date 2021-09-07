##############################################################################################################
## Retile the overall coupling and trends in rolling coupling 
## Put the results one raster/column instead of multiple
##############################################################################################################
library(rgdal)
library(gdalUtils)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)

# start a cluster 
UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project" # set a home directory 
source(paste0(home, "/Analysis/scripts/analysis_funcs.R")) # bring in the analysis functions 

# templates for individual tiles and the full roi
template <- paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif")
tile_template <- data.table(tile_file = list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), pattern = ".tif$", full.names=T))
tile_template$tile <- substr(tile_template$tile_file, nchar(tile_template$tile_file)-15, nchar(tile_template$tile_file)-4)

# list all of the .csv files of sens slope and overall coupling results for each tile 
files <- list.files(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling"), full.names=T, pattern = ".csv$")
rolling_files <- data.table(file = files[grep("Coupling_x", files)])
rolling_files$tile <- substr(rolling_files$file, nchar(rolling_files$file) - 15, nchar(rolling_files$file)-4)
sens_files <- data.table(file = files[grep("sens", files)])
sens_files$tile <- substr(sens_files$file, nchar(sens_files$file)-15, nchar(sens_files$file)-4)

# join the tile templates with the results. Each result file is associated with a tile. 
rolling_files <- merge(rolling_files, tile_template, by = "tile", all.x=T) # coupling 
sens_files <- merge(sens_files, tile_template, by = "tile", all.x=T) # trend in coupling 

# Retile the trend in coupling files (sens_files) 
# first convert each file to a raster using the template associated with each file 
# drop rows that do not have a significant result to return a raster with only significant trends in coucpling 
foreach(i = 1:nrow(sens_files))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(sens_files$file[i])
  dtIN <- dtIN[(!is.na(pvalue)==T) & pvalue<=0.05,]
  if(nrow(dtIN) == 0)return(NULL) # if this result had no signifiacnt trends just move onto the next file 
  
  sel_cols <- c("cellnum","senSlope")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(sens_files$tile_file[i])
  fout <- gsub(".csv", ".tif", sens_files$file[i])
  DTtoRast(dtIN, ExampleRaster, fout, open=F) # convert to a raster 
}

# take all of rasters with significant slopes and retile into one final raster 
tif_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling"), full.names=T, pattern = ".tif$")
rast_list <- gsub(".csv", ".tif", sens_files$file)
rast_list <- rast_list[rast_list %in% tif_files]
out <- paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/sig_sens_coupling_final.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)


# Now do overall coupling 
# Convert overall coupling .csv files to rasters using the tile templates and then retile to one raster 
foreach(i = 1:nrow(rolling_files))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(rolling_files$file[i])
  # these I only kept overall coupling that was significant so no need to filter
  sel_cols <- c("cellnum","all")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(rolling_files$tile_file[i])
  fout <- gsub(".csv", ".tif", rolling_files$file[i])
  DTtoRast(dtIN, ExampleRaster, fout, open=F)
  
}
rast_list <- gsub(".csv", ".tif", rolling_files$file)
out <- paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/sig_overall_coupling_final.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)


