##############################################################################################
## Retile the Landsat average no drought ET anomaly and residual .csv files 
##############################################################################################
library(rgdal)
library(gdalUtils)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)

# start cluster
n <- 10 
UseCores <- n
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project" 
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

# get the templates for whole ROI and for individual tiles 
template <- paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif")
tile_template <- list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), pattern = ".tif$", full.names=T)

# get the no drought ET files 
inFiles <- list.files(paste0(home, "/Analysis/outputs/Landsat"), full.names=T, pattern = ".csv$")
inFiles <- na_files[grep("avgNonDroughtETResponse", na_files)]
output_dir <- paste0(home,  "/Analysis/outputs/Landsat/")

# put the files in a DT and create a tile column to subset easily 
dt <- data.table(file = inFiles, 
                 tile = substr(inFiles, nchar(inFiles)-15, nchar(inFiles)-4))
# create another column that subs .csv to .tif 
dt$tif_file_anom <- gsub(".csv", "_anom.tif", dt$file)
dt$tif_file_res <- gsub(".csv", "_res.tif", dt$file)

# convert tile filenames to DT and create column for tile ID 
tile_template <- data.table(tile_file = tile_template, 
                            tile = substr(tile_template, nchar(tile_template) - 15, nchar(tile_template)-4))
# merge no drought DT and tile DT to associate tile template to each no drought file 
dt <- merge(dt, tile_template, by = "tile", all.x=T)

### Convert each file to a raster with the same name 
### Regrid all of the new rasters to the overall template 
### Convert the new full raster to a data.table with columns cellnum and etAnom and etRes
# first do anomalies 
foreach(i = 1:nrow(dt))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(dt$file[i])
  sel_cols <- c("cellnum","etAnomAvg")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(dt$tile_file[i])
  fout <- dt$tif_file_anom[i]
  DTtoRast(dtIN, ExampleRaster, fout, open=F)
}
# put the anoms into one big raster 
rast_list <- dt$tif_file_anom
out <- paste0(output_dir, "avgNonDroughtETResponse_anoms.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)


# next to residuals 
foreach(i = 1:nrow(dt))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(dt$file[i])
  sel_cols <- c("cellnum","etResAvg")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(dt$tile_file[i])
  fout <- dt$tif_file_res[i]
  DTtoRast(dtIN, ExampleRaster, fout, open=F)
}

# put the residuals into one big raster 
rast_list <- dt$tif_file_res
out <- paste0(output_dir, "avgNonDroughtETResponse_res.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)


# put them both back into one csv 
anom_dt <- as.data.frame(raster(paste0(output_dir, "avgNonDroughtETResponse_anoms.tif")))
res_dt <- as.data.frame(raster(paste0(output_dir, "avgNonDroughtETResponse_res.tif")))
cellnum <- as.data.frame(seq(1, nrow(res_dt)))
all_dt <- cbind(cellnum, anom_dt, res_dt)
colnames(all_dt) <- c("cellnum", "etAnomAvg", "etResAvg")
fwrite(all_dt, paste0(output_dir, "avgNonDroughtETResponse.csv"))
