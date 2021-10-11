#############################################################################################
## Landsat tiles are too big once mosaiced over the full ROI so retile into 42 smaller tiles
#############################################################################################
library(gdalUtils)
library(raster)
library(sp)
library(rgdal)
library(doParallel)
library(foreach)
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project/Data"
#home <- "G:/My Drive/Chapter1_ET_Project/Data"

# initialize cluster 
UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)
print(cl)

# read in alexi shapefile that is a small 4km pixel around the flux tower
alexi <- readOGR(paste0(home, "/ALEXI/coweeta_alexi_pixel.shp"))

# list the landsat et files 
efile <- list.files(paste0(home, "/Landsat_ET/tifs_resampled"), full.names = T, pattern = ".tif$") # get the filepaths for et rasters
efile_short <- list.files(paste0(home, "/Landsat_ET/tifs_resampled"), full.names = F, pattern = '.tif$') # get the filepaths for et rasters

# convert crop shapefile to match crs of ET 
temp <- raster(efile[1])
alexi <- spTransform(alexi, crs(temp))
ex <- extent(alexi)

# crop the raster to shapefile 
foreach(i = 1:length(efile))%dopar%{
  library(raster)
  library(gdalUtils)
  gdal_translate(efile[i], 
                 paste0(home, "/ALEXI/Landsat_subset/", efile_short[i]), 
                 projwin=c(xmin(ex), ymax(ex), xmax(ex), ymin(ex)))
}

