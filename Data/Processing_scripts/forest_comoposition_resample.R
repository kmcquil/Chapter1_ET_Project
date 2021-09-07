############################################################################################
## resample forest composition % diffuse, ring, tracheid to match landsat and modis 
## calculate the majority type of BA in each cell using the resampled diff, ring, and tracheid and save 
## resample ty wilson maps to match MODIS and Landsat
## resample riley maps to match MODIS and Landsat 
############################################################################################
library(raster)
library(rgdal)
library(gdalUtils)
home <- "G:/My Drive/Chapter1_ET_Project/Data" # base directory
source(paste0(home, "/Processing_scripts/helpers.R")) # source a helper function 

# function to find index of majority 
maj_fun <- function(x){
   if(sum(is.na(x[1:3])) == 3){
      y <- NA
   }else{
      y <- which.max(x[1:3])
   }
}

#----------------------- Start with Ty Wilson maps ----------------------------------------------#
# Resample to MODIS resolution 
## MODIS is ~500m resolution while the wilson maps are 250m resolution 
## resample the % diffuse raster to match MODIS 
modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) #MODIS template to resample data to
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/diffuse_percent.tif"), 
      modis_template, 
      paste0(home, "/forest_composition/wilson_modis_diffuse_percent.tif"), 
      'bilinear')

## resample the % ring raster to match MODIS 
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/ring_percent.tif"), 
       modis_template, 
       paste0(home, "/forest_composition/wilson_modis_ring_percent.tif"), 
       'bilinear')

## resample the % tracheid raster to match MODIS 
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/tracheid_percent.tif"), 
       modis_template, 
       paste0(home, "/forest_composition/wilson_modis_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/wilson_modis_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/wilson_modis_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/wilson_modis_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/wilson_modis_classified.tif"), format = "GTiff", overwrite = T)


# Resample for Landsat resolution
# resample the % diffuse raster to match Lansdat 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/diffuse_percent.tif"), 
      landsat_template, 
      paste0(home, "/forest_composition/wilson_landsat_diffuse_percent.tif"), 
      'bilinear')

# resample the % ring raster to match Lansdat 
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/ring_percent.tif"), 
       landsat_template, 
       paste0(home, "/forest_composition/wilson_landsat_ring_percent.tif"), 
       'bilinear')

# resample the % tracheid raster to match Lansdat 
warpMn(paste0(home, "/forest_composition/Wilson_FC_SBR/tracheid_percent.tif"), 
       landsat_template, 
       paste0(home, "/forest_composition/wilson_landsat_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/wilson_landsat_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/wilson_landsat_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/wilson_landsat_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/wilson_landsat_classified.tif"), format = "GTiff", overwrite = T)


#----------------------- Now do Riley maps ----------------------------------------------------#

## Resample to MODIS resolution
## Resample the Riley % diffuse to the MODIS resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/diffuse_percent.tif"), 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_diffuse_percent.tif"), 
       'bilinear')

## Resample the Riley % ring to the MODIS resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/ring_percent.tif"), 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_ring_percent.tif"), 
       'bilinear')

## Resample the Riley % tracheid to the MODIS resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/tracheid_percent.tif"), 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/riley_modis_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/riley_modis_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/riley_modis_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/riley_modis_classified.tif"), format = "GTiff", overwrite = T)


## Resample to Landsat resolution 
## Resample the Riley % diffuse to the Landsat resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/diffuse_percent.tif"), 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_diffuse_percent.tif"), 
       'bilinear')

## Resample the Riley % ring to the Landsat resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/ring_percent.tif"), 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_ring_percent.tif"), 
       'bilinear')

## Resample the Riley % tracheid to the Landsat resolution 
warpMn(paste0(home, "/forest_composition/Riley_FC_SBR/tracheid_percent.tif"), 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/riley_landsat_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/riley_landsat_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/riley_landsat_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack,  maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/riley_landsat_classified.tif"), format = "GTiff", overwrite = T)


