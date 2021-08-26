############################################################################################
## resample forest composition to match landsat and modis 
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

wilson_fc <- "G:/My Drive/Dissertation/Forest Composition/FC_SBR/FC_classified_simple_majority.tif"

# MODIS 
## MODIS is ~500m resolution while the wilson maps are 250m resolution 
## resample the % diffuse, % ring, and % tracheid rasters to match MODIS 
## reclassify as diffuse, ring, tracheid by as imple majority BA in each cell 

modis_template <- raster(paste0(home, "/MODIS_ET/clean/et_20000101.tif")) #MODIS template to resample data to

# diffuse 
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/diffuse_percent.tif", 
      modis_template, 
      paste0(home, "/forest_composition/wilson_modis_diffuse_percent.tif"), 
      'bilinear')

# ring
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/ring_percent.tif", 
      modis_template, 
      paste0(home, "/forest_composition/wilson_modis_ring_percent.tif"), 
      'bilinear')

# tracheid
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/tracheid_percent.tif", 
      modis_template, 
      paste0(home, "/forest_composition/wilson_modis_tracheid_percent.tif"), 
      'bilinear')

diff <- raster(paste0(home, "/forest_composition/wilson_modis_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/wilson_modis_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/wilson_modis_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/wilson_modis_classified.tif"), format = "GTiff", overwrite = T)


# Now do the same for landsat 
landsat_template <- raster(paste0(home, "/Landsat_ET/landsat_template/landsat_template_ROI.tif"))

# diffuse 
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/diffuse_percent.tif", 
      landsat_template, 
      paste0(home, "/forest_composition/wilson_landsat_diffuse_percent.tif"), 
      'bilinear')

# ring
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/ring_percent.tif", 
      landsat_template, 
      paste0(home, "/forest_composition/wilson_landsat_ring_percent.tif"), 
      'bilinear')

# tracheid
warpMn("G:/My Drive/Dissertation/Forest Composition/FC_SBR/tracheid_percent.tif", 
      landsat_template, 
      paste0(home, "/forest_composition/wilson_landsat_tracheid_percent.tif"), 
      'bilinear')

diff <- raster(paste0(home, "/forest_composition/wilson_landsat_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/wilson_landsat_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/wilson_landsat_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/wilson_landsat_classified.tif"), format = "GTiff", overwrite = T)








#----------------------- Now do Riley maps ----------------------------------------------#

## MODIS first 
# diffuse 
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent.tif", 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_diffuse_percent.tif"), 
       'bilinear')

# ring
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent.tif", 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_ring_percent.tif"), 
       'bilinear')

# tracheid
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent.tif", 
       modis_template, 
       paste0(home, "/forest_composition/riley_modis_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/riley_modis_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/riley_modis_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/riley_modis_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack, maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/riley_modis_classified.tif"), format = "GTiff", overwrite = T)





## Now do Landsat 
# diffuse 
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent.tif", 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_diffuse_percent.tif"), 
       'bilinear')

# ring
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent.tif", 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_ring_percent.tif"), 
       'bilinear')

# tracheid
warpMn("G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent.tif", 
       landsat_template, 
       paste0(home, "/forest_composition/riley_landsat_tracheid_percent.tif"), 
       'bilinear')

diff <- raster(paste0(home, "/forest_composition/riley_landsat_diffuse_percent.tif"))
ring <- raster(paste0(home, "/forest_composition/riley_landsat_ring_percent.tif"))     
trach <- raster(paste0(home, "/forest_composition/riley_landsat_tracheid_percent.tif"))
class_stack <- stack(diff, ring, trach)
classified_maj <- calc(class_stack,  maj_fun)
writeRaster(classified_maj, paste0(home, "/forest_composition/riley_landsat_classified.tif"), format = "GTiff", overwrite = T)

