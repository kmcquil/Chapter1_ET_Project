####################################################################################3
## Resample the TWI computed at 30m in ArcGIS Pro to match Landsat data 
## Aggregate to 500m to match MODIS data 
library(rgdal)
library(raster)

# Bring in sample MODIS raster 
home <- "G:/My Drive/Chapter1_ET_Project/Data/Topography/"
modis_elevation <- raster(paste0(home, "usgsNED_elevation/elevation500m.tif"))
 
# Bring in sample Landsat raster  
landsat_elevation <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/elevation30m.tif")

# Bring in boundary to crop to 
sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")

# Bring in the Lansat TWI and mask and crop it to the boundary and resample to match other data 
landsat_twi <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/TWI/TWI.tif")
landsat_twi <- mask(crop(landsat_twi, sbr), sbr)
landsat_twi_resamp <- resample(landsat_twi, landsat_elevation)
writeRaster(landsat_twi_resamp, "G:/My Drive/Chapter1_ET_Project/Data/Topography/TWI/TWI_landsat.tif", 
            format="GTiff", overwrite = T)


# resample the landsat TWI to MODIS resolution 
modis_twi <- resample(landsat_twi, modis_elevation)
writeRaster(modis_twi, "G:/My Drive/Chapter1_ET_Project/Data/Topography/TWI/TWI_modis.tif",
            format = "GTiff", overwrite = T)


