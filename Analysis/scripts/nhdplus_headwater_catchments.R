#################################################################################
# Identify all headwater catchments that intersect the SBR 
library(rgdal)
library(sf)
library(raster)
library(stars)
home <- "/Volumes/GoogleDrive/My Drive/Chapter1_ET_Project/"

# See all of the layers 
st_layers(paste0(home, "Data/nhdPlus/03/NHDPLUS_H_0301_HU4_GDB.gdb"))

# use the SBR polygon to check if the nhd huc4s intersect our area of itnerest 
sbr <- st_read(paste0(home, "Data/NA_CEC_Eco_Level3/blue_ridge.shp"))

# Find which HUC 4 watersheds from each basin intersect with the SBR
# If they do intersect, join the VAA table, and save to a new folder of SBR catchments 
hucs3 <- list.files(paste0(home, "Data/nhdPlus/03"), recursive = F, full.names=T, pattern=".gdb$")
hucs5 <- list.files(paste0(home, "Data/nhdPlus/05"), recursive = F, full.names=T, pattern=".gdb$")
hucs6 <- list.files(paste0(home, "Data/nhdPlus/06"), recursive = F, full.names=T, pattern=".gdb$")
hucs <- c(hucs3, hucs5, hucs6)
isnull_list <- c()
completed <- 0
for(i in 1:length(hucs)){
  wbdhuc4 <- st_read(hucs[i], layer = "WBDHU4")
  wbdhuc4 <- st_transform(wbdhuc4, crs(sbr))
  intt = st_intersection(st_geometry(wbdhuc4), sbr)
  isnull <- is.na(st_bbox(intt)[[1]])
  isnull_list <- c(isnull_list, isnull)
  if(isnull == TRUE){
    next
  }else{
    catchment <- st_read(hucs[i], layer = "NHDPlusCatchment")
    flowlineVAA <- st_read(hucs[i], layer = "NHDPlusFlowlineVAA")
    catchment <- merge(catchment, flowlineVAA, by = "NHDPlusID")
    
    catchment <- catchment[(catchment$HWType == 0) & (!is.na(catchment$HWType == T)), ]
    catchment <- st_transform(catchment, crs(sbr))
    final_catchment <- catchment[unlist(st_intersects(sbr, catchment)),]

    outname <- substr(hucs[i], 67, nchar(hucs[i])-4)
    st_write(final_catchment, paste0(home, "Data/nhdPlus/sbr_headwater_catchments/", outname, ".shp"), append=FALSE)
    completed <- completed + 1
  }
  print(i)
}


headwaters_files <- list.files(paste0(home, "Data/nhdPlus/sbr_headwater_catchments"), full.names=T, pattern = ".shp$")
headwaters <- st_read(headwaters_files[1])
for(i in 2:length(headwaters_files)){
  headd <- st_read(headwaters_files[i])
  headwaters <- rbind(headwaters, headd)
}

# so there are 114,448 headwater catchments identified by NHDplus (i don't think this is the high resolution version)

# The idea would be to calculate TWI for each catchment 
# Get an annual EVI composite for each catchment from JJA from 1984 - 2020 (maybe even 1980?)
# Calculate the slope along the TWI gradient in each year
# And then see if that slope has changed over time 

# Find all gages in the region that include the headwater catchments we're looking at 
# Calculate trends in RR and lowflow baseflow recession 

# This would basically look at where lateral flow is changing the most and connect it directly to changing hydrology 

# Or, another idea, take just the 2007 - 2009 drought 
# Use the ARS ET and calculate the gradient of ET across TWI 
# Calculate on a monthly basis 

