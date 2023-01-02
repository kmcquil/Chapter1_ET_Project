#################################################################################
# Calculate the length of headwater streams in the high elevation areas 
# Calculate the percent of total streams in the sbr in high elevation areas 
# This will help us say more about the impact on streamflow 
library(rgdal)
library(sf)
library(raster)
library(stars)
home <- "/Volumes/GoogleDrive/My Drive/Chapter1_ET_Project/"

# See all of the layers 
st_layers(paste0(home, "Data/nhdPlus/03/NHDPLUS_H_0301_HU4_GDB.gdb"))

# use the SBR polygon to check if the nhd huc4s intersect our area of itnerest 
sbr <- st_read(paste0(home, "Data/NA_CEC_Eco_Level3/blue_ridge.shp"))

# Using the elevation raster, create a shapefile of all pixels > 1000m 
elevation <- raster(paste0(home, "Data/Topography/usgsNED_elevation/elevation30m.tif"))
high_elevation <- elevation
high_elevation[high_elevation < 1000] <- NA
high_elevation[!is.na(high_elevation)] <- 1000
high_elevation_stars <- st_as_stars(high_elevation)
high_elevation_sf <- st_as_sf(high_elevation_stars, as_points=FALSE, merge=TRUE)
high_elevation_sf <- st_make_valid(high_elevation_sf) 

# Find which HUC 4 watersheds from each basin intersect with the SBR
# If they do intersect, join the VAA table, and save to a new folder of SBR burnlines 
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
    burnline <- st_read(hucs[i], layer = "NHDPlusBurnLineEvent")
    flowlineVAA <- st_read(hucs[i], layer = "NHDPlusFlowlineVAA")
    burnline <- merge(burnline, flowlineVAA, by = "NHDPlusID")
    burnline <- st_transform(burnline, crs(sbr))
    final_lines <- st_intersection(burnline, sbr)
    outname <- substr(hucs[i], 67, nchar(hucs[i])-4)
    st_write(final_lines, paste0(home, "Data/nhdPlus/sbr_burnlines/", outname, ".shp"), append=FALSE)
    completed <- completed + 1
  }
}


# Calculate total km of all streams
# Calculate total km of all headwater streams 
# Calculate the total km of all streams in high elevations 
# Calculate the total km of all headwater streams at high elevation 
burnline_files <- list.files(paste0(home, "Data/nhdPlus/sbr_burnlines"), full.names=T, pattern = ".shp$")
burnlines <- st_read(burnline_files[1])
for(i in 2:length(burnline_files)){
  burn <- st_read(burnline_files[i])
  burnlines <- rbind(burnlines, burn)
}

burnlines_all <- burnlines
burnlines_all$length <- st_length(burnlines)

burnlines_high <- st_intersection(burnlines, high_elevation_sf)
burnlines_high$length <- st_length(burnlines_high)

headwaters_all <- burnlines[(burnlines$HWType == 0) & (!is.na(burnlines$HWType == T)), ]
headwaters_all$length <- st_length(headwaters_all)

headwaters_high <- st_intersection(burnlines[(burnlines$HWType == 0) & (!is.na(burnlines$HWType == T)), ], high_elevation_sf)
headwaters_high$length <- st_length(headwaters_high) 

# How many km of headwater streams are located at high elevation? = 5941.727 km 
length_high_elevation_headwaters <- sum(as.numeric(headwaters_high$length))/1000

# What percentage of high elevations streams are headwater streams? = 66.29%
pct_high_headwaters_all <- (sum(as.numeric(headwaters_high$length))/sum(as.numeric(burnlines_high$length)))*100

# What percentage of all streams in the sbr are headwater streams? = 49.93%
pct_headwaters_all <- (sum(as.numeric(headwaters_all$length))/sum(as.numeric(burnlines_all$length)))*100

# What percentage of all streams in the sbr are high elevation headwater streams = 7.33%
pct_highheadwaters_all <- (sum(as.numeric(headwaters_high$length))/sum(as.numeric(burnlines_all$length)))*100








