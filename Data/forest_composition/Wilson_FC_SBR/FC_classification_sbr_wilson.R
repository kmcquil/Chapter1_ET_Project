###########################################################################################################
###########################################################################################################
## Identify the 50 most abundant species on the landscape using the Ty Wilson 250m species BA maps 
## Then calculate the % diffuse porous BA in each pixel based on the top 50 species 
## These maps will be used as a reference to compare with Riley maps 
###########################################################################################################
###########################################################################################################
library(raster)
library(rgdal)
library(readr)
library(foreach)
library(doParallel)
home <- "G:/My Drive/Chapter1_ET_Project/Data/forest_composition/Wilson_FC_SBR"

UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

# bring in the original rasters 
sp_files <- list.files(path = paste0(home, "/Input/RasterMaps"),pattern = "*.img", all.files = T, full.names = T, recursive = T)
apps <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")
apps <- spTransform(apps, crs(raster(sp_files[1])))

# crop them to the SBR study region 
foreach(i = 1:length(sp_files)) %dopar%{
  library(raster)
  library(tools)
  rast <- raster(sp_files[i])
  cut <- mask(crop(rast, apps), apps)
  writeRaster(cut, filename = paste0(home, "/RasterMaps_Cut/",
                                     file_path_sans_ext(basename(sp_files[i])), ".tif"), format="GTiff")
}

# Calculate total BA of each species across the landscape  
sp_files <- list.files(path =  paste0(home,"/RasterMaps_Cut"), pattern = "*.tif", all.files = T, full.names = T, recursive = T)
sp_abundance <- foreach(i = 1:length(sp_files), .combine = rbind)%dopar%{
  library(raster)
  r <- raster(sp_files[i])
  x <- sum(getValues(r), na.rm = TRUE)
  y <- c(names(r), x)
}

# create a DF of the species >0 BA 
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
sp_ab_df <- as.data.frame(sp_abundance)
sp_ab_df$V2 <- as.numeric.factor(sp_ab_df$V2)
sp_ab_df_nonZero <- sp_ab_df[!(sp_ab_df$V2 == 0), ]
ordered <- sp_ab_df_nonZero[order(-sp_ab_df_nonZero$V2),]
ordered$SPCD <- parse_number(as.character(ordered$V1))
 
total_ba <- sum(ordered$V2) # calculate the total BA of all species on landsacape 
landscape_species_basal_area <- data.table(SPCD = ordered$SPCD, total_basal = ordered$V2)
landscape_species_basal_area$perc_total <- (landscape_species_basal_area$total_basal/total_ba)*100 # find the percent of each species of the total 
landscape_species_basal_area$cummulative <- cumsum(landscape_species_basal_area$perc_total) # calculate the cumulative percent explained 
fwrite(landscape_species_basal_area,  paste0(home,"/landscape_species_basal_area_wilson.csv"))

# now create a new .csv based on the landscape_species_basal_area that classifies the top 50 species (accounting for ~98% of basal area) as diffuse, ring, tracheid

# For each category (diffuse, ring, tracheid), sum the basal area on a pixel basis 
# start by grabbing the species classifications 
sp_ba <- fread( paste0(home,"/top50_species_ba_wilson.csv"))
sp_files_short <- list.files(path =  paste0(home,"/RasterMaps_Cut"), pattern = "*.tif", all.files = T, full.names = F, recursive = T)
sp_files <- data.table(file = sp_files, SPCD = parse_number(sp_files_short))

# get the spcd that correspond to diffuse porous species and use to subset the sp_files and create a total diffuse basal area raster 
diff_spcd <- sp_ba[Type == "diffuse"]$SPCD
diff_files <- sp_files[SPCD %in% diff_spcd]$file
diff_stack <- do.call("stack", lapply(diff_files, raster))
diff_sum <- sum(diff_stack, na.rm= T)

# do the same for ring 
ring_spcd <- sp_ba[Type == "ring"]$SPCD
ring_files <- sp_files[SPCD %in% ring_spcd]$file
ring_stack <- do.call("stack", lapply(ring_files, raster))
ring_sum <- sum(ring_stack, na.rm= T)

# do the same for tracheid 
trach_spcd <- sp_ba[Type == "tracheid"]$SPCD
trach_files <- sp_files[SPCD %in% trach_spcd]$file
trach_stack <- do.call("stack", lapply(trach_files, raster))
trach_sum <- sum(trach_stack, na.rm= T)

# sum all of them to get total basal area
all_stack <- stack(diff_sum, ring_sum, trach_sum)
all_sum <- sum(all_stack, na.rm=T)

# now get the percent BA of diffuse, ring, and tracheid in each cell 
full_stack <- stack(all_sum, diff_sum, ring_sum, trach_sum)

diff_prc <- calc(full_stack, function(x) { x[2]/x[1] })
writeRaster(diff_prc,  paste0(home,"/diffuse_percent.tif"), format = "GTiff", overwrite = T)

ring_prc <- calc(full_stack, function(x) { x[3]/x[1] })
writeRaster(ring_prc,  paste0(home,"/ring_percent.tif"), format = "GTiff", overwrite = T)

trach_prc <- calc(full_stack, function(x) { x[4]/x[1] })
writeRaster(trach_prc,  paste0(home,"/tracheid_percent.tif"), format = "GTiff", overwrite = T)


