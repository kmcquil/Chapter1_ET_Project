##################################################################################################
##################################################################################################
## Use Riley product to calculate the percent diffuse porous BA in each pixel on the landscape 
## and make maps of the top 50 species BA across the landscape 
##################################################################################################
##################################################################################################
library(raster)
library(gdalUtils)
library(data.table)
library(rgeos)
library(rgdal)
library(sp)
library(tidyr)
library(doParallel)

home <- "G:/My Drive/Chapter1_ET_Project/Data/forest_composition/Riley_FC_SBR"   # home directory 

# bring in teh raster with the FIA plot id in each cell 
tl_id_rast <- raster(paste0(home, "/national_c2014_tree_list.tif"))

# in and out file names and the shapefile for cropping 
srfile <- paste0(home, "/national_c2014_tree_list.tif")
dsfile <- paste0(home, "national_c2014_tree_list_cropped.tif")
spfile <- "G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp"
sp <- readOGR(spfile)
sp1 <- spTransform(sp, crs(tl_id_rast))

# crop and mask the FIA raster to the ROI
new_tl_id_rast <- mask(crop(tl_id_rast, sp1), sp1)
writeRaster(new_tl_id_rast, paste0(home, "/national_c2014_tree_list_cropped.tif"), overwrite = TRUE, format = "GTiff")

# Now get all unique tl_id values from the cropped raster 
new_tl_id_rast <- raster(paste0(home, "/national_c2014_tree_list_cropped.tif"))
un_tl_id <- unique(new_tl_id_rast)

# Bring in the national_c2014_tree_list.txt as a datatable and subset to just the records which correspond to plots 
# located in our study region. Keep only the most recent years records for each tl_id of interest. 
tree_table <- fread(paste0(home, "/Tree_table_CONUS.txt"))
app_tree_table <- tree_table[tree_table$tl_id %in% un_tl_id,]

#make sure that it got all of the plots 
check <- sort(unique(app_tree_table$tl_id))
check1 <- sort(un_tl_id)
check2 <- sum(check == check1)

# go through each tl_id, and find the number of years of data it has represented in this data table 
years <- app_tree_table[,.(number_of_distinct_years = uniqueN(INVYR)), by = tl_id]
un_years <- unique(years$number_of_distinct_years)
# we find that each tl_id includes only one (im assuming the most recent when they were creating this) year of data. 
# This means we don't need to do any further subsetting

# Now on the subset table, use DIA and TP_UNADJ to calculate basal area of each species for each plot. 
# im creating a basal area column at the end of the table. Use the DIA to BA coefficient 0.005454 andn teh TPA_UNADJ expansion factor to findd basal area from DIA
app_tree_table$basal_area <- ((app_tree_table$DIA^2)*0.005454)*app_tree_table$TPA_UNADJ

#now create a new table summing species basal area by plot id with tl_id, species, total basal area 
tl_id_species_sums <- app_tree_table[,.(basal_sum = sum(basal_area, na.rm = TRUE)), by = list(tl_id, SPCD)]
fwrite(tl_id_species_sums, paste0(home, "/tl_id_species_sums.csv"))

# Find which species are most abundant on the landscape
species_sums <- fread(paste0(home, "/tl_id_species_sums.csv"))
species_rast <- raster(paste0(home, "/national_c2014_tree_list_cropped.tif"))
freqs <- freq(species_rast, progress = 'text') 
freqs <- as.data.table(freqs)

# merge the freqs table with the species sums to find how many times each plot id (and associate species basal area) 
# shows up on the landscape. Multiply the basal area by frequency and then sum total basal area per species across the landscape. 
# Find the relative abundance of each species. 
setkey(species_sums, tl_id)
setkey(freqs, value)

species_sums_freq <- merge(species_sums, freqs,by.x = "tl_id", by.y = "value",  all.x = TRUE)
species_sums_freq$total <- species_sums_freq$basal_sum * species_sums_freq$count

landscape_species_basal_area <- species_sums_freq[,.(total_basal = sum(total)), by = SPCD]
total_ba <- sum(landscape_species_basal_area$total_basal)
landscape_species_basal_area$perc_total <- landscape_species_basal_area$total_basal/total_ba
landscape_species_basal_area <- landscape_species_basal_area[order(-landscape_species_basal_area$perc_total),]
landscape_species_basal_area$cummulative <- cumsum(landscape_species_basal_area$perc_total)
fwrite(landscape_species_basal_area, paste0(home, "/landscape_species_basal_area_riley.csv"))

# read in the top 50 species and subset the speices sums to just the top 50 species and then merge them together 
sub_unique_species <- fread(paste0(home, "/top50_species_ba_riley.csv"))

tl_id_species_sums <- fread(paste0(home, "/tl_id_species_sums.csv"))
tl_id_species_sums <- tl_id_species_sums[tl_id_species_sums$SPCD %in% sub_unique_species$SPCD,]

setkey(tl_id_species_sums, SPCD)
setkey(sub_unique_species, SPCD)
tl_id_merged <- merge(tl_id_species_sums, sub_unique_species, all.x = TRUE)

# sum basal area by tl_id
tl_sum <- tl_id_merged[, .(basal_area_sum = sum(basal_sum, na.rm = TRUE)), by = tl_id]

# sum basal area by diffuse, ring, or tracheid 
tl_diff <- tl_id_merged[Type == "diffuse", .(diff_sum = sum(basal_sum, na.rm=T)), by = list(tl_id)]
tl_ring <- tl_id_merged[Type == "ring", .(ring_sum = sum(basal_sum, na.rm=T)), by = list(tl_id)]
tl_trach <- tl_id_merged[Type == "tracheid", .(trach_sum = sum(basal_sum, na.rm=T)), by = list(tl_id)]


# merge with the full tl_id 
setkey(tl_sum, tl_id)
setkey(tl_diff, tl_id)
setkey(tl_ring, tl_id)
setkey(tl_trach, tl_id)

tl_sum_merged <- merge(tl_sum, tl_diff, all.x=T)
tl_sum_merged <- merge(tl_sum_merged, tl_ring, all.x=T)
tl_sum_merged <- merge(tl_sum_merged, tl_trach, all.x=T)

tl_sum_merged$perc_diff <- tl_sum_merged$diff_sum/tl_sum_merged$basal_area_sum
tl_sum_merged$perc_ring <- tl_sum_merged$ring_sum/tl_sum_merged$basal_area_sum
tl_sum_merged$perc_trach <- tl_sum_merged$trach_sum/tl_sum_merged$basal_area_sum

fwrite(tl_sum_merged, paste0(home, "/tl_sum_porosity.csv"))

### fix so i have the full list ids
ugh <- fread(paste0(home, "/tl_sum_porosity.csv"))
full_ids <- fread(paste0(home, "/tl_id_species_sums.csv"))
full_tlid <- as.data.table(unique(full_ids$tl_id))
colnames(full_tlid) <- "tl_id"

full_subset <- merge(full_tlid, ugh, by.x = "tl_id", by.y = "tl_id", all.x = TRUE)
fwrite(full_subset, paste0(home, "/tl_sum_porosity.csv"))

# Associate these values back to the tl_id in each raster cell and create a % diffuse, ring, and tracheid raster 
tl_sum_merged <- fread(paste0(home, "/tl_sum_porosity.csv"))
tl_sum_merged <- as.data.frame(tl_sum_merged)

new_tl_id_rast <- raster(paste0(home, "/national_c2014_tree_list_cropped.tif"))

diff_rast <- reclassify(new_tl_id_rast, as.matrix(tl_sum_merged)[,c(1,6)])
writeRaster(diff_rast, paste0(home,"/diffuse_percent.tif"), overwrite = T, format = "GTiff")

ring_rast <- reclassify(new_tl_id_rast, as.matrix(tl_sum_merged)[,c(1,7)])
writeRaster(ring_rast, paste0(home,"/ring_percent.tif"), overwrite = T, format = "GTiff")

tracheid_rast <- reclassify(new_tl_id_rast, as.matrix(tl_sum_merged)[,c(1,8)])
writeRaster(tracheid_rast, paste0(home,"/tracheid_percent.tif"), overwrite = T, format = "GTiff")



#########################################################################################
###########################################################################################
# Make a raster for each of the top 50 species with the basal area in each cell 
#########################################################################################
sub_unique_species <- fread(paste0(home, "/top50_species_ba_riley.csv"))
full_ids <- fread(paste0(home, "/tl_id_species_sums.csv"))
full_tlid <- as.data.table(unique(full_ids$tl_id))
colnames(full_tlid) <- "tl_id"

tl_id_species_sums <- fread(paste0(home, "/tl_id_species_sums.csv"))
tl_id_species_sums <- tl_id_species_sums[tl_id_species_sums$SPCD %in% sub_unique_species$SPCD,]

# make it the wide way 
tl_id_species_sums_wide <- tl_id_species_sums %>% pivot_wider(names_from = SPCD, values_from = basal_sum)
tl_id_species_sums_wide <- as.data.table(tl_id_species_sums_wide)

# merge back with full_ids so that all of the tl_id are accounted for but we still have only the species of interest 
# this is important when reclassifying the tl_id as basal sum based on spcd 
tl_id_species_sums_wide <- merge(full_tlid, tl_id_species_sums_wide, by = "tl_id", all.x=T, all.y=T)

beginCluster(6)
for(i in 2:ncol(tl_id_species_sums_wide)){
  spec_rast <- reclassify(new_tl_id_rast, as.matrix(tl_id_species_sums_wide)[,c(1,i)])
  writeRaster(spec_rast, paste0(home, "/species/s", colnames(tl_id_species_sums_wide[,..i]), ".tif"), overwrite = T, format = "GTiff")
  print(i)
}
endCluster()


#########################################################################################
###########################################################################################
# Check how many trees are in a plot (and pixel) 
#########################################################################################

# how many records per fia plot 
records_per_plot <- app_tree_table[,.(count = .N), .(tl_id)]
records_per_plot <- records_per_plot[order(records_per_plot$count),]

# match this back to the raster 
records_rast <- reclassify(new_tl_id_rast, as.matrix(records_per_plot))

# resample to match the forest mask layer and mask out pixels that are not continuously forested 
forest_mask_landsat <- raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
forest_mask_landsat[forest_mask_landsat == 0] <- NA
records_rast <- projectRaster(records_rast, forest_mask_landsat, method = "ngb")
records_rast_fm <- mask(records_rast, forest_mask_landsat)

# look at the number of trees in each pixel 
num_records <- values(records_rast_fm)
num_records <- num_records[!is.na(num_records)]

hist(num_records)
min(num_records)
max(num_records)

num_records_df <- as.data.frame(table(num_records))

# how many forested pixels have less than 10 trees in the pixel 
percent_pixels_less10trees <- (sum(num_records_df[1:9, 2])/sum(num_records_df[, 2]))*100
percent_pixels_less20trees <- (sum(num_records_df[1:20, 2])/sum(num_records_df[, 2]))*100
percent_pixels_less50trees <- (sum(num_records_df[1:50, 2])/sum(num_records_df[, 2]))*100


