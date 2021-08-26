###########################################################################################
## Validate the forest composition maps 
## compare classification and species 
###########################################################################################
library(raster)
library(gdalUtils)
library(data.table)


# Resample the Riley 30m FC classification to match the ty wilson FC classification
# here are the ty wilson products 
diff_rast_w1 <- raster("G:/My Drive/Dissertation/Forest Composition/FC_SBR/diffuse_percent.tif")
ring_rast_w1 <- raster("G:/My Drive/Dissertation/Forest Composition/FC_SBR/ring_percent.tif")
trach_rast_w1 <- raster("G:/My Drive/Dissertation/Forest Composition/FC_SBR/tracheid_percent.tif")
classified_w1 <- raster("G:/My Drive/Dissertation/Forest Composition/FC_SBR/FC_classified.tif")

res <- res(diff_rast_w1)
res_out <- crs(diff_rast_w1)
t1 <- c(xmin(diff_rast_w1), ymin(diff_rast_w1), 
        xmax(diff_rast_w1), ymax(diff_rast_w1))

gdalwarp("G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent.tif", 
         "G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent_resampled_tw.tif",
         t_srs=res_out, tr=res,te = t1, r='bilinear', overwrite = T)

gdalwarp("G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent.tif", 
         "G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent_resampled_tw.tif",
         t_srs=res_out, tr=res,te = t1, r='bilinear', overwrite = T)

gdalwarp("G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent.tif", 
         "G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent_resampled_tw.tif",
         t_srs=res_out, tr=res,te=t1, r='bilinear', overwrite = T)


# bring in the resampled riley rasters and classify as diffuse, ring, tracheid 
diff_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent_resampled_tw.tif")
ring_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent_resampled_tw.tif")
trach_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent_resampled_tw.tif")
class_stack <- stack(diff_rast, ring_rast, trach_rast)

por_fun <- function(x){
  ifelse(x[1] >= 0.5, 1, 
         ifelse(x[2] >= 0.5, 2, 
                ifelse(x[3] >= 0.5, 3, 
                       ifelse(is.na(x[1:3]) == T, NA, 
                              ifelse(x[1] < 0.5 & x[2] < 0.5 & x[3] < 0.5, 4, NA)))))
  
}

# this is the riley product 
classified_r2 <- calc(class_stack, por_fun)

# this is the wilson product 
classified_stack <- stack(classified_w1, classified_r2)
cl_st_df <- as.data.frame(classified_stack)
colnames(cl_st_df) <- c("Wilson", "Riley")

cl_st_df$Riley <- as.factor(cl_st_df$Riley)
cl_st_df$Wilson <- as.factor(cl_st_df$Wilson)
library(caret)
cl_cm <- confusionMatrix(cl_st_df$Riley, cl_st_df$Wilson)



# plot the original riley classification, resampled riley (250), and the wilson 
classified_r1 <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/FC_classified.tif")
par(mfrow = c(1,3))
plot(classified_w1, main = "Wilson")
plot(classified_r1, main = "Riley (30m)")
plot(classified_r2, main = "Riley (250m)")


# how many of the top 50 species do they have in common 
riley_top50 <- fread("G:/My Drive/Dissertation/species_data/FC_SBR/top50_species_ba_riley.csv")
wilson_top50 <- fread("G:/My Drive/Dissertation/Forest Composition/FC_SBR/top50_species_ba_wilson.csv")

### how many of the riley are in the wilson top 50 
w_in <- wilson_top50$SPCD[wilson_top50$SPCD %in% riley_top50$SPCD]
length(w_in)

w10 <- wilson_top50[1:10,]
r10 <- riley_top50[1:10,]

w_in10 <- w10$SPCD[w10$SPCD %in% r10$SPCD]
length(w_in10)



# reclassify riley resampled (250m) product using a simple majority 
maj_fun <- function(x){
  if(sum(is.na(x[1:3])) == 3){
    y <- NA
  }else{
    y <- which.max(x[1:3])
  }
}


# bring in the resampled riley rasters and classify as diffuse, ring, tracheid 
diff_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/diffuse_percent_resampled_tw.tif")
ring_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/ring_percent_resampled_tw.tif")
trach_rast <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/tracheid_percent_resampled_tw.tif")
class_stack_r <- stack(diff_rast, ring_rast, trach_rast)

classified_r2_simple_majority <- calc(class_stack_r, maj_fun)
classified_w1_simple_majority <- raster("G:/My Drive/Dissertation/Forest Composition/FC_SBR/FC_classified_simple_majority.tif")

classified_stack_sm <- stack(classified_w1_simple_majority, classified_r2_simple_majority)
cl_st_df_sm <- as.data.frame(classified_stack_sm)
colnames(cl_st_df_sm) <- c("Wilson", "Riley")

cl_st_df_sm$Riley <- as.factor(cl_st_df_sm$Riley)
cl_st_df_sm$Wilson <- as.factor(cl_st_df_sm$Wilson)
library(caret)
cl_cm_sm <- confusionMatrix(cl_st_df_sm$Riley, cl_st_df_sm$Wilson)

r1_sm <- raster("G:/My Drive/Dissertation/species_data/FC_SBR/FC_classified_simple_majority.tif")
par(mfrow = c(1,3))
plot(classified_w1_simple_majority, main = "Wilson") 
plot(r1_sm, main = "Riley (30m)")
plot(classified_r2_simple_majority, main = "Riley (250m)")



# plot the difference in percent ba for each group 

d_diff <- diff_rast_w1 - diff_rast
d_ring <- ring_rast_w1 - ring_rast
d_trach <- trach_rast_w1 - trach_rast

library(RColorBrewer)
colss <- brewer.pal(100, "RdYlBu")
par(mfrow = c(1,3))
plot(d_diff, main = "Diffuse", col = colss) 
plot(d_ring, main = "Ring", col = colss)
plot(d_trach, main = "Tracheid", col = colss)

