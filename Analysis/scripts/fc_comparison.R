#########################################################################################################
## Validate the Riley forest composition maps (% diffuse porous BA) using Ty Wilson maps as a reference  
## 1. Compare the predicted percent diffuse porous BA on a pixel basis  
## 2. Compare the top 50 and top 10 species on the landscape 
## 3. Compare the total BA to estimates from the literature 
#########################################################################################################
library(raster)
library(gdalUtils)
library(data.table)
home <- "G:/My Drive/Chapter1_ET_Project/Data/forest_composition"
##########################################################################################################
## Compare % diffuse basal area between Riley and Wilson 
###########################################################################################################

# Bring in the Wilson % diffuse porous BA raster
diff_rast_w1 <- raster(paste0(home, "/Wilson_FC_SBR/diffuse_percent.tif"))

# resample the Riley 30m % diffuse porous BA raster to match the Wilson raster to facilitate comparison
# this means resampling to 250 m resolution using bilinear method 
res <- res(diff_rast_w1)
res_out <- crs(diff_rast_w1)
t1 <- c(xmin(diff_rast_w1), ymin(diff_rast_w1), 
        xmax(diff_rast_w1), ymax(diff_rast_w1))

gdalwarp(paste0(home, "/Riley_FC_SBR/diffuse_percent.tif"), 
         paste0(home, "/Riley_FC_SBR/diffuse_percent_resampled_tw.tif"),
         t_srs=res_out, tr=res,te = t1, r='bilinear', overwrite = T)

# bring in the resampled riley % diffuse porous raster 
diff_rast <- raster(paste0(home, "/Riley_FC_SBR/diffuse_percent_resampled_tw.tif"))

# Find the difference, convert to epsg 32617 to make the map consistent with other figures, and multiply by 100 to get percent 
d_diff <- diff_rast_w1 - diff_rast
d_diff <- projectRaster(d_diff, crs = "EPSG:32617")
d_diff <- d_diff*100

sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp") # bring in the blue ridge ecoregion outline to use as an outline in the figure 
quants <- quantile(c(values(d_diff)), c(0.02, 0.98), na.rm=T) # to saturate the map when plotting use a 2% stretch 

# create a color palette that diverges at 0 
col5 <- colorRampPalette(c("#d73027", 'gray96', "#313695")) #create color ramp starting from blue to red
color_levels=50  #the number of colors to use
max_absolute_value=max(abs(c(quants[[1]], quants[[2]]))) #what is the maximum absolute value of raster?
color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)

# function to add a label to each sub plot 
add_label_legend <- function(x, y, label, ...) {
  legend(x,y, label, bty = "n", ...)
}

# plot the difference in percent diffuse porous BA. Add the legend last
tiff("G:/My Drive/Chapter1_ET_Project/Figures/FigS1.tiff", units="in", width=5.5, height=9, res=800)
par(mfrow = c(2,1), mai = c(0, 0, 0, 0))
plot(d_diff, col=col5(n=color_levels), breaks=color_sequence, 
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), labels=pretty(quants[[1]]:quants[[2]])), 
     zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, 
     cex.lab = 1.65, 
     cex.axis=1.5, 
     font = 2,
     font.lab = 2)
plot(sbr, add = T)
add_label_legend(180000, 4120000, "a", cex = 1.25, text.font = 2)

par(mai=c(0,2, 0.3, 0),new = F)# add the legend 
plot(d_diff, legend.only=TRUE, legend.shrink=0.9, legend.width=1, 
     zlim=c(quants[[1]], quants[[2]]),
     col=col5(n=color_levels), breaks=color_sequence,
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), labels=pretty(quants[[1]]:quants[[2]])),
     legend.args=list(text='', side=4, font=2,line=2.3))


# make a histogram of the difference in % diff porous BA 
par(mai=c(0.85,0.75,0.25,1.35))
hist(values(d_diff), xlab = "Wilson % BA - Riley % BA", ylab = "",main = "", col= "#313695", xlim = c(-50, 50), 
     cex.lab = 1, 
     cex.axis=1, 
     font = 1,
     font.lab = 1)
add_label_legend(-60, 95000, "b", cex = 1.25, text.font = 2)

dev.off()

##########################################################################################################
## Find how many of Riley top 10 and top 50 species are also identified by Wilson 
###########################################################################################################
riley_top50 <- fread(paste0(home, "/Riley_FC_SBR/top50_species_ba_riley.csv"))
wilson_top50 <- fread(paste0(home, "/Wilson_FC_SBR/top50_species_ba_wilson.csv"))

# how many of the riley are in the wilson top 50 
w_in <- wilson_top50$SPCD[wilson_top50$SPCD %in% riley_top50$SPCD]
length(w_in) # 45 

# how many are in the top 10 
w10 <- wilson_top50[1:10,]
r10 <- riley_top50[1:10,]
w_in10 <- w10$SPCD[w10$SPCD %in% r10$SPCD]
length(w_in10) # 8 



##########################################################################################################
## Calculate the total BA (of top 50 species) for Wilson at 250m resolution and Riley at 30m 
## and compare to each other and visually to values from the literature 
###########################################################################################################

# get all of the Wilson rasters that correspond to species in the top 50 
wilson_sp <- list.files(paste0(home, "/Wilson_FC_SBR/species"), full.names=T, pattern = ".tif$")
top50_wilson <- fread(paste0(home, "/Wilson_FC_SBR/top50_species_ba_wilson.csv"))
top50_wilson$name <- paste0("s", top50_wilson$SPCD, ".tif")
wilson_sp <- wilson_sp[grep(paste(top50_wilson$name, collapse = "|"), wilson_sp)]
# sum the BA of the top 50 wilson species 
wilson_stk <- do.call("stack", lapply(wilson_sp, raster))
wilson_sum <- sum(wilson_stk, na.rm=T)
writeRaster(wilson_sum, paste0(home, "/Wilson_FC_SBR/wilson_sum_BA.tif"), overwrite = T, format ="GTiff")


# get the riley rasters that correspond to species in the top 50 
riley_sp30 <- list.files(paste0(home, "/Riley_FC_SBR/species"), full.names=T, pattern = ".tif$")
# sum them using a cluster because they're a lot bigger and loop because i can't bring all of that into ram at once 
r30_sum <- raster(riley_sp30[1])
beginCluster(5)
for(i in 2:length(riley_sp30)){
  r <- raster(riley_sp30[i])
  r30_sum <- sum(r30_sum, r, na.rm=T)
  print(i)
}
writeRaster(r30_sum, paste0(home, "/Riley_FC_SBR/riley_30_sum_BA.tif"), overwrite = T, format ="GTiff")

