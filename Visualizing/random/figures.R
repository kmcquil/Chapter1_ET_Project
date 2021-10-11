#############################################################################################
## Make figures that combine MODIS and Landsat results in the same plot 
##############################################################################################
library(rasterVis)
library(raster)
library(pals)
library(RColorBrewer)
library(trend)
library(data.table)
library(viridis)
library(ggplot2)
library(rgdal)
library(cowplot)
home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))


####################################################################################################################
#### Plot anomalies at drought peak for all, moderate, severe, and extreme droughts for MODIS and Landsat 
#####################################################################################################################
# Bring in the SBR shapefile to plot as outline of region on map 
sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")

# Stack MODIS anomalies at drought peak (all, moderate, severe, extreme)
et_tifs_m <- list.files(paste0(home, "/Analysis/outputs/MODIS/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs_m <- et_tifs_m[grep("ETanom", et_tifs_m)]
anom_stk_m <- do.call("stack", lapply(anom_tifs_m, raster))

# Do the same for Landsat 
et_tifs_l <- list.files(paste0(home, "/Analysis/outputs/Landsat/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs_l <- et_tifs_l[grep("ETanom", et_tifs_l)]
anom_stk_l <- do.call("stack", lapply(anom_tifs_l, raster))
# Non-forested pixels haven't beens creened out of the .tifs yet for MODIS so do that here 
forest_mask_landsat <- raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
forest_mask_landsat[forest_mask_landsat == 0] <- NA
anom_stk_l <- mask(anom_stk_l, forest_mask_landsat)

# use a 2% stretch from modis and landsat anomaly values to make maps look nicer/easier to read 
quants <- quantile(c(values(anom_stk_m), values(anom_stk_l)), c(0.02, 0.98), na.rm=T)

# create color palette 
col5 <- colorRampPalette(c("#d73027", 'gray96', "#313695"))  #create color ramp starting from blue to red
color_levels=50 #the number of colors to use
max_absolute_value=max(abs(c(quants[[1]], quants[[2]]))) #what is the maximum absolute value of raster?
color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)

# plot the MODIS and then Landsat maps without a legend 
par(mfrow = c(2,4), mai = c(0, 0, 0.2, 0), omi=c(0,0.1,0.2,1))
plot(anom_stk_m[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "All", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Moderate", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Severe", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Extreme", cex.main = 2)
plot(sbr, add = T)

plot(anom_stk_l[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)

# add a big legend on the side 
par(mfrow=c(1, 1), mai=c(0.5, 6, 0.5, 0),omi=c(0,0.8,0,0), new=FALSE)
plot(anom_stk_m[[1]], legend.only=TRUE, legend.shrink=1, legend.width=1, 
     zlim=c(quants[[1]], quants[[2]]),
     col=col5(n=color_levels), breaks=color_sequence,
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), font = 2,labels=pretty(quants[[1]]:quants[[2]])),
     legend.args=list(text='', side=4, font=2, line=2.3))


####################################################################################################################
#### Plot Elevation, HAND, and % diffuse porous BA at 30m resolution 
#####################################################################################################################
elevation <- raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation30m.tif"))
hand <- raster(paste0(home, "/Data/Topography/HeightAboveNearestDrainage/hand_landsat.tif"))
pctDiffR <- raster(paste0(home, "/Data/forest_composition/riley_landsat_diffuse_percent.tif"))
pctDiffR <- pctDiffR*100 # convert from decimal to percent to match other figures 
par(mfrow = c(1,3), mai = c(0, 0.25, 0, 0), omi=c(0,0.1,0,0.25))
plot(elevation, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))
plot(hand, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))
plot(pctDiffR, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))



####################################################################################################################
#### Plot the percent of anomalies > 0 at all, moderate, severe, and extreme drought levels  
#####################################################################################################################
# get file where the % greater than 0 has already been calcualted for MODIS and Landsat 
modis_files <- list.files(paste0(home, "/Visualizing/MODIS/"), full.names=T, pattern = ".csv$")
landsat_files <- list.files(paste0(home, "/Visualizing/Landsat/"), full.names=T, pattern = ".csv$")
landsat_files <- landsat_files[grep("pct_greater0", landsat_files)]

# bring together the first file for landsat and modis to plot the % of anomalies greater than 0 
pg0_l <- fread(landsat_files[1])
pg0_l$sensor <- rep("Landsat", 5)
pg0_m <- fread(modis_files[1])
pg0_m$sensor <- rep("MODIS", 5)
pg0 <- setDT(rbind(pg0_l, pg0_m))
pg0$DroughtSeverity <- factor(pg0$DroughtSeverity, ordered=TRUE, 
                                 levels = c("Moderate", "Severe", "Extreme", "All", "None"))
pg0_sub <- pg0[DroughtSeverity == "Moderate" |DroughtSeverity == "Severe" | DroughtSeverity == "Extreme", ]
pg0_sub$pg0 <- pg0_sub$pg0*100
# Make a bar plot 
ggplot(pg0_sub, aes(x = DroughtSeverity, y = pg0, fill = sensor)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic()+
  xlab("Drought Severity") + 
  ylab(expression(bold(`%`~ET["dp"]~`>`~0))) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=12))+
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65))
  


#################################################################################################################
## Make plots of the binned overall correlation and the binned trend 
## across elevation, HAND, and % diffuse porous gradients for MODIS and Landsat 
#################################################################################################################
dir_l <- "G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat"
dir_m <- "G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS"

# Start with overall coupling 
### first do elevation
elev_corr_in <- prep(paste0(dir_m, "/elevation_allCoupling.csv"), paste0(dir_l, "/elevation_allCoupling.csv"), "elevation_bin")
elev_correlation <- ggplot(data =elev_corr_in[[1]], aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("Elevation (m)") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 725, y = c(-0.3, -0.45), label = c(elev_corr_in[[2]], elev_corr_in[[3]]))
  

### now do HAND 
hand_corr_in <- prep(paste0(dir_m, "/hand_allCoupling.csv"), paste0(dir_l, "/hand_allCoupling.csv"), "HAND_bin")
hand_correlation <- ggplot(data = hand_corr_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("HAND (m)") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 125, y = c(-0.17, -0.27), label = c(hand_corr_in[[2]], hand_corr_in[[3]]))

### now do percent diffuse porous 
diff_corr_in <- prep(paste0(dir_m, "/pctDiffR_allCoupling.csv"), paste0(dir_l, "/pctDiffR_allCoupling.csv"), "pctDiff_R_bin")
diff_correlation <- ggplot(data = diff_corr_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("% Diffuse Porous") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.05, -0.12), label = c(diff_corr_in[[2]], diff_corr_in[[3]]))


# start the sen's slope plots 
elev_trend_in <- prep(paste0(dir_m, "/elevation_sensCoupling.csv"), paste0(dir_l, "/elevation_sensCoupling.csv"), "elevation_bin")
elev_trend <- ggplot(data = elev_trend_in[[1]], aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("Elevation (m)") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 600, y = c(-0.02, -0.028), label = c(elev_trend_in[[2]], elev_trend_in[[3]]))



### next do HAND 
hand_trend_in <- prep(paste0(dir_m, "/HAND_sensCoupling.csv"), paste0(dir_l, "/HAND_sensCoupling.csv"), "HAND_bin")
hand_trend <- ggplot(data = hand_trend_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("HAND (m)") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 100, y = c(-0.022, -0.028), label = c(hand_trend_in[[2]], hand_trend_in[[3]]))


### next do % diffuse porous 
diff_trend_in <- prep(paste0(dir_m, "/pctDiffR_sensCoupling.csv"), paste0(dir_l, "/pctDiffR_sensCoupling.csv"), "pctDiff_R_bin")
diff_trend <- ggplot(data = diff_trend_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("%Diffuse Porous") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.023, -0.029), label = c(diff_trend_in[[2]], diff_trend_in[[3]]))

plot_grid(elev_correlation, elev_trend, 
          hand_correlation, hand_trend, 
          diff_correlation, diff_trend, 
          labels = "AUTO", nrow = 3)





############################################################################################
## Make plots of the ETdp across gradients broken down by drought severity 
############################################################################################
# start with landsat get the .csv and make the plots individually 
l_e <- contPlots_1(paste0(dir_l, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
l_d <- contPlots_1(paste0(dir_l, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
l_h <- contPlots_1(paste0(dir_l, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "HAND (m)")

# then do modis 
m_e <- contPlots_1(paste0(dir_m, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
m_d <- contPlots_1(paste0(dir_m, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
m_h <- contPlots_1(paste0(dir_m, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "HAND (m)")

# Add the sen's slope of M, S, and E drought with a * for < 0.05 and ** < 0.01, and *** < 0.001
# do elevation first for MODIS and Landsat 
le_labels <- makeLabels(l_e[[2]])
LE <- l_e[[1]] + annotate("text", x = 1500, y = c(-0.61, -0.85, -1.07), label = le_labels)
me_labels <- makeLabels(m_e[[2]])
ME <- m_e[[1]] + annotate("text", x = 1250, y = c(-1, -1.2, -1.4), label = me_labels)

# Next do Hand 
lh_labels <- makeLabels(l_h[[2]])
LH <- l_h[[1]] + annotate("text", x = 350, y = c(-0.33, -0.48, -0.63), label = lh_labels)
mh_labels <- makeLabels(m_h[[2]])
MH <- m_h[[1]] + annotate("text", x = 110, y = c(-0.88, -1.03, -1.2), label = mh_labels)


# Finally do % Diff porous BA
ld_labels <- makeLabels(l_d[[2]])
LD <- l_d[[1]] + annotate("text", x = 70, y = c(-0.3, -0.45, -0.60), label = ld_labels)
md_labels <- makeLabels(m_d[[2]])
MD <- m_d[[1]] + annotate("text", x = 55, y = c(-0.9, -1.1, -1.3), label = md_labels)

plot_grid(ME, LE, 
          MH, LH, 
          MD, LD,
          labels = "AUTO", nrow = 3)


