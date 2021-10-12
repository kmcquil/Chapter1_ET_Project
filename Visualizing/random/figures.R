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
library(ggpubr)
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

# initiate the save to a .tif file for the full plot 
tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_3.tiff", units="in", width=5.5, height=7, res=800)

# plot the MODIS and then Landsat maps without a legend 
par(mfrow = c(4,2), mai = c(0, 0.25, 0.2, 0), omi=c(0.05,0,0,1))
#all
image(anom_stk_m[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, box = F, axes=F)
mtext(side = 3, line = 0, "MODIS", cex=1)
mtext(side = 2, line = 0.2,  "All", cex = 1)plot(sbr, add = T)
image(anom_stk_l[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
mtext(side = 3, line = 0, "Landsat", cex=1)
plot(sbr, add = T)

#moderate
image(anom_stk_m[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
mtext(side = 2, line = 0.2,  "Moderate", cex = 1)
image(anom_stk_l[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)

#severe
image(anom_stk_m[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
mtext(side = 2, line = 0.2,  "Severe", cex = 1)
image(anom_stk_l[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)

#extreme
image(anom_stk_m[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
mtext(side = 2, line = 0.2,  "Extreme", cex = 1)
image(anom_stk_l[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)

# add a big legend on the side 
par(mfrow=c(1, 1), mai=c(0.1, 6, 0.1, 0),omi=c(0,0.8,0,0), new=FALSE)
plot(anom_stk_m[[1]], legend.only=TRUE, legend.shrink=1, legend.width=1, 
     zlim=c(quants[[1]], quants[[2]]),
     col=col5(n=color_levels), breaks=color_sequence,
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), font = 1,cex = 1, labels=pretty(quants[[1]]:quants[[2]])),
     legend.args=list(text='', side=4, font=1, line=2.3))

dev.off()

####################################################################################################################
#### Plot Elevation, TWI, and % diffuse porous BA at 30m resolution 
#####################################################################################################################
elevation <- raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation30m.tif"))
TWI <- raster(paste0(home, "/Data/Topography/TWI/TWI_landsat.tif"))
pctDiffR <- raster(paste0(home, "/Data/forest_composition/riley_landsat_diffuse_percent.tif"))
pctDiffR <- pctDiffR*100 # convert from decimal to percent to match other figures 

# function to add a label to each sub plot 
add_label_legend <- function(x, y, label, ...) {
  legend(x,y, label, bty = "n", ...)
}

tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_S2.tiff", units="in", width=7, height=3, res=800)
par(mfrow = c(1,3), mai = c(0, 0.2, 0, 0.2), omi=c(0,0.1,0,0.1))
plot(elevation, col = viridis(100), axes=F, box = F, legend.width = 1.5)
add_label_legend(90500, 4180565, "A", cex = 1.25, text.font = 2)
plot(TWI, col = viridis(100), axes=F, box = F, legend.width = 1.5)
add_label_legend(90500, 4180565, "B", cex = 1.25, text.font = 2)
plot(pctDiffR, col = viridis(100), axes=F, box = F, legend.width = 1.5)
add_label_legend(90500, 4180565, "C", cex = 1.25, text.font = 2)
dev.off()


####################################################################################################################
#### Plot the percent of anomalies > 0 at all, moderate, severe, and extreme drought levels  
#####################################################################################################################
# get file where the % greater than 0 has already been calcualted for MODIS and Landsat 
modis_files <- list.files(paste0(home, "/Visualizing/MODIS/"), full.names=T, pattern = ".csv$")
modis_files <- modis_files[grep("pct_greater0", modis_files)]
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
tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_4.tiff", units="in", width=6, height=5, res=800)
ggplot(pg0_sub, aes(x = DroughtSeverity, y = pg0, fill = sensor)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic()+
  xlab("Drought Severity") + 
  ylab(expression(bold(`%`~ET["dp"]~`>`~0))) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=12), 
        legend.position = c(0.93, 0.92))+
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12, face ='bold')) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65))
  
dev.off()

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
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
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
  

### now do HAND ---- really TWI 
hand_corr_in <- prep(paste0(dir_m, "/hand_allCoupling.csv"), paste0(dir_l, "/hand_allCoupling.csv"), "HAND_bin")
hand_correlation <- ggplot(data = hand_corr_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("TWI") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 20, y = c(0, 0.05), label = c(hand_corr_in[[2]], hand_corr_in[[3]]))

### now do percent diffuse porous 
diff_corr_in <- prep(paste0(dir_m, "/pctDiffR_allCoupling.csv"), paste0(dir_l, "/pctDiffR_allCoupling.csv"), "pctDiff_R_bin")
diff_correlation <- ggplot(data = diff_corr_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
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
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
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



### next do HAND --- really TWI
hand_trend_in <- prep(paste0(dir_m, "/HAND_sensCoupling.csv"), paste0(dir_l, "/HAND_sensCoupling.csv"), "HAND_bin")
hand_trend <- ggplot(data = hand_trend_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("TWI") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 20, y = c(-0.018, -0.022), label = c(hand_trend_in[[2]], hand_trend_in[[3]]))


### next do % diffuse porous 
diff_trend_in <- prep(paste0(dir_m, "/pctDiffR_sensCoupling.csv"), paste0(dir_l, "/pctDiffR_sensCoupling.csv"), "pctDiff_R_bin")
diff_trend <- ggplot(data = diff_trend_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
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

## make a plot with just the legend -- no plot 
legend <- ggplot(data = diff_trend_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2, color=NA) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("%Diffuse Porous") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.023, -0.029), label = c(diff_trend_in[[2]], diff_trend_in[[3]]))
legend <- ggpubr::get_legend(legend)

tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_7.tiff", units="in", width=8, height=8, res=800)
ggarrange(elev_correlation, elev_trend, 
          hand_correlation, hand_trend, 
          diff_correlation, diff_trend, 
          labels = "AUTO", nrow = 3, ncol = 2, 
          common.legend = TRUE, 
          legend.grob = legend,
          legend = "bottom")
dev.off()


############################################################################################
## Make plots of the ETdp across gradients broken down by drought severity 
############################################################################################
# start with landsat get the .csv and make the plots individually 
l_e <- contPlots_1(paste0(dir_l, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
l_d <- contPlots_1(paste0(dir_l, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
l_h <- contPlots_1(paste0(dir_l, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "TWI")

# then do modis 
m_e <- contPlots_1(paste0(dir_m, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
m_d <- contPlots_1(paste0(dir_m, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
m_h <- contPlots_1(paste0(dir_m, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "TWI")

# Add the sen's slope of M, S, and E drought with a * for < 0.05 and ** < 0.01, and *** < 0.001
# do elevation first for MODIS and Landsat 
le_labels <- makeLabels(l_e[[2]])
LE <- l_e[[1]] + annotate("text", x = 1500, y = c(-0.61, -0.85, -1.07), label = le_labels)
me_labels <- makeLabels(m_e[[2]])
ME <- m_e[[1]] + annotate("text", x = 1250, y = c(-1, -1.2, -1.4), label = me_labels)

# Next do Hand 
lh_labels <- makeLabels(l_h[[2]])
LH <- l_h[[1]] + annotate("text", x = 20, y = c(-0.55, -0.75, -0.95), label = lh_labels)
mh_labels <- makeLabels(m_h[[2]])
MH <- m_h[[1]] + annotate("text", x = 11, y = c(-0.9, -1.06, -1.22), label = mh_labels)


# Finally do % Diff porous BA
ld_labels <- makeLabels(l_d[[2]])
LD <- l_d[[1]] + annotate("text", x = 70, y = c(-0.3, -0.45, -0.60), label = ld_labels)
md_labels <- makeLabels(m_d[[2]])
MD <- m_d[[1]] + annotate("text", x = 55, y = c(-0.9, -1.1, -1.3), label = md_labels)

m_h_legend <- contPlots_legend(paste0(dir_m, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "HAND (m)")
m_h_legend <- ggpubr::get_legend(m_h_legend)

tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_5.tiff", units="in", width=9, height=7, res=800)
ggarrange(ME, LE,  
          MH, LH, 
          MD, LD,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend.grob = m_h_legend,
          legend="bottom", 
          labels = "AUTO")
dev.off()



#########################################################################################################################
# Figure showing the percent of forested area that has a significant overall positive vs negative coupling vs non significant 
# And another figure showing the same thing but for a significant positive vs negative change in coupling vs non significant 
##########################################################################################################################

sensor <- c(rep("MODIS", 6), rep("Landsat", 6))
class <- rep(c("Non-Significant", "Positive", "Negative"), 4)
metric <- c(rep("R(SPI~ET)", 3), rep("Trend in R(SPI~ET)", 3), rep("R(SPI~ET)", 3), rep("Trend in R(SPI~ET)", 3))
value <- c(53.82, 45.70, 0.47, 
           86.78, 4.29, 8.93, 
           83.94, 14.56, 1.50, 
           98.83, 0.32, 0.85)
perc_df <- data.frame(sensor = sensor, class = class, metric = metric, value = value)
perc_df$group2 <- paste0(perc_df$sensor, "+", perc_df$metric)
perc_df$class <- as.factor(perc_df$class)
perc_df$class <- ordered(perc_df$class, levels = c("Positive", "Negative", "Non-Significant"))

# make a bar plot where each bar is broken into a percentage corresponding to non-significant , positive, or negative 
pal <- c("#30678D","#36B677", "#440154")

tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_6.tiff", units="in", width=9, height=7, res=800)
ggplot(perc_df, aes(fill = class, y = value, x = sensor)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = pal)+
  theme_classic()+
  facet_grid(~metric, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x") + 
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.background = element_rect(fill = "white", color = 'white')) +
  xlab("") + 
  ylab("% Forested Area") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=12), 
        legend.position = "bottom", 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-20,0,0,0))+
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0))
dev.off()






#############################################################################################################
## Figure of the % of each drought category for each month from 1984 - 2020 
#############################################################################################################

spi_new <- data.table(file = list.files("G:/My Drive/Chapter1_ET_Project/Data/SPI/SPI90_MODIS", full.names = T, pattern = ".tif$"), 
                      date = as.Date(paste0(substr(list.files("G:/My Drive/Chapter1_ET_Project/Data/SPI/SPI90_MODIS", full.names = F, pattern = ".tif$"), 1, 6), "01"), format = "%Y%m%d"))
spi_new$mean <- as.numeric(rep(NA, nrow(spi_new)))
spi_new$sd_below <- as.numeric(rep(NA, nrow(spi_new)))
spi_new$sd_above <- as.numeric(rep(NA, nrow(spi_new)))

for(i in 1:nrow(spi_new)){
  r <- as.data.frame(raster(spi_new$file[i]))
  r <- r[complete.cases(r),]
  spi_new$mean[i] <- mean(r, na.rm=T)
  sdd <- sd(r, na.rm=T)
  spi_new$sd_below[i] <- spi_new$mean[i]-sdd
  spi_new$sd_above[i] <- spi_new$mean[i]+sdd
  print(i)
}

spi_new_gs <- spi_new[month(date) > 3 & month(date) < 10,]
tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_2.tiff", units="in", width=8, height=5, res=800)
ggplot(spi_new_gs, aes(x = date, y = mean, ymax = sd_above, ymin = sd_below))+ 
  geom_line(color="black", size = 0.55) + 
  geom_ribbon(alpha=0.5, color = NA, fill = "#800026", size = 1.5) + 
  xlab("Date") + 
  ylab("SPI") +
  theme_classic()+
  theme(legend.title = element_blank(),legend.position = "bottom") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2.5, 2.5)) 
dev.off()



avg_top10_drought_months <- mean(do.call("stack", lapply(spi_new_gs[order(mean),][1:3,]$file, raster)), na.rm=T)

