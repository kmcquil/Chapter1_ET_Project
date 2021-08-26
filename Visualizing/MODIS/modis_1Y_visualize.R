######################################################################################
## Visualizing results from MODIS portion of analysis 
######################################################################################
library(data.table)
library(raster)
library(rasterVis)
library(pals)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

# bring in the datatable with drought composite info and fc and topo attributes
data <- fread(paste0(home, "/Analysis/outputs/MODIS_1Y/final_drought_attribute_dt.csv"), na.strings = (""))

# factorize the Fc, hand classifications, and elevation classifications 
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)             
data$HAND_class <- as.factor(data$HAND_class)  
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))



#### make maps of the rasters for ET anomalies and ET residuals at each level of drought 
et_tifs <- list.files(paste0(home, "/Analysis/outputs/MODIS_1Y/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs <- et_tifs[grep("ETanom", et_tifs)]
anom_stk <- do.call("stack", lapply(anom_tifs, raster))
residual_tifs <- et_tifs[grep("residual", et_tifs)]
residual_stk <- do.call("stack", lapply(residual_tifs, raster))
plot_names <- c("All", "Moderate", "Severe", "Extreme")


# stretch the anom raster 
qs <- quantile(values(anom_stk), c(0.01, 0.99), na.rm=T)
anom_stk[anom_stk > qs[2]] <- qs[2]
anom_stk[anom_stk < qs[1]] <- qs[1]
anom_plot <- levelplot(anom_stk, names.attr=plot_names, main = "MODIS: ET Drought Composite z-score anomalies (2% stretch)")
diverge0(anom_plot, ramp="RdBu")

ss <- colorRampPalette(rev(kovesi.rainbow(100)))
cw <- colorRampPalette(coolwarm(100))

# stretch the residual raster 
qs <- quantile(values(residual_stk), c(0.02, 0.98), na.rm=T)
residual_stk[residual_stk > qs[2]] <- qs[2]
residual_stk[residual_stk < qs[1]] <- qs[1]
res_plot <- levelplot(residual_stk, names.attr=plot_names, main = "MODIS: ET Drought Composite Residuals (2% stretch)",
                      col.regions = ss)
res_plot



#######################################################################################################
## Now I want to visualize and assess how ET responses vary across continuous variables 
#######################################################################################################

# create 5% bins of percent diffuse porous for Riley fc product 
# min is 0.087 and max is 0.45
# so the bins can start at 5% and end at 85% 
pd_bins <- seq(0.05, 0.85, 0.05)
pd_names <- pd_bins[1:length(pd_bins)-1]*100
data$pctDiff_R_bin <- cut(data$pctDiff_R, breaks = pd_bins, labels = pd_names)
data$pctDiff_W_bin <- cut(data$pctDiff_W, breaks = pd_bins, labels = pd_names)


# create 50m bins of elevation 
# min is 218 and max is 1950
# so the bins can start at 200 and go to 1925
e_bins <- seq(200, 1950, 50)
e_names <- e_bins[1:length(e_bins)-1]
data$elevation_bin <- cut(data$elevation, breaks = e_bins, labels = e_names)


# create 10m bins of HAND 
# min is 1.9 and max is 261
# so bins start at 0 and go to 270 
h_bins<- seq(0, 270, 10)
h_names <- h_bins[1:length(h_bins)-1]
data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)

data_anom <- data[,c(1,2,3,4,5,10,11,12,13,14,15,16,17,18,20,21,22, 23)]
data_anom_long <- gather(data_anom, DroughtSeverity, ET_Anomaly, c(ET_all:ET_extreme, avgNoDroughtETanom), factor_key = F)
data_anom_long$DroughtSeverity <- ifelse(data_anom_long$DroughtSeverity == "ET_all", "All", 
                                         ifelse(data_anom_long$DroughtSeverity == "ET_moderate", "Moderate", 
                                                ifelse(data_anom_long$DroughtSeverity == "ET_severe", "Severe", 
                                                       ifelse(data_anom_long$DroughtSeverity == "ET_extreme", "Extreme",
                                                              ifelse(data_anom_long$DroughtSeverity == "avgNoDroughtETanom", "None", NA)))))


data_res <- data[, c(1,6,7,8,9,10,11,12,13, 14, 15,16,17,19,20,21,22, 23)]
data_res_long <- gather(data_res, DroughtSeverity, ET_Residual, c(Res_all:Res_extreme, avgNoDroughtETres), factor_key = F)
data_res_long$DroughtSeverity <- ifelse(data_res_long$DroughtSeverity == "Res_all", "All", 
                                        ifelse(data_res_long$DroughtSeverity == "Res_moderate", "Moderate", 
                                               ifelse(data_res_long$DroughtSeverity == "Res_severe", "Severe", 
                                                      ifelse(data_res_long$DroughtSeverity == "Res_extreme", "Extreme", 
                                                             ifelse(data_res_long$DroughtSeverity == "avgNoDroughtETres", "None", NA)))))


data_long <- merge(data_anom_long, data_res_long, by = colnames(data_anom_long)[1:14], all.x=T, all.y=T)
data_long <- setDT(data_long)
data_long$DroughtSeverity <- ordered(data_long$DroughtSeverity, levels = c("None", "Moderate", "Severe", "Extreme", "All"))

data_long <- data_long[!(is.na(ET_Anomaly) & is.na(ET_Residual)),]
data_long$pctDiff_R_bin <- as.numeric(as.character(data_long$pctDiff_R_bin))
data_long$pctDiff_W_bin <- as.numeric(as.character(data_long$pctDiff_W_bin))
data_long$elevation_bin <- as.numeric(as.character(data_long$elevation_bin))
data_long$HAND_bin <- as.numeric(as.character(data_long$HAND_bin))



# aggregate and find the mean and sd of ET response (anomaly and residual) grouped by elevation bin and drought type 
data_long_sub <- data_long[!DroughtSeverity == "All",]
contPlots(data_long_sub, "elevation_bin", "DroughtSeverity", "Binned Elevation")

# now group by HAND and drought type 
contPlots(data_long_sub, "HAND_bin", "DroughtSeverity", "Binned HAND")

# now group by percent diffuse bin and drought type 
contPlots(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "Binned % Diffuse Porous")


# now subset first by drought type, group by % diffuse porous binned and elevation 
data_long_all <- data_long[DroughtSeverity == "All", ]
contPlots(data_long_all, "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")

data_long_mod <- data_long[DroughtSeverity == "Moderate", ]
contPlots(data_long_mod, "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")

data_long_severe <- data_long[DroughtSeverity == "Severe", ]
contPlots(data_long_severe, "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")

data_long_extreme <- data_long[DroughtSeverity == "Extreme", ]
contPlots(data_long_extreme, "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")

contPlots(data_long_all, "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots(data_long_mod, "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots(data_long_severe, "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots(data_long_extreme, "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")



###################################################################################################################
## Visualize results from the rolling coupling of monthly GS SPI and ET anoms and the trends of the coupling on a 
##  pixel and whole region basis 
###################################################################################################################

# bring this in to make sure everything was masked to only forested pixels 
forest_mask <- as.data.frame(raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif"))
forest_mask <- cbind(seq(1, nrow(forest_mask)), forest_mask)
colnames(forest_mask) <- c("cellnum", "fmask")
forest_mask1 <- forest_mask[forest_mask$fmask== 1 & !is.na(forest_mask$fmask),]

# bring in the pixel wise sens slope of the 5 year rolling coupling and mask out sens slope for pixels that are not significant 
coupling_slope <- fread(paste0(home, "/Analysis/outputs/MODIS_1Y/rollingCoupling/coupling_sensslope.csv"))
coupling_slope <- coupling_slope[coupling_slope$cellnum %in% forest_mask1$cellnum,]
coupling_slope_sig <- coupling_slope[coupling_slope$pvalue < 0.05,]

# convert the significant sens slope to a raster
sig_slope <- DTtoRast(coupling_slope_sig[,c(1,2)], anom_stk[[1]], 
                      paste0(home, "/Analysis/outputs/MODIS_1Y/rollingCoupling/sig_slope.tif"), open = T)
sig_slope_plot <- levelplot(sig_slope, margin=F, main = "Significant slopes", par.settings=list(panel.background=list(col="black")))
diverge0(sig_slope_plot, ramp="RdBu")

# calculate the % area that is significantly changing in coupling 
percent_area_changing <- (nrow(coupling_slope_sig)/nrow(forest_mask1))*100

# find that 5% area changing 


# bring in the 5 year rolling coupling and the overall coupling 
rollingCoupling <- fread(paste0(home, "/Analysis/outputs/MODIS_1Y/rollingCoupling/rollingCoupling.csv"))
rollingCoupling <- rollingCoupling[rollingCoupling$cellnum %in% forest_mask1$cellnum,]

# convert the overall coupling to a raster 
DTtoRast(rollingCoupling[,c("cellnum", "all")], anom_stk[[1]], 
         paste0(home, "/Analysis/outputs/MODIS_1Y/rollingCoupling/overallCoupling.tif"), open = T)
coupling_rast <- raster(paste0(home, "/Analysis/outputs/MODIS_1Y/rollingCoupling/overallCoupling.tif"))
coupling_plot <- levelplot(coupling_rast, margin=F, main = "R (ET ~ SPI)", par.settings=list(panel.background=list(col="black")))
diverge0(coupling_plot, ramp="RdBu")

# mostly positively coupled 

# merge the overall coupling and the sens slope (only significant) of the rolling coupling with the dataframe so i can look at diferences 
# along elevation, hand, and diffuse porous gradients 
data <- merge(data, coupling_slope_sig[,c("cellnum", "senSlope")], by = "cellnum", all.x=T)
data <- merge(data, rollingCoupling[, c("cellnum", "all")], by = "cellnum", all.x=T)


# merge the sens slope of ET anomalies 
sensAnom <- as.data.frame(raster(paste0(home, "/Analysis/outputs/MODIS/sig_drought_trend.tif")))
sensAnom <- cbind(cellnum = seq(1, nrow(sensAnom)), sensAnom = sensAnom$sig_drought_trend)
data <- merge(data, sensAnom, by = "cellnum", all.x=T)



# get the percent area significantly coupled at each time step 
# percent greater than 0 
g0 <- apply(rollingCoupling[,2:(ncol(rollingCoupling)-1)], MARGIN=2, FUN = function(x) {  sum(x>0, na.rm=T) })
g0f <- g0/nrow(forest_mask1)
nams <- colnames(rollingCoupling[,2:(ncol(rollingCoupling)-1)])
date_labels <- as.Date(paste0(substr(nams, 2, 5),  substr(nams, 6, 7), "01"), format = "%Y%m%d")

# percent less than 0
l0 <- apply(rollingCoupling[,2:(ncol(rollingCoupling)-1)], MARGIN=2, FUN = function(x) {  sum(x<0, na.rm=T) })
l0f <- l0/nrow(forest_mask1)
plot(date_labels, l0f, 'l')

# non-significantly coupled 
n0f <- 1-(g0f + l0f)

# plot all together 
areas <- data.table(Date = date_labels, greater = g0f, less = l0f, no = n0f)
areas_long <- gather(areas, direction, percent_area, greater:no)
areas_long$direction <- ifelse(areas_long$direction == "greater", "R > 0", 
                               ifelse(areas_long$direction == "less", "R<0", 
                                      ifelse(areas_long$direction == "no", "Not Significant", NA)))
areas_long$percent_area <- areas_long$percent_area*100
ggplot(data = areas_long, aes(x = Date, y = percent_area, color = direction)) + 
  geom_line(size = 1.7) + 
  xlab("Date") + 
  ylab("% of Forested Area") + 
  scale_color_manual(values = viridis(4))+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

