######################################################################################
## Visualizing results from Landsat portion of analysis 
######################################################################################
library(data.table)
library(raster)
library(tidyr)
library(ggplot2)
library(parallel)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
cores <- detectCores() - 1
data <- fread(paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"), na.strings = (""), nThread = cores)

# factorize the Fc, hand classifications, and elevation classifications
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)
data$HAND_class <- as.factor(data$HAND_class)
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))

# pearson correlation between elevation, HAND, pctDiff_R
preds <- data[,.(elevation, HAND, pctDiff_R)]
landsat_corr <- cor(preds, use="pairwise.complete.obs")



#######################################################################################################
## Visualize and assess how ET responses vary across continuous variables
# create 5% bins of percent diffuse porous for Riley fc product
# min is 0.087 and max is 0.95
# so the bins can start at 5% and end at 85%
pd_bins <- seq(0, 1, 0.05)
pd_names <- pd_bins[1:length(pd_bins)-1]*100
data$pctDiff_R_bin <- cut(data$pctDiff_R, breaks = pd_bins, labels = pd_names)
data$pctDiff_W_bin <- cut(data$pctDiff_W, breaks = pd_bins, labels = pd_names)
# create 50m bins of elevation
# min is 199 and max is 2034
# so the bins can start at 200 and go to 1925
e_bins <- c(199, seq(250, 2050, 50))
e_names <- e_bins[1:length(e_bins)-1]
data$elevation_bin <- cut(data$elevation, breaks = e_bins, labels = e_names)
# create 10m bins of HAND
# min is 0 and max is 592
# so bins start at 0 and go to 600
h_bins<- seq(0, 600, 10)
h_names <- h_bins[1:length(h_bins)-1]
data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)

# reshape so that it is long ways and there is a column that indicates drought severity or no drought
# make separate long dfs for anoms and residuals
data_anom <- data[,c(seq(1,5), seq(10, 27))]
data_anom_long <- gather(data_anom, DroughtSeverity, ET_Anomaly, c(ET_all:ET_extreme, avgNoDroughtETanom), factor_key = F)
data_anom_long$DroughtSeverity <- ifelse(data_anom_long$DroughtSeverity == "ET_all", "All",
                                         ifelse(data_anom_long$DroughtSeverity == "ET_moderate", "Moderate",
                                                ifelse(data_anom_long$DroughtSeverity == "ET_severe", "Severe",
                                                       ifelse(data_anom_long$DroughtSeverity == "ET_extreme", "Extreme",
                                                              ifelse(data_anom_long$DroughtSeverity == "avgNoDroughtETanom", "None", NA)))))
data_res <- data[, c(1,seq(6,27))]
data_res_long <- gather(data_res, DroughtSeverity, ET_Residual, c(Res_all:Res_extreme, avgNoDroughtETres), factor_key = F)
data_res_long$DroughtSeverity <- ifelse(data_res_long$DroughtSeverity == "Res_all", "All",
                                        ifelse(data_res_long$DroughtSeverity == "Res_moderate", "Moderate",
                                               ifelse(data_res_long$DroughtSeverity == "Res_severe", "Severe",
                                                      ifelse(data_res_long$DroughtSeverity == "Res_extreme", "Extreme",
                                                             ifelse(data_res_long$DroughtSeverity == "avgNoDroughtETres", "None", NA)))))
          
data_anom_long <- setDT(data_anom_long)
data_res_long <- setDT(data_res_long)
setkeyv(data_anom_long, colnames(data_anom_long)[1:19])
setkeyv(data_res_long, colnames(data_res_long)[1:19])
data_long <- merge(data_anom_long, data_res_long, all.=T, all.y=T)
data_long <- setDT(data_long)
data_long$DroughtSeverity <- ordered(data_long$DroughtSeverity, levels = c("None", "Moderate", "Severe", "Extreme", "All"))
data_long <- data_long[!(is.na(ET_Anomaly) & is.na(ET_Residual)),]
data_long$pctDiff_R_bin <- as.numeric(as.character(data_long$pctDiff_R_bin))
data_long$pctDiff_W_bin <- as.numeric(as.character(data_long$pctDiff_W_bin))
data_long$elevation_bin <- as.numeric(as.character(data_long$elevation_bin))
data_long$HAND_bin <- as.numeric(as.character(data_long$HAND_bin))



data_long_sub <- data_long[!DroughtSeverity == "All",]
agg_DF(data_long_sub, "elevation_bin", "DroughtSeverity", "/elevation_DS.csv")
agg_DF(data_long_sub, "HAND_bin", "DroughtSeverity", "/HAND_DS.csv")
agg_DF(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "/pctDiffR_DS.csv")


data_long_all <- data_long[DroughtSeverity == "All", ]
agg_DF(data_long_all, "pctDiff_R_bin", "Elev_class", "/pctDiffR_ElevClass_all.csv") 

data_long_mod <- data_long[DroughtSeverity == "Moderate", ]
agg_DF(data_long_mod, "pctDiff_R_bin", "Elev_class", "/pctDiffR_ElevClass_moderate.csv") 

data_long_severe <- data_long[DroughtSeverity == "Severe", ]
agg_DF(data_long_severe, "pctDiff_R_bin", "Elev_class", "/pctDiffR_ElevClass_severe.csv") 

data_long_extreme <- data_long[DroughtSeverity == "Extreme", ]
agg_DF(data_long_extreme, "pctDiff_R_bin", "Elev_class", "/pctDiffR_ElevClass_extreme.csv")


### now subset first by drought type, group by % diffuse porous binned and HAND
agg_DF(data_long_all, "pctDiff_R_bin", "HAND_class", "/pctDiffR_HANDClass_all.csv")
agg_DF(data_long_mod, "pctDiff_R_bin", "HAND_class", "/pctDiffR_HANDClass_moderate.csv")
agg_DF(data_long_severe, "pctDiff_R_bin", "HAND_class", "/pctDiffR_HANDClass_severe.csv")
agg_DF(data_long_extreme, "pctDiff_R_bin", "HAND_class", "/pctDiffR_HANDClass_extreme.csv")



# just use regular data 
agg_DF_XY(data, "elevation_bin", "sensCoupling", "/elevation_sensCoupling.csv")
agg_DF_XY(data, "HAND_bin", "sensCoupling", "/hand_sensCoupling.csv")

agg_DF_XY(data, "pctDiff_R_bin", "sensCoupling", "/pctDiffR_sensCoupling.csv")
agg_DF_XY(data, "elevation_bin", "allCoupling","/elevation_allCoupling.csv")
agg_DF_XY(data, "HAND_bin", "allCoupling", "/hand_allCoupling.csv")
agg_DF_XY(data, "pctDiff_R_bin", "allCoupling", "/pctDiffR_allCoupling.csv")

agg_DF_XY(data, "elevation_bin", "sensAnoms","/elevation_sensAnoms.csv")
agg_DF_XY(data, "HAND_bin", "sensAnoms", "/hand_sensAnoms.csv")
agg_DF_XY(data, "pctDiff_R_bin", "sensAnoms", "/pctDiffR_sensAnoms.csv")


data$elevation_bin <- as.numeric(as.character(data$elevation_bin))
data$pctDiff_R_bin <- as.numeric(as.character(data$pctDiff_R_bin))
agg_DF_XYZ(data, "pctDiff_R_bin", "allCoupling", "Elev_class", "/pctDiffR_allCoupling_elevClass.csv")
agg_DF_XYZ(data, "pctDiff_R_bin", "sensCoupling", "Elev_class", "/pctDiffR_sensCoupling_elevClass.csv")
agg_DF_XYZ(data, "pctDiff_R_bin", "sensAnoms", "Elev_class", "/pctDiffR_sensAnoms_elevClass.csv")
agg_DF_XYZ(data, "pctDiff_R_bin", "allCoupling", "HAND_class", "/pctDiffR_allCoupling_HANDClass.csv")
agg_DF_XYZ(data, "pctDiff_R_bin", "sensCoupling", "HAND_class", "/pctDiffR_sensCoupling_HANDClass.csv")
agg_DF_XYZ(data, "pctDiff_R_bin", "sensAnoms", "HAND_class", "/pctDiffR_sensAnoms_HANDClass.csv")








#########################################################################################################
## Make actual figures using data summarized on hpc 
#########################################################################################################
library(rasterVis)
library(pals)
library(RColorBrewer)
library(trend)
library(data.table)
library(ggplot2)
home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

#### make maps of the rasters for ET anomalies and ET residuals at each level of drought 
et_tifs <- list.files(paste0(home, "/Analysis/outputs/Landsat/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs <- et_tifs[grep("ETanom", et_tifs)]
anom_stk <- do.call("stack", lapply(anom_tifs, raster))
residual_tifs <- et_tifs[grep("residual", et_tifs)]
residual_stk <- do.call("stack", lapply(residual_tifs, raster))
plot_names <- c("All", "Moderate", "Severe", "Extreme")

forest_mask_landsat <- raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
forest_mask_landsat[forest_mask_landsat == 0] <- NA
plot_names <- c("All", "Moderate", "Severe", "Extreme")

anom_stk <- mask(anom_stk, forest_mask_landsat)
qs <- quantile(values(anom_stk), c(0.01, 0.99), na.rm=T)
anom_stk[anom_stk > qs[2]] <- qs[2]
anom_stk[anom_stk < qs[1]] <- qs[1]
anom_plot <- levelplot(anom_stk, names.attr=plot_names, main = "Landsat: ET Drought Composite z-score anomalies (2% stretch)")
diverge0(anom_plot, ramp="RdYlBu")
diverge0(anom_plot, ramp="RdBu")

residual_stk <- mask(residual_stk, forest_mask_landsat)
qs <- quantile(values(residual_stk), c(0.01, 0.99), na.rm=T)
residual_stk[residual_stk > qs[2]] <- qs[2]
residual_stk[residual_stk < qs[1]] <- qs[1]
res_plot <- levelplot(residual_stk, names.attr=plot_names, main = "Landsat: ET Drought Composite residuals (2% stretch)", 
                      par.settings=list(panel.background=list(col="black")) )
diverge0(res_plot, ramp = "RdBu")






## plots of anomalies and residuals over gradients 
dir <- "G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat"
contPlots_1(paste0(dir, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Binned Elevation")
contPlots_1(paste0(dir, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "Binned HAND")


contPlots_1(paste0(dir, "/pctDiffR_ElevClass_all.csv"), "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/pctDiffR_ElevClass_moderate.csv"), "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/pctDiffR_ElevClass_severe.csv"), "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous") 
contPlots_1(paste0(dir, "/pctDiffR_ElevClass_extreme.csv"), "pctDiff_R_bin", "Elev_class", "Binned % Diffuse Porous") 

contPlots_1(paste0(dir, "/pctDiffR_HANDClass_all.csv"), "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/pctDiffR_HANDClass_moderate.csv"), "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/pctDiffR_HANDClass_severe.csv"), "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")
contPlots_1(paste0(dir, "/pctDiffR_HANDClass_extreme.csv"), "pctDiff_R_bin", "HAND_class", "Binned % Diffuse Porous")


# trend in sen slope with elevation, diffuse porous, and hand 
contPlots_single_1(paste0(dir, "/elevation_sensCoupling.csv"), "elevation_bin", "senSlope", xlab = "Binned Elevation", "Sen's Slope of Coupling")
contPlots_single_1(paste0(dir, "/HAND_sensCoupling.csv"), "HAND_bin", "senSlope", xlab = "Binned HAND", "Sen's Slope of Coupling")
contPlots_single_1(paste0(dir, "/pctDiffR_sensCoupling.csv"), "pctDiff_R_bin", "senSlope", xlab = "Binned % Diffuse Porous", "Sen's Slope of coupling")

# trend in the average coupling with elevation, diffuse porous, and hand 
contPlots_single_1(paste0(dir, "/elevation_allCoupling.csv"), "elevation_bin", "all", xlab = "Binned Elevation", "R (SPI ~ ET anom)")
contPlots_single_1(paste0(dir, "/HAND_allCoupling.csv"), "HAND_bin", "all", xlab = "Binned HAND", "R (SPI ~ ET anom)")
contPlots_single_1(paste0(dir, "/pctDiffR_allCoupling.csv"), "pctDiff_R_bin", "all", xlab = "Binned % Diffuse Porous", "R (SPI ~ ET anom)")

# trend in sens slope of anomalies with diffuse porous, elevation, and hand 
contPlots_single_1(paste0(dir, "/elevation_sensAnoms.csv"), "elevation_bin", "senSlope", xlab = "Binned Elevation", "Sen's Slope of Anomalies")
contPlots_single_1(paste0(dir, "/HAND_sensAnoms.csv"), "HAND_bin", "senSlope", xlab = "Binned HAND", "Sen's Slope of Anomalies")
contPlots_single_1(paste0(dir, "/pctDiffR_sensAnoms.csv"), "pctDiff_R_bin", "senSlope", xlab = "Binned % Diffuse Porous", "Sen's Slope of Anomalies")





# now break down further by elevation class 
contPlots_rolling_1(paste0(dir, "/pctDiffR_allCoupling_elevClass.csv"), "pctDiff_R_bin", "all", "Elev_class", "Binned % Diffuse Porous", "R (ET~SPI)")
contPlots_rolling_1(paste0(dir, "/pctDiffR_sensCoupling_elevClass.csv"), "pctDiff_R_bin", "senSlope", "Elev_class", "Binned % Diffuse Porous", "Sen's Slope of Coupling")
contPlots_rolling_1(paste0(dir, "/pctDiffR_sensAnoms_elevClass.csv"), "pctDiff_R_bin", "senSlope", "Elev_class", "Binned % Diffuse Porous", "Sen's Slope of Anomalies")


# now by HAND 
contPlots_rolling_1(paste0(dir, "/pctDiffR_allCoupling_HANDClass.csv"), "pctDiff_R_bin", "all", "HAND_class", "Binned % Diffuse Porous", "R (ET~SPI)")
contPlots_rolling_1(paste0(dir, "/pctDiffR_sensCoupling_HANDClass.csv"), "pctDiff_R_bin", "senSlope", "HAND_class", "Binned % Diffuse Porous", "Sen's Slope of coupling")
contPlots_rolling_1(paste0(dir, "/pctDiffR_sensAnoms_HANDClass.csv"), "pctDiff_R_bin", "senSlope", "HAND_class", "Binned % Diffuse Porous", "Sen's Slope of Anomalies")






# I want to calculate the percent of forested area that the rolling coupling is significnatly changing 
percent_sensCoupling <- (sum(!is.na(data$sensCoupling))/nrow(data))*100  # 1.17%
percent_overallCoupling <- (sum(!is.na(data$allCoupling))/nrow(data))*100  # 16.06%
percent_sensAnoms <- (sum(!is.na(data$sensAnoms))/nrow(data))*100 # 3% 





# get the average elevation 
landsat_elevation <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/elevation30m.tif")
meanElevation <- mean(values(landsat_elevation), na.rm=T)

landsat_forest <- raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
v_lf <- values(landsat_forest)
cell0 <- length(v_lf[v_lf == 0 & !is.na(v_lf) ==T])
cell1 <- length(v_lf[v_lf == 1 & !is.na(v_lf) ==T])

avg_forest <- cell1/(cell0 + cell1)






#################################################################################################
## Find the percent of the drought anomalies that are positive for different levels of drought severity
## And then further break that down across elevation gradient, HAND gradient, and % diff porous 
##################################################################################################

pct_greater_0(data_long, home, "Landsat")








