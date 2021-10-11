#################################################################################################
# MODIS 
# Look at the variation in ETdp, R(SPI-ET), and dR(SPI-ET) across % diffuse porous BA gradients 
# binned by HAND and elevation 

library(data.table)
library(raster)
library(rasterVis)
library(pals)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(parallel)

home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

data <- fread(paste0(home, "/Analysis/outputs/MODIS/final_drought_attribute_dt.csv"), na.strings = (""))

# factorize the Fc, hand classifications, and elevation classifications 
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)             
data$HAND_class <- as.factor(data$HAND_class)  
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))

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
e_bins <- seq(200, 2000, 600)
e_names1 <- e_bins[1:length(e_bins)-1]
e_names2 <- e_bins[2:length(e_bins)]
e_names <- paste0(e_names1, "-",e_names2)
data$elevation_bin <- cut(data$elevation, breaks = e_bins, labels = e_names)


# create 10m bins of HAND 
# min is 1.9 and max is 261
# so bins start at 0 and go to 270 
h_bins<- seq(0, 270, 10)
h_names <- h_bins[1:length(h_bins)-1]
data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)

# reshape so that it is long ways and there is a column that indicates drought severity or no drought
# make separate long dfs for anoms and residuals 
data_anom <- data[,c(1,2,3,4,5,10,11,12,13,14,15,16,17,18,20,21, 22, 23)]
data_anom_long <- gather(data_anom, DroughtSeverity, ET_Anomaly, c(ET_all:ET_extreme, avgNoDroughtETanom), factor_key = F)
data_anom_long$DroughtSeverity <- ifelse(data_anom_long$DroughtSeverity == "ET_all", "All", 
                                         ifelse(data_anom_long$DroughtSeverity == "ET_moderate", "Moderate", 
                                                ifelse(data_anom_long$DroughtSeverity == "ET_severe", "Severe", 
                                                       ifelse(data_anom_long$DroughtSeverity == "ET_extreme", "Extreme",
                                                              ifelse(data_anom_long$DroughtSeverity == "avgNoDroughtETanom", "None", NA)))))


data_res <- data[, c(1,6,7,8,9,10,11,12,13, 14, 15,16,17,19,20,21, 22, 23)]
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
#data_long$elevation_bin <- as.numeric(as.character(data_long$elevation_bin))
data_long$HAND_bin <- as.numeric(as.character(data_long$HAND_bin))
data_long_sub <- data_long[!DroughtSeverity == "All",]

#######################################################################################################
# get the mean and sd of ETdp grouped by drought severity, % diff bin + elevation bin 
agg_DF_new <- function(DF, X, Y, Z, folder, out){
  agg <- DF[, .(meanAnom = mean(ET_Anomaly, na.rm=T), 
                meanRes = mean(ET_Residual, na.rm=T), 
                sdAnom = sd(ET_Anomaly, na.rm=T), 
                sdRes = sd(ET_Residual, na.rm=T), 
                count = .N), by=c(X, Z, Y)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,] # only keep bins with at least 100 observations 
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  agg$res_plus_sd <- agg$meanRes + agg$sdRes
  agg$res_minus_sd <- agg$meanRes - agg$sdRes
  
  
  fwrite(agg, paste0(home, "/Visualizing/",folder,"/", out))
}

home <- "G:/My Drive/Chapter1_ET_Project"
agg_DF_new(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "elevation_bin", "/MODIS/trial","pctDiffR_DS1.csv")
agg_DF_new(data_long_sub, "HAND_bin", "DroughtSeverity", "elevation_bin", "/MODIS/trial","HAND_DS1.csv")



############################################################################################################
# get the mean and sd of R and dR grouped by % diff bin and HAND + elevation bin 
forest_mask <- as.data.frame(raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif"))
forest_mask <- cbind(seq(1, nrow(forest_mask)), forest_mask)
colnames(forest_mask) <- c("cellnum", "fmask")
forest_mask1 <- forest_mask[forest_mask$fmask== 1 & !is.na(forest_mask$fmask),]
coupling_slope <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/coupling_sensslope.csv"))
coupling_slope <- coupling_slope[coupling_slope$cellnum %in% forest_mask1$cellnum,]
coupling_slope_sig <- coupling_slope[coupling_slope$pvalue < 0.05,]
rollingCoupling <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/rollingCoupling.csv"))
rollingCoupling <- rollingCoupling[rollingCoupling$cellnum %in% forest_mask1$cellnum,]
data <- merge(data, coupling_slope_sig[,c("cellnum", "senSlope")], by = "cellnum", all.x=T)
data <- merge(data, rollingCoupling[, c("cellnum", "all")], by = "cellnum", all.x=T)

agg_DF_XY_new <- function(DF, X, Z, Y, folder, out){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X, Z)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,] # don't keep groups with less than 100 observations
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  agg <- agg[order(agg[[1]]),]
  fwrite(agg, paste0(home, "/Visualizing/",folder, out))
}

agg_DF_XY_new(data, "pctDiff_R_bin", "elevation_bin","senSlope", "MODIS/trial","/pctDiffR_sensCoupling1.csv")
agg_DF_XY_new(data, "pctDiff_R_bin", "elevation_bin","all", "MODIS/trial","/pctDiffR_allCoupling1.csv")

agg_DF_XY_new(data, "HAND_bin", "elevation_bin","senSlope", "MODIS/trial","/HAND_sensCoupling1.csv")
agg_DF_XY_new(data, "HAND_bin", "elevation_bin","all", "MODIS/trial","/HAND_allCoupling1.csv")

