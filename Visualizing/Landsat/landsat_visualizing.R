######################################################################################
## Prep Landsat ET results for visualization 
######################################################################################
library(data.table)
library(raster)
library(tidyr)
library(ggplot2)
library(parallel)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
cores <- detectCores() - 1
source(paste0(home, "/Visualizing/MODIS/helpers.R"))
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
## Calculate average binned ET anomalies at drought peak across elevation, HAND, and % diff porous 
## for all drought levels (all, moderate, severe, extreme)
#######################################################################################################
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
# h_bins<- seq(0, 600, 10)
# h_names <- h_bins[1:length(h_bins)-1]
# data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)
## WE ACUALLY CHANGED THIS TO TWI BUT THE VAR NAME IS STILL HAND 
h_bins<- seq(5, 35, 0.5)
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



##############################################################################################################
## Calculate average binned overall coupling and sen's slope of rolling coupling across evironmental gradients
##############################################################################################################
agg_DF_XY(data, "elevation_bin", "sensCoupling", "/elevation_sensCoupling.csv")
agg_DF_XY(data, "HAND_bin", "sensCoupling", "/hand_sensCoupling.csv")
agg_DF_XY(data, "pctDiff_R_bin", "sensCoupling", "/pctDiffR_sensCoupling.csv")

agg_DF_XY(data, "elevation_bin", "allCoupling","/elevation_allCoupling.csv")
agg_DF_XY(data, "HAND_bin", "allCoupling", "/hand_allCoupling.csv")
agg_DF_XY(data, "pctDiff_R_bin", "allCoupling", "/pctDiffR_allCoupling.csv")

agg_DF_XY(data, "elevation_bin", "sensAnoms","/elevation_sensAnoms.csv")
agg_DF_XY(data, "HAND_bin", "sensAnoms", "/hand_sensAnoms.csv")
agg_DF_XY(data, "pctDiff_R_bin", "sensAnoms", "/pctDiffR_sensAnoms.csv")


###############################################################################################################
## Calculate a couple of  metrics for the study area section of the paper 
###############################################################################################################

# Calculate the percent of forested area with significant overall coupling and with a significant change in overall coupling
percent_sensCoupling <- (sum(!is.na(data$sensCoupling))/nrow(data))*100  # 1.17%
percent_sensCoupling_positive <- (sum(!is.na(data$sensCoupling) & data$sensCoupling>0)/nrow(data))*100 # 0.3254743
percent_sensCoupling_negative <- (sum(!is.na(data$sensCoupling) & data$sensCoupling<0)/nrow(data))*100 #  0.8484882

percent_overallCoupling <- (sum(!is.na(data$allCoupling))/nrow(data))*100  # 16.06%
percent_overallCoupling_positive <- (sum(!is.na(data$allCoupling) & data$allCoupling>0)/nrow(data))*100 # 14.56378
percent_overallCoupling_negative <- (sum(!is.na(data$allCoupling) & data$allCoupling<0)/nrow(data))*100 # 1.497039

# new 2/27 -- percent and acres that are becoming less constrained 
nrow(data[(data$sensCoupling > 0 & !is.na(data$sensCoupling) & data$allCoupling <0) | 
            (data$sensCoupling < 0 & !is.na(data$sensCoupling) & data$allCoupling > 0),])/nrow(data) * 100 

nrow(data[(data$sensCoupling > 0 & !is.na(data$sensCoupling) & data$allCoupling <0) | 
            (data$sensCoupling < 0 & !is.na(data$sensCoupling) & data$allCoupling > 0),]) * 30 * 30 * 0.000247105



# get the average elevation across the landscape 
landsat_elevation <- raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation30m.tif"))
meanElevation <- mean(values(landsat_elevation), na.rm=T)

# Quantify the percent forested area of the landscape 
landsat_forest <- raster(paste0(home, "/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif"))
v_lf <- values(landsat_forest)
cell0 <- length(v_lf[v_lf == 0 & !is.na(v_lf) ==T])
cell1 <- length(v_lf[v_lf == 1 & !is.na(v_lf) ==T])
avg_forest <- cell1/(cell0 + cell1)


#############################################################################################################
## Find the percent of the drought anomalies that are positive for different levels of drought severity
## And then further break that down across elevation gradient, HAND gradient, and % diff porous 
#############################################################################################################
pct_greater_0(data_long, home, "Landsat")





#######################################################################################################################
# Calculate the percent of forested pixels in each 50m elevation bin and 0.5 TWI bin 
# Calculate the cumulative percentage of the area as elevation/TWI bin increases 
total <- nrow(data)
percent_forest_binned_elevation <- data[order(elevation_bin),.(percent_of_forest = (.N/total)*100), elevation_bin]
percent_forest_binned_elevation$cummulative <- cumsum(percent_forest_binned_elevation$percent_of_forest)
percent_forest_binned_elevation <- percent_forest_binned_elevation[complete.cases(percent_forest_binned_elevation),]
percent_forest_binned_elevation$elevation_bin <- as.numeric(as.character(percent_forest_binned_elevation$elevation_bin))
fwrite(percent_forest_binned_elevation, paste0(home, "/Analysis/outputs/Landsat/forested_area_by_elevation.csv"))

landsat_elev <- fread(paste0(home, "/Analysis/outputs/Landsat/forested_area_by_elevation.csv"))
landsat_elev$sensor <- rep("Landsat", nrow(landsat_elev))
sum(landsat_elev[elevation_bin<500, percent_of_forest])
sum(landsat_elev[elevation_bin>=1000, percent_of_forest])
sum(landsat_elev[elevation_bin>=500 & elevation_bin<1000, percent_of_forest])

# TWI
pfbt <- data[order(HAND_bin), .(percent_forest = (.N/total)*100), HAND_bin]
pfbt$cummulative <- cumsum(pfbt$percent_forest)
pfbt <- pfbt[complete.cases(pfbt),]
pfbt$HAND_bin <- as.numeric(as.character(pfbt$HAND_bin))
fwrite(pfbt, paste0(home, "/Analysis/outputs/Landsat/forested_area_by_TWI.csv"))

landsat_TWI <- fread(paste0(home, "/Analysis/outputs/Landsat/forested_area_by_TWI.csv"))
landsat_TWI$sensor <- rep("Landsat", nrow(landsat_TWI))
sum(landsat_TWI[HAND_bin <=7.5,]$percent_forest)
sum(landsat_TWI[HAND_bin >=10,]$percent_forest)
sum(landsat_TWI[HAND_bin >=10,]$percent_forest)


