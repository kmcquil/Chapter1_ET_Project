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


pd_bins <- seq(0, 1, 0.05)
pd_names <- pd_bins[1:length(pd_bins)-1]*100
data$pctDiff_R_bin <- cut(data$pctDiff_R, breaks = pd_bins, labels = pd_names)
data$pctDiff_W_bin <- cut(data$pctDiff_W, breaks = pd_bins, labels = pd_names)


e_bins <- c(199, seq(250, 2050, 50))
e_names <- e_bins[1:length(e_bins)-1]
data$elevation_bin <- cut(data$elevation, breaks = e_bins, labels = e_names)

# create 10m bins of HAND
# min is 0 and max is 592
# so bins start at 0 and go to 600
## Actually sub in TWI 
twi=values(raster(paste0(home, "/Data/Topography/usgsNED_elevation/Landsat_TWI/TWI_resampled.tif")))
twi <- data.table(twi=twi, 
                  cellnum = seq(1, length(twi), 1))
data[,HAND:=as.numeric(rep(NA, nrow(data)))]
setkey(data, cellnum)
setkey(twi, cellnum)
data[twi, HAND:=twi]
h_bins<- seq(5, 35, 1)
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
#data_long$HAND_bin <- as.numeric(as.character(data_long$HAND_bin))
data_long_sub <- data_long[!DroughtSeverity == "All",]


source(paste0(home, "/Visualizing/MODIS/helpers.R"))

agg_DF(data_long_sub, "elevation_bin", "DroughtSeverity", "/TWI/elevation_DS.csv")
agg_DF(data_long_sub, "HAND_bin", "DroughtSeverity", "/TWI/HAND_DS.csv")
agg_DF(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "/TWI/pctDiffR_DS.csv")


agg_DF_XY(data, "elevation_bin", "sensCoupling", "/TWI/elevation_sensCoupling.csv")
agg_DF_XY(data, "HAND_bin", "sensCoupling", "/TWI/hand_sensCoupling.csv")
agg_DF_XY(data, "pctDiff_R_bin", "sensCoupling", "/TWI/pctDiffR_sensCoupling.csv")

agg_DF_XY(data, "elevation_bin", "allCoupling","/TWI/elevation_allCoupling.csv")
agg_DF_XY(data, "HAND_bin", "allCoupling", "/TWI/hand_allCoupling.csv")
agg_DF_XY(data, "pctDiff_R_bin", "allCoupling", "/TWI/pctDiffR_allCoupling.csv")




#########################################################################################33
## Find whether diffuse porous BA differs across HAND bins 
rm(data_long_sub)
h_bins <- c(5,10,15,30)
h_names1 <- h_bins[1:length(h_bins)-1]
h_names2 <- h_bins[2:length(h_bins)]
h_names <- paste0(h_names1, "-",h_names2)
data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)
data_long$HAND_bin <- cut(data_long$HAND, breaks = h_bins, labels = h_names)
data_long_sub <- data_long[!DroughtSeverity == "All",]


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

dir.create("/share/klmarti3/kmcquil/Chapter1_ET_Project/Visualizing/Landsat/TWI_bin")
agg_DF_new(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "HAND_bin", "/Landsat/TWI_bin","pctDiffR_DS1.csv")



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

agg_DF_XY_new(data, "pctDiff_R_bin", "HAND_bin","sensCoupling", "Landsat/TWI_bin","/pctDiffR_sensCoupling.csv")
agg_DF_XY_new(data, "pctDiff_R_bin", "HAND_bin","allCoupling", "Landsat/TWI_bin","/pctDiffR_allCoupling.csv")
