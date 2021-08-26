############################################################################################
## Validate the MODIS and Landsat ET data 
############################################################################################
library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(doParallel)
library(data.table)
library(Metrics)
library(ggplot2)

home <- "G:/My Drive/Chapter1_ET_Project"

# bring in data frame with the coweeta daily ET calculated from the flux tower 
ec_et <- fread("G:/My Drive/Dissertation/ameriflux/et_mod_obs_df.csv")
ec_et <- ec_et[V1 >= as.Date("2011-01-01") & V1 < as.Date("2016-01-01")]

# bring in the coweeta flux tower point 
cwt <- readOGR(paste0(home, "/Data/coweeta/cwt_flx_twr.shp"))

# bring in the subset ET 8-day data 
tifs <- list.files(paste0(home, "/Data_Testing/MODIS_ET/clean"), full.names=T, pattern = ".tif$")
tifs_short <- list.files(paste0(home, "/Data_Testing/MODIS_ET/clean"), full.names=F, pattern = ".tif$")

# make a data.table of the file name, start date, end date, and empty numeric column of the MODIS ET and the flux tower ET 
et_df <- data.table(file = tifs, 
                    start_date = as.Date(substr(tifs_short, 4, 11), "%Y%m%d"), 
                    end_date = as.Date(rep(NA, length(tifs))), 
                    MODIS_ET = as.numeric(rep(NA, length(tifs))), 
                    OBS_ET = as.numeric(rep(NA, length(tifs))),
                    countNA = as.numeric(rep(NA, length(tifs))))

for(i in 1:nrow(et_df)-1){
  et_df$end_date[i] <- et_df$start_date[i+1]-1
}

et_df <- et_df[start_date >= as.Date("2011-01-01") & end_date < as.Date("2016-01-01")]

# the flux tower sits towards teh upper right of teh MODIS pixel. 
# add a 300m buffer around the flux tower and use that to extract the average MODIS value 
cwt_buff <- buffer(cwt, 300)

for(i in 1:nrow(et_df)){
  et_df$MODIS_ET[i] <- raster::extract(raster(et_df$file[i]), cwt_buff, fun = mean, na.rm=T, df=T)[1,2]
  
  et_df$OBS_ET[i] <- sum(ec_et[V1 >= et_df$start_date[i] & V1 <= et_df$end_date[i]]$ET_obs_coweeta, na.rm=T)
  et_df$countNA[i] <- sum(is.na(ec_et[V1 >= et_df$start_date[i] & V1 <= et_df$end_date[i]]$ET_obs_coweeta))
  print(i)
}



# calculate R2, slope, MAE, RMSE, RE (%) where RE = MAE/mean(observed) and return as a data.table
metrics_fun <- function(observed, predicted, time){
  n = length(observed)
  et_lm <- lm(predicted ~ observed)
  slope <- et_lm$coefficients[2]
  r2 <- summary(et_lm)$r.squared
  RMSE <- rmse(observed, predicted) # remember this is for 8 days 
  MAE <- mae(observed, predicted)
  RE <- (MAE/mean(observed, na.rm=T))*100
  
  met_dt <- data.table(TimeStep = time, Count = n, RMSE = RMSE, MAE = MAE, RE = RE, Slope = slope, R2 = r2)
  return(met_dt)
}

Mod8day <- metrics_fun(et_df$OBS_ET, et_df$MODIS_ET, "8Day")



## now look get monthly averages and do the metrics again 
monthly_et_df <- et_df[, .(MODIS_ET_avg = mean(MODIS_ET), OBS_ET_avg = mean(OBS_ET)), .(year(start_date), month(start_date))]
ModMonthly <- metrics_fun(monthly_et_df$OBS_ET_avg, monthly_et_df$MODIS_ET_avg, "Monthly")



par(mfrow=c(1,4), omi=c(0,0.25,0,0.25), mar = c(5.1,6.1,4.1,1.1))

plot(et_df$OBS_ET, et_df$MODIS_ET, 
     pch = 16, cex = 1, col = "#313695",main="", 
     xlab = "Obs. ET (mm/8day)", 
     ylab = "MODIS ET (mm/8day)", 
     xlim = c(0,48), 
     ylim = c(0,48), 
     cex.lab = 1.65, 
     font = 1,
     font.lab = 2,
     cex.axis=1.5)
abline(lm(et_df$MODIS_ET ~ et_df$OBS_ET), col = "#313695", lwd = 3)
abline(0, 1, col = "#d73027", lty=2, lwd=3)
legend(-3, 49.9, legend=c("Fitted (Pred ~ Obs)", "1:1"),
       col=c("#313695", "#d73027"), lty=1:3, cex=1.6, bty='n', x.intersp =0.2, seg.len=0.75, text.font = 2)


plot(monthly_et_df$OBS_ET, monthly_et_df$MODIS_ET, 
     pch = 16, cex = 1, col = "#313695", main = "", 
     xlab = "Obs. ET (avg monthly mm/8day)", 
     ylab = "MODIS ET (avg monthly mm/8day)", 
     xlim = c(0,48), 
     ylim = c(0,48), 
     font = 1,
     font.lab = 2,
     cex.lab = 1.65, 
     cex.axis=1.5)
abline(lm(monthly_et_df$MODIS_ET_avg ~ monthly_et_df$OBS_ET_avg), col = "#313695", lwd = 3)
abline(0, 1, col = "#d73027", lty=2, lwd=3)














# bring in the Landsat data and do that as well 
cwt <- readOGR(paste0(home, "/Data/coweeta/cwt_flx_twr.shp"))
cwt_buff <- buffer(cwt, 300)

ec_et <- fread("G:/My Drive/Dissertation/ameriflux/et_mod_obs_df.csv")
ec_et <- ec_et[V1 >= as.Date("2011-01-01") & V1 < as.Date("2016-01-01")]

landsat_tifs <- list.files(paste0(home, "/Data_Testing/Landsat_ET/tifs_resampled"), full.names=T, pattern = ".tif$")
landsat_tifs_short <- list.files(paste0(home, "/Data_Testing/Landsat_ET/tifs_resampled"), full.names=F, pattern = ".tif$")


et_df <- data.table(file = landsat_tifs, 
                    start_date = as.Date(substr(landsat_tifs_short, 4, 11), "%Y%m%d"),
                    Landsat_ET = as.numeric(rep(NA, length(landsat_tifs))), 
                    countNA = as.numeric(rep(NA, length(landsat_tifs))))

et_df <- et_df[start_date >= as.Date("2011-01-01") & start_date <= as.Date("2015-01-01"), ]

for(i in 1:nrow(et_df)){
  ext_df <- raster::extract(raster(et_df$file[i]), cwt_buff, df=T)
  et_df$Landsat_ET[i] <- mean(ext_df[,2], na.rm=T)
  et_df$countNA[i] <- sum(is.na(ext_df[,2]))/nrow(ext_df)
  print(i)
}

et_df_sub <- et_df[!et_df$countNA == 1,]
et_df_sub_50 <- et_df_sub[et_df_sub$countNA < 0.2,]

et_df_sub_50 <- merge(et_df_sub_50, ec_et[,c("V1", "ET_obs_coweeta")], by.x="start_date", by.y = "V1", all.x=T)

metrics_fun(et_df_sub_50$ET_obs_coweeta, et_df_sub_50$Landsat_ET, "Overpass")



plot(et_df_sub_50$ET_obs_coweeta, et_df_sub_50$Landsat_ET, 
     pch = 16, cex = 1, col = "#313695", main = "", 
     xlab = "Obs. ET (mm/day)", 
     ylab = "Landsat ET (mm/day)", 
     xlim=c(0,7), 
     ylim=c(0,7), 
     cex.lab = 1.65, 
     font = 1,
     font.lab = 2,
     cex.axis=1.5)
abline(lm(et_df_sub_50$Landsat_ET ~ et_df_sub_50$ET_obs_coweeta), col = "#313695", lwd = 3)
abline(0, 1, col = "#d73027", lty=2, lwd=3)





# look at monthly averages 
et_df_sub_50$MY <- paste0(year(et_df_sub_50$start_date), "_", month(et_df_sub_50$start_date))
et_df_sub_50_mo <- et_df_sub_50[, .(Landsat_ET_mean = mean(Landsat_ET, na.rm=T), ET_obs_coweeta_mean = mean(ET_obs_coweeta, na.rm=T)), .(MY)]
metrics_fun(et_df_sub_50_mo$ET_obs_coweeta_mean, et_df_sub_50_mo$Landsat_ET_mean, "Monthly Average")


plot(et_df_sub_50_mo$ET_obs_coweeta_mean, et_df_sub_50_mo$Landsat_ET_mean, 
     pch = 16, cex = 1, col = "#313695", main = "", 
     xlab = "Obs. ET (avg monthly mm/day)", 
     ylab = "Landsat ET (avg monthly mm/day)", 
     xlim=c(0,7), 
     ylim=c(0,7), 
     cex.lab = 1.65, 
     cex.axis=1.5, 
     font = 1,
     font.lab = 2)
abline(lm(et_df_sub_50$Landsat_ET ~ et_df_sub_50$ET_obs_coweeta), col = "#313695", lwd = 3)
abline(0, 1, col = "#d73027", lty=2, lwd=3)

