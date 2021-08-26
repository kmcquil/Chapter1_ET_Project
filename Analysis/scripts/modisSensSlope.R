library(data.table)
library(trend)
library(tidyr)
library(ggplot2)

home <- "G:/My Drive/Chapter1_ET_Project"
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

# bring in residual, ETanom, and SPI data tables
res_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanResiduals.csv"))
et_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanET.csv"))
spi_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanSPI.csv"))


res_dt_long <- gather(res_dt[,1:(ncol(res_dt)-1)], Year, ET_Residual, `2000`:`2019`)
et_dt_long <- gather(et_dt[,1:(ncol(et_dt)-1)], Year, ET_anom, `2000`:`2019`)
spi_dt_long <- gather(spi_dt[,1:(ncol(spi_dt)-1)], Year, SPI, `2000`:`2019`)

dt_long <- merge(res_dt_long, et_dt_long, by=c("cellnum", "Year"), all.x=T, all.y=T)
dt_long <- merge(dt_long, spi_dt_long, by = c("cellnum", "Year"), all.x=T, all.y=T)

dt_long <- setDT(dt_long)
dt_long$Year <- as.numeric(dt_long$Year)

#########################################################################################################
## Calculate sens slope and man kendall test of significance for ET residuals from 2000 - 2019 for all 
## droughts annual, 3yr, and 5yr

all_drought_sen <- sens_fun(dt_long)

dt_long_mod <- dt_long[SPI <= -1.3 & SPI > -1.6, ]
mod_drought_sen <- sens_fun(dt_long_mod)

dt_long_severe <- dt_long[SPI <= -1.6 & SPI > -2, ]
severe_drought_sen <- sens_fun(dt_long_severe)


dt_long_extreme <- dt_long[SPI <= -2,]
extreme_drought_sen <- sens_fun(dt_long_extreme)



##########################################################################################################
## Calculate sens slope and man kendall test of significane for every cell that has had at least 5 drought peaks 

countDroughtInCell <- dt_long[!is.na(SPI), .N, .(cellnum)]
cellsOfInterest <- countDroughtInCell[N >=5]$cellnum

dt_long_sub <- dt_long[cellnum %in% cellsOfInterest, .(cellnum, Year, ET_Residual)]

# put it back the wide way 
dt_wide_sub <- spread(dt_long_sub, Year, ET_Residual)
idx <- apply(dt_wide_sub[,2:ncol(dt_wide_sub)], 1, FUN = function(x) { sum(is.na(x))})
dt_wide_sub <- dt_wide_sub[!idx== (ncol(dt_wide_sub) - 1),]

lm_fun <- function(x){
  cellnum <- x[1]
  x <- as.numeric(x[2:length(x)])
  x <- x[!is.na(x)]
  s <- sens.slope(x)
  ss <- c(cellnum, s$estimates[1], s$p.value[1])

}

cellWiseSense <- apply(dt_wide_sub, 1, lm_fun)
cellWiseSense <- t(cellWiseSense)

sig_trend <- as.data.table(cellWiseSense[cellWiseSense[,3] < 0.05,])
colnames(sig_trend) <- c("cellnum", "slope", "pval")
ExampleRaster <- raster(paste0(home, "/Data/MODIS_ET/clean/et_20000101.tif"))
sel_cols <- c("cellnum", "slope")
sig_trend <- sig_trend[, ..sel_cols]
dr_trend_rast <- DTtoRast(sig_trend, ExampleRaster, paste0(home, "/Analysis/outputs/MODIS/sig_drought_trend.tif"), open = T)




