library(rgdal)
library(raster)
library(lubridate)
library(trend)
library(data.table)
library(parallel)
library(Rmpi)
library(parallel)
library(snow)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))


coupling_fun <- function(x){
  library(data.table)
  thresh <- length(x)/2
  spis <- x[1:thresh]
  ets <- x[(thresh+1):length(x)]
  sub <- data.table(spi=spis, et=ets)
  
  # initialize a vector for rho and p-value
  rho <- c()
  pval <- c()
  for(i in 1:(nrow(sub)-29)){
    subsub <- sub[i:(i+29),]
    if(sum(!is.na(subsub$et)) < 15){
      rho[i] <- NA
      pval[i] <- NA
    }else{
      subsubc <- subsub[complete.cases(subsub),]
      cortest <- cor.test(subsubc$et, subsubc$spi, method = "spearman")
      rho[i] <- cortest$estimate
      pval[i] <- cortest$p.value
    }
  }
  
  full <- sub[complete.cases(sub),]
  if(nrow(full) < 15){
    rho <- c(rho, NA)
    pval <- c(pval, NA)
    
  }else{
    cortest <- cor.test(full$et, full$spi, method = "spearman")
    rho <- c(rho, cortest$estimate)
    pval <- c(pval, cortest$p.value)
  }
  final_rho <- rho
  
  return(final_rho)
}


# this function calcualtes sens slope for each row with at least 25% of the coupling values significiant 
sens_row <- function(x){
  library(trend)
  x <- x[1:length(x)-1] # because the last row is the overall coupling 
  len <- round(length(x)/4, 0)
  x <- x[complete.cases(x)]
  print(length(x))
  
  if(length(x) == 0){
    c(NA, NA)
  }else if(length(x) < len){
    c(NA, NA)
  }else{
    yy <- sens.slope(x)
    c(yy$estimates, yy$p.value)
  }
}


rollingCoupling <- function(drought_files, ET_files, dates, fmask, cl, sensor, tl){
  
  spi_dt <- parSapply(cl, drought_files, rastToDF)
  spi_dt <- as.data.frame(do.call(cbind, spi_dt))
  print("finished creating spi dataframe")
  
  et_dt <- parSapply(cl, ET_files, rastToDF)
  et_dt <- as.data.frame(do.call(cbind, et_dt))
  print("finished creating et dataframe")
  
  cellnum <- data.frame(cellnum = seq(1, nrow(spi_dt)))
  spi_dt <- cbind(cellnum, spi_dt)
  et_dt <- cbind(cellnum, et_dt)
  
  # convert forest mask to a dataframe
  f_dt <- as.data.frame(raster(fmask))
  colnames(f_dt) <- c("fmask")
  
  # drop rows from each dataframe where forest mask == NA
  spi_dt <- spi_dt[!(is.na(f_dt$fmask) | f_dt$fmask == 0),]
  et_dt <- et_dt[!(is.na(f_dt$fmask) | f_dt$fmask == 0),]
  print(dim(spi_dt))
  print(dim(et_dt))
  
  combine_dt <- cbind(spi_dt[,2:ncol(spi_dt)], et_dt[,2:ncol(et_dt)])
  r <- round(nrow(combine_dt)/4, 0)
  ends <- c(r, r*2,r*3)
  combine_dt1 <- combine_dt[1:ends[1],]
  combine_dt2 <- combine_dt[(ends[1] + 1):ends[2],]
  combine_dt3 <- combine_dt[(ends[2] + 1):ends[3],]
  combine_dt4 <- combine_dt[(ends[3] + 1):nrow(combine_dt),]
  
  if(nrow(combine_dt) == 0){
    return()
  }else{
    beta_dt1 <- parApply(cl, combine_dt1, FUN=coupling_fun, MARGIN = 1)
    beta_dt2 <- parApply(cl, combine_dt2, FUN=coupling_fun, MARGIN=1)
    beta_dt3 <- parApply(cl, combine_dt3, FUN=coupling_fun, MARGIN = 1)
    beta_dt4 <- parApply(cl, combine_dt4, FUN=coupling_fun, MARGIN=1)
    beta_dt1 <- t(beta_dt1)
    beta_dt2 <- t(beta_dt2)
    beta_dt3 <- t(beta_dt3)
    beta_dt4 <- t(beta_dt4)
    
    beta_dt <- rbind(beta_dt1, beta_dt2, beta_dt3, beta_dt4)
  }
  beta_dt <- cbind(spi_dt$cellnum, beta_dt)
  dates <- dates[1:(length(dates) - 29)]
  dates <- paste0("d",substr(dates, 1, 4), substr(dates, 6, 7))
  colnames(beta_dt) <- c("cellnum", dates, "all")
  print(paste0("The size of the DT being written is ", nrow(beta_dt), "    ",ncol(beta_dt)))
  
  if(is.na(tl) == T){
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/BIV_rollingCoupling.csv"))
  }else{
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/BIV_rollingCoupling_", tl,".csv"))
  }
  
  
  sens_dt1 <- t(parApply(cl, beta_dt1, FUN=sens_row, MARGIN = 1))
  sens_dt2 <- t(parApply(cl, beta_dt2, FUN=sens_row, MARGIN = 1))
  sens_dt3 <- t(parApply(cl, beta_dt3, FUN=sens_row, MARGIN = 1))
  sens_dt4 <- t(parApply(cl, beta_dt4, FUN=sens_row, MARGIN = 1))
  
  sens_dt <- rbind(sens_dt1, sens_dt2, sens_dt3, sens_dt4)
  sens_dt <- cbind(spi_dt$cellnum, sens_dt)
  colnames(sens_dt) <- c("cellnum", "senSlope", "pvalue")
  
  if(is.na(tl) == T){
    fwrite(sens_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/BIV_coupling_sensslope.csv"))
  }else{
    fwrite(sens_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/BIV_coupling_sensslope_", tl,".csv"))
  }
}




workers <- mpi.universe.size() -1
cl <- makeMPIcluster(workers, type='MPI')

# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data/SPI/SPI90_MODIS"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data/SPI/SPI90_MODIS"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data/MODIS_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 

# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,1,4), substr(spi_short,5,6), "01"), "%Y%m%d"))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 4, 7), substr(et_short, 8, 9), "01"), "%Y%m%d"))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = 'date')
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

fmask <- paste0(home, "/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif")

rollingCoupling(files_dt_gs$spi_file, files_dt_gs$et_file, files_dt_gs$date, fmask, cl, "MODIS", tl = NA)

stopCluster(cl)
mpi.quit()





home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

et_tifs <- list.files(paste0(home, "/Analysis/outputs/MODIS/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs <- et_tifs[grep("ETanom", et_tifs)]
anom_stk <- do.call("stack", lapply(anom_tifs, raster))



# bring in the pixel wise sens slope of the 5 year rolling coupling and mask out sens slope for pixels that are not significant 
coupling_slope <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_coupling_sensslope.csv"))
coupling_slope <- coupling_slope[coupling_slope$cellnum %in% forest_mask1$cellnum,]

# convert the significant sens slope to a raster
sig_slope <- DTtoRast(coupling_slope[,c(1,2)], anom_stk[[1]], 
                      paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_slope.tif"), open = T)

# bring in the 5 year rolling coupling and the overall coupling 
rollingCoupling <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_rollingCoupling.csv"))
rollingCoupling <- rollingCoupling[rollingCoupling$cellnum %in% forest_mask1$cellnum,]

# convert the overall coupling to a raster 
overall <- DTtoRast(rollingCoupling[,c("cellnum", "all")], anom_stk[[1]], 
         paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_overallCoupling.tif"), open = T)

