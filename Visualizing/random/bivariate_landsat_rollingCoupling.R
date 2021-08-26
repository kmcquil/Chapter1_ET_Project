########################################################################################
## Landsat calculate the overall coupling, rolling 5 year coupling, and sen's slope of rolling coupling
## between growing season monthly ET anomalies and SPI 
########################################################################################

library(rgdal)
library(raster)
library(lubridate)
library(data.table)

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




########################################################################################
# Identify droughts and retrieve ET anomalies at drought peak 
########################################################################################

# make sure to jsut delete spi files that are outside of the modis scope 
# modis data is from 2000-01 : 2019-12

# bring in SPI and ET files. 
spi <- list.files(paste0(home, "/Data/SPI/retile"), full.names=T, pattern = ".tif$") 
spi_short <- list.files(paste0(home, "/Data/SPI/retile"), full.names=F, pattern = ".tif$") 
et <- list.files(paste0(home, "/Data/Landsat_ET/monthlyAnom"), full.names=T, pattern = ".tif$") 
et_short <- list.files(paste0(home, "/Data/Landsat_ET/monthlyAnom"), full.names=F, pattern = ".tif$") 


# convert each to a DT and get the month year 
spi <- data.table(spi_file = spi, date = as.Date( paste0(substr(spi_short,14,17), substr(spi_short,18,19), "01"), "%Y%m%d"), 
                  tile = substr(spi_short, 1, 12))
et <- data.table(et_file = et, date = as.Date( paste0(substr(et_short, 14, 17), substr(et_short, 18, 19), "01"), "%Y%m%d"), 
                 tile = substr(et_short, 1, 12))
# join the tables so that we only keep dates that both sources have 
files_dt <- merge(spi, et, by = c('date', 'tile'))

# get a list of the years 
years <- unique(year(files_dt$date))
tiles <- unique(files_dt$tile)


########################################################################################
# Fit the linear model between ETanom ~ SPI for each grid cell
########################################################################################

# subset the spi and et files to just growing season 
files_dt_gs <- files_dt[month(date) >3 & month(date) < 10]

# grab the forest mask file path for the modis data 
fmask <- list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), full.names = T, pattern = ".tif")

# calculate the linear models and store B0 and B1 coefficients in DT for future use 
library(Rmpi)
library(parallel)
library(snow)
workers <- mpi.universe.size() -1
cl <- makeMPIcluster(workers, type='MPI')

for(i in 1:length(tiles)){
  print("start")
  ptm <- proc.time()
  fmask_sub <- fmask[grep(tiles[i], fmask)]
  files_dt_gs_tl <- files_dt_gs[tile == tiles[i]]
  if(cellStats(raster(files_dt_gs_tl$spi_file[1]), min) == Inf ){next}
  
  rollingCoupling(files_dt_gs_tl$spi_file, files_dt_gs_tl$et_file, files_dt_gs_tl$date, fmask_sub, cl, "Landsat", tiles[i])
  
  print(proc.time() - ptm)
  print(i)
}

print("The coupling calcs are done!")


stopCluster(cl)
mpi.quit()






##############################################################################################################
## Get the results of the landsat trends in rolling coupling into one raster/column instead of multiple tiles 
##############################################################################################################
library(rgdal)
library(gdalUtils)
library(raster)
library(lubridate)
library(data.table)
library(doParallel)
library(foreach)
UseCores <- 10
cl <- makeCluster(UseCores)
registerDoParallel(cl)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))

# bring in the templates for individual tiles and the full roi
template <- paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif")
tile_template <- data.table(tile_file = list.files(paste0(home, "/Data/landcover/LANDSAT_FOREST/retile"), pattern = ".tif$", full.names=T))
tile_template$tile <- substr(tile_template$tile_file, nchar(tile_template$tile_file)-15, nchar(tile_template$tile_file)-4)

# list all of the .csv files correspondign to tile for sens slope and overall coupling 
files <- list.files(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling"), full.names=T, pattern = ".csv$")
rolling_files <- data.table(file = files[grep("BIV_rollingCoupling_x", files)])
rolling_files$tile <- substr(rolling_files$file, nchar(rolling_files$file) - 15, nchar(rolling_files$file)-4)
sens_files <- data.table(file = files[grep("BIV_coupling_sensslope_x", files)])
sens_files$tile <- substr(sens_files$file, nchar(sens_files$file)-15, nchar(sens_files$file)-4)

# join the tile templates with each of the dt with file of interest
rolling_files <- merge(rolling_files, tile_template, by = "tile", all.x=T)
sens_files <- merge(sens_files, tile_template, by = "tile", all.x=T)

# start with the sens_files since those will go quicker 
foreach(i = 1:nrow(sens_files))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(sens_files$file[i])
  dtIN <- dtIN[!is.na(pvalue)==T,]
  if(nrow(dtIN) == 0)return(NULL)
  
  sel_cols <- c("cellnum","senSlope")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(sens_files$tile_file[i])
  fout <- gsub(".csv", ".tif", sens_files$file[i])
  DTtoRast(dtIN, ExampleRaster, fout, open=F)
  
}

# save all of the sig sens slopes into one final .tif 
tif_files <- list.files(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling"), full.names=T, pattern = ".tif$")
rast_list <- gsub(".csv", ".tif", sens_files$file)
rast_list <- rast_list[rast_list %in% tif_files]
out <- paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/BIV_sig_sens_coupling_final.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)



# now grab the overall coupling and put into a raster  
foreach(i = 1:nrow(rolling_files))%dopar%{
  library(data.table)
  library(raster)
  home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
  source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
  
  dtIN <- fread(rolling_files$file[i])
  sel_cols <- c("cellnum","all")
  dtIN <- dtIN[,..sel_cols]
  ExampleRaster <- raster(rolling_files$tile_file[i])
  fout <- gsub(".csv", ".tif", rolling_files$file[i])
  DTtoRast(dtIN, ExampleRaster, fout, open=F)
  
}

rast_list <- gsub(".csv", ".tif", rolling_files$file)
out <- paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/BIV_sig_overall_coupling_final.tif")
align_rasters(rast_list, template, dstfile = out, overwrite = T)


