#########################################################################################################################
## Functions for the analysis 
#########################################################################################################################

#########################################################################################################################
## Function: Identify the drought peak and associated anomaly for a given year if drought peak < -1.3 ( at least moderate drought)
## Inputs: 
#### drought_files: vector of monthly spi raster files 
#### ET_files: vector of monthly ET anomaly files
#### n: the number of cores you want to use
#### home: the home path 
#### sensor: "MODIS", or "Landsat"
## Outputs: 
#### write raster of the month index of drought peak 
#### write raster of the value of SPI at drought peak 
#### write raster of the value of ET anomaly at drought peak 
##########################################################################################################################

droughtID <- function(drought_files, ET_files, n, home, y, sensor){
  beginCluster(n = n)
  # bring in spi files list as a raster stack 
  spi_stack <- do.call("stack", lapply(drought_files, raster))
  # bring in et files list as a raster stack -- must match the date order as spi
  et_stack <- do.call("stack", lapply(ET_files, raster))
  
  # function to calculate the index of the drought peak 
  peakFun <- function(x){ if(sum(!is.na(x)) == 0){NA}else{z<-which(x==min(x)); z[length(z)]}}
  # function to calcualte the magnitude SPI at drought peak 
  peakMagFun <- function(x){ if(sum(!is.na(x)) == 0){NA}else{z<-which(x==min(x)); idx = z[length(z)]; x[idx]}}
  
  
  # get the index of drought peak 
  peakIdx <- clusterR(spi_stack, fun=calc, args=list(fun=peakFun))
  # get the spi value at drought peak 
  peakMag <- clusterR(spi_stack, fun=calc, args=list(fun=peakMagFun))
  # if the spi value at drought peak > -1.3, it is not severe enough and we will not include it 
  peakMag[peakMag > -1.3] <- NA
  # use the peakMag raster to mask out not severe enough droughts from the index raster
  peakIdx <- mask(peakIdx, peakMag)
  
  # to find ET anom at peak start by stacking index raster with ET anom raster 
  combo <- stack(peakIdx, et_stack)
  peakETFun <- function(x){ if(is.na(x[1] == T)){NA}else{idx=x[1]+1; x[idx] }}
  peakET <- clusterR(combo, fun = calc, args=list(fun=peakETFun))
  
  # write out the results
  writeRaster(peakMag, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakSPI_", y, ".tif"), format = "GTiff", overwrite = T)
  writeRaster(peakIdx, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakIndex_", y, ".tif"), format = "GTiff", overwrite = T)
  writeRaster(peakET, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakETanom_", y, ".tif"), format = "GTiff", overwrite = T)
  
  endCluster()
  }



#########################################################################################################################
## Function: Identify the drought peak and associated anomaly for a given year if drought peak < -1.3 (moderate drought)
## This is essentiallly the same as droughtID function but for use on the landsat data broken down into tiles since
## they were too large to process all at once 
## Inputs: 
#### drought_files: vector of monthly spi raster files 
#### ET_files: vector of monthly ET anomaly files
#### n: the number of cores you want to use
#### home: the home path 
#### y: year
#### tl: tile
#### sensor: "MODIS", or "Landsat"
## Outputs: 
#### write raster of the month index of drought peak for a given year and tile 
#### write raster of the value of SPI at drought peak for a given year and tile
#### write raster of the value of ET anomaly at drought peak for a given year and tile 
##########################################################################################################################

droughtID_tl <- function(drought_files, ET_files, n, home, y, tl, sensor){
  beginCluster(n = n)
  # bring in spi files list as a raster stack 
  spi_stack <- do.call("stack", lapply(drought_files, raster))
  # bring in et files list as a raster stack -- must match the date order as spi
  et_stack <- do.call("stack", lapply(ET_files, raster))
  
  # function to calculate the index of the drought peak 
  peakFun <- function(x){ if(sum(!is.na(x)) == 0){NA}else{z<-which(x==min(x)); z[length(z)]}}
  # function to calcualte the magnitude SPI at drought peak 
  peakMagFun <- function(x){ if(sum(!is.na(x)) == 0){NA}else{z<-which(x==min(x)); idx = z[length(z)]; x[idx]}}
  
  
  # get the index of drought peak 
  peakIdx <- clusterR(spi_stack, fun=calc, args=list(fun=peakFun))
  # get the spi value at drought peak 
  peakMag <- clusterR(spi_stack, fun=calc, args=list(fun=peakMagFun))
  # if the spi value at drought peak > -1.3, it is not severe enough and we will not include it 
  peakMag[peakMag > -1.3] <- NA
  # use the peakMag raster to mask out not severe enough droughts from the index raster
  peakIdx <- mask(peakIdx, peakMag)
  
  # to find ET anom at peak start by stacking index raster with ET anom raster 
  combo <- stack(peakIdx, et_stack)
  peakETFun <- function(x){ if(is.na(x[1] == T)){NA}else{idx=x[1]+1; x[idx] }}
  peakET <- clusterR(combo, fun = calc, args=list(fun=peakETFun))
  
  # write out the results
  writeRaster(peakMag, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakSPI_", y, "_", tl,".tif"), format = "GTiff", overwrite = T)
  writeRaster(peakIdx, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakIndex_", y, "_", tl,".tif"), format = "GTiff", overwrite = T)
  writeRaster(peakET, paste0(home, "/Analysis/outputs/", sensor, "/droughtPeak/droughtPeakETanom_", y, "_", tl,".tif"), format = "GTiff", overwrite = T)
  
  endCluster()
}



#########################################################################################################################
## Function: Fit a linear model of ET anom ~ SPI for growing season months for each grid cell 
## Inputs: 
#### drought_files: spi monthly raster files 
#### ET_files: et anomaly monthly raster files 
#### forest mask: file for the forest mask 
#### n = number of cores to use 
#### sensor = "MODIS" or "Landsat"
## Outputs: 
#### write out the beta coefficients with cell number for ID
##########################################################################################################################

modFit <- function(drought_files, ET_files, fmask, n, sensor){
  library(foreach)
  library(doParallel)
  UseCores <- n
  cl <- makeCluster(UseCores)
  registerDoParallel(cl)
  
  # bring all SPI raster files into a dataframe where each row is a cell and each column is a growing season month 
  spi_dt <- foreach(i = 1:length(drought_files), .combine = cbind)%dopar%{
    library(raster)
    as.data.frame(raster(drought_files[i]))
  }
  print("finished creating spi dataframe")
  
  # do the same for the ET raster files 
  et_dt <- foreach(i = 1:length(ET_files), .combine = cbind)%dopar%{
    library(raster)
    as.data.frame(raster(ET_files[i]))
  }
  print("finished creating et dataframe")
  
  # add a column called cell number 
  cellnum <- data.frame(cellnum = seq(1, nrow(spi_dt)))
  spi_dt <- cbind(cellnum, spi_dt)
  et_dt <- cbind(cellnum, et_dt)
  
  # convert forest mask to a dataframe
  f_dt <- as.data.frame(raster(fmask))
  colnames(f_dt) <- c("fmask")
  
  # keep rows where the forest mask == 1, indicating forest 
  spi_dt <- spi_dt[f_dt$fmask == 1,]
  et_dt <- et_dt[f_dt$fmask == 1,]
  
  # fit a linearmodel for each corresonding row (pixel) and put the model B0 and B1 coefficients into a dt
  beta_dt <- foreach(i = 1:nrow(spi_dt), .combine=rbind)%dopar%{
    library(data.table)
    # get the et and spi rows which correspond to the same cell and put into a dt and keep only complete cases 
    # we only want complete cases because we will use them to fit a linear model
    sub <- data.table(spi=as.numeric(spi_dt[i,]), et=as.numeric(et_dt[i,])) 
    sub <- sub[complete.cases(sub),]
    
    # fit linear model between ET ~ SPI and put B0 and B1 into dt 
    # some pixels have no data across the time series, in that case just insert NA values 
    if(dim(sub)[1] == 0){
      betas <- c(NA, NA)
    }else{
      linmod <- lm(sub$et ~ sub$spi)
      betas <- c(linmod$coefficients[1], linmod$coefficients[2])
    }
    
  }
  
  # add the cellnumber to the beta coefficient DT so that we can match the coefficients to the correct cell later 
  beta_dt <- cbind(spi_dt$cellnum, beta_dt)
  colnames(beta_dt) <- c("cellnum", "B0", "B1")
  rownames(beta_dt) <- c()
  # save the beta coefficient DT 
  fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/model_beta_coef.csv"), nThread = n)
  
  stopCluster(cl)
  
  rm(spi_dt)
  rm(et_dt)
  return()
}

##########################################################################################################################
## Function: Take a .tif file and convert it to a DF
## Input: 
#### tif: the file path of the raster to be converted to a data frame 
## Output: 
#### the new dataframe 
##########################################################################################################################
rastToDF <- function(tif){
  library(raster)
  as.data.frame(raster(tif))
}

#########################################################################################################################
## Function: Fit a linear model of ET anom ~ SPI for growing season months for each grid cell to the Landsat data that 
## has been broken into tiles. This is essentially the same as function above but broken into tiles because the landsat 
## rasters were too big to process all at once 
## Inputs: 
#### drought_files: spi monthly raster files 
#### ET_files: et anomaly monthly raster files 
#### fmask: file for the forest mask 
#### cl = the parallel cluster -- use Rmpi 
#### sensor = "MODIS" or "Landsat"
#### tl = tile ID used to pull in correct files 
## Outputs: 
#### write out the beta coefficeints with cell number and ID by tile 
##########################################################################################################################

##########################################################################################################################
## Function: Calculate the beta coefficients within the modFit_tl function 
## Input: 
#### x: one long vector that has the SPI values and then ET values in chronological order for a pixel. Put ET and SPI into one vector to work with Rmpi 
## Output: 
#### vector with the B0 and B1 coefficeints for that particular pixel 
##########################################################################################################################
betas_fun <- function(x){
  library(data.table)
  thresh <- length(x)/2 # identify the threshold of SPI vs ET values 
  spis <- x[1:thresh] # grab the spi 
  ets <- x[(thresh+1):length(x)] # grab the ET anomalies 
  sub <- data.table(spi=spis, et=ets) # put into a DT and then drop incomplete rows. This will be used to fit linear model
  sub <- sub[complete.cases(sub),]
  
  # just return NAs if there are no complete rows. Otherwise fit the model and return coefficients 
  if(dim(sub)[1] == 0){
    c(NA, NA)
  }else{
    linmod <- lm(sub$et ~ sub$spi)
    c(linmod$coefficients[1], linmod$coefficients[2])
  }
}


modFit_tl <- function(drought_files, ET_files, fmask, cl, sensor, tl){
  
  # create a DF of SPI values. Each row is a cell and each column is a month 
  spi_dt <- parSapply(cl, drought_files, rastToDF)
  spi_dt <- as.data.frame(do.call(cbind, spi_dt))
  print("finished creating spi dataframe")
  
  # same for ET 
  et_dt <- parSapply(cl, ET_files, rastToDF)
  et_dt <- as.data.frame(do.call(cbind, et_dt))
  print("finished creating et dataframe")
  
  # add a column of cell number to each DF 
  cellnum <- data.frame(cellnum = seq(1, nrow(spi_dt)))
  spi_dt <- cbind(cellnum, spi_dt)
  et_dt <- cbind(cellnum, et_dt)
  
  # convert forest mask to a dataframe
  f_dt <- as.data.frame(raster(fmask))
  colnames(f_dt) <- c("fmask")
  
  # mask out non-forest pixels 
  spi_dt <- spi_dt[!(is.na(f_dt$fmask) | f_dt$fmask == 0),]
  et_dt <- et_dt[!(is.na(f_dt$fmask) | f_dt$fmask == 0),]
  
  print(dim(spi_dt))
  print(dim(et_dt))
  
  # combine the data tables but don't include cell number 
  combine_dt <- cbind(spi_dt[,2:ncol(spi_dt)], et_dt[,2:ncol(et_dt)])
  # since rmpi gets mad when I try to send the full combine_dt, break it up into four separate DFs
  r <- round(nrow(combine_dt)/4, 0)
  ends <- c(r, r*2,r*3)
  combine_dt1 <- combine_dt[1:ends[1],]
  combine_dt2 <- combine_dt[(ends[1] + 1):ends[2],]
  combine_dt3 <- combine_dt[(ends[2] + 1):ends[3],]
  combine_dt4 <- combine_dt[(ends[3] + 1):nrow(combine_dt),]
  
  # fit the linear model and return the B0 and B1 coefficeints using betas_fun
  # if the initial Dt was empty just return a blank 
  if(nrow(combine_dt) == 0){
    return()
  }else{
    beta_dt1 <- parApply(cl, combine_dt1, FUN=betas_fun, MARGIN = 1) 
    beta_dt2 <- parApply(cl, combine_dt2, FUN=betas_fun, MARGIN=1)
    beta_dt3 <- parApply(cl, combine_dt3, FUN=betas_fun, MARGIN = 1)
    beta_dt4 <- parApply(cl, combine_dt4, FUN=betas_fun, MARGIN=1)
    # parApply returns to rows. Use transpose to convert to columns (B0 and B1)
    beta_dt1 <- t(beta_dt1)
    beta_dt2 <- t(beta_dt2)
    beta_dt3 <- t(beta_dt3)
    beta_dt4 <- t(beta_dt4)
    
    # bring them all together. 
    beta_dt <- rbind(beta_dt1, beta_dt2, beta_dt3, beta_dt4)
  }
  print("finished calculating betas")
  
  # add cellnumber to make it possible to match beta coefficients to appropriate cell in next step 
  beta_dt <- cbind(spi_dt$cellnum, beta_dt)
  colnames(beta_dt) <- c("cellnum", "B0", "B1")
  rownames(beta_dt) <- c()
  
  # save the beta coefficeints DF and name using tile ID 
  fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/model_beta_coef", "_", tl,".csv"))
  
}


#########################################################################################################################
## Function: Calculate the residuals at drought peak 
## Inputs: 
#### spi_file: spi at drought peak 
#### et_file: et anom at drought peak 
#### beta_file: file for the beta coefficient data frame
#### y = year
#### sensor = "MODIS" or "Landsat"
## Outputs: 
#### a csv of ET anomalies at drought peaak for each year 
#### a csv of SPI at drought peaak for each year 
#### a csv of ET residuals at drought peak for each year 
##########################################################################################################################

calcResiduals <- function(spi_file, et_file, beta_file, y, sensor){
  # take spi file and et file and convert them to rasters, then dataframes, add cellnumber and merge together
  spi_dt <- as.data.frame(raster(spi_file))
  et_dt <- as.data.frame(raster(et_file))
  cellnum <- data.table(cellnum = seq(1, nrow(spi_dt)))
  combo <- cbind(cellnum, spi_dt, et_dt)
  
  # bring in the beta file as a csv 
  betas <- fread(beta_file)
  
  # join them by cell number 
  combo <- merge(combo, betas, by ="cellnum")
  colnames(combo) <- c("cellnum", "spi", "et", "B0", "B1")
  
  # calculate the predicted ET anom using the beta coefficients 
  combo$pred <- as.numeric(rep(NA, nrow(combo)))
  combo[, pred:=(spi*B1)+B0]

  # now calculate the residual (actual - predicted)
  # the residual is the ET response with drought severity impact removed 
  combo$residual <- combo$et - combo$pred
  
  # create a data.table with the cell number and residual and save as a .csv 
  fin <- combo[,c("cellnum", "residual")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/residuals/", y, ".csv"))
  
  # create a data.table with the cell number and ETanomaly and save as a .csv 
  fin <- combo[,c("cellnum", "et")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/ETanom/", y, ".csv"))
  
  # create a data.table with the cell number and SPI value and save as a .csv 
  fin <- combo[,c("cellnum", "spi")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/spi/", y, ".csv"))
}



#########################################################################################################################
## Function: Calculate the residuals at drought peak for Landsat data that was broken into tiles 
## Inputs: 
#### spi_file: spi at drought peak 
#### et_file: et anom at drought peak 
#### beta_file: file for the beta coefficient data frame
#### y = year
#### tl = tile 
#### sensor = "MODIS" or "Landsat"
## Outputs: 
#### a csv of ET anomalies at drought peaak for each year 
#### a csv of SPI at drought peaak for each year 
#### a csv of ET residuals at drought peak for each year 
##########################################################################################################################

calcResiduals_tl <- function(spi_file, et_file, beta_file, y,tl, sensor){
  # take spi file and et file and convert them to rasters, then dataframes, add cellnumber and merge together
  spi_dt <- as.data.frame(raster(spi_file))
  et_dt <- as.data.frame(raster(et_file))
  cellnum <- data.table(cellnum = seq(1, nrow(spi_dt)))
  combo <- cbind(cellnum, spi_dt, et_dt)
  
  # bring in the beta file as a csv 
  betas <- fread(beta_file)
  
  # join them by cell number 
  combo <- merge(combo, betas, by ="cellnum")
  colnames(combo) <- c("cellnum", "spi", "et", "B0", "B1")
  
  # calculate the predicted ET anom using the beta coefficients 
  combo$pred <- as.numeric(rep(NA, nrow(combo)))
  combo[, pred:=(spi*B1)+B0]
  
  # now calculate the residual (actual - predicted)
  # the residual is the ET response with drought severity impact removed 
  combo$residual <- combo$et - combo$pred
  
  # create a data.table with the cell number and residual and save as a .csv 
  fin <- combo[,c("cellnum", "residual")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/residuals/", y, "_", tl, ".csv"))
  
  # create a data.table with the cell number and ETanomaly and save as a .csv 
  fin <- combo[,c("cellnum", "et")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/ETanom/", y, "_", tl, ".csv"))
  
  # create a data.table with the cell number and SPI value and save as a .csv 
  fin <- combo[,c("cellnum", "spi")]
  colnames(fin) <- c("cellnum", y)
  fwrite(fin, paste0(home, "/Analysis/outputs/", sensor, "/spi/", y, "_", tl, ".csv"))
}


#########################################################################################################################
## Function: Take all of the residuals/ET anomalies/SPI and find the average by cell to get the average drought impact 
## Inputs: 
#### resid_files: files of residuals/ET anomalies/SPI for every year with drought
#### n = number of threads (just to open and close all of those csv files)
#### sensor = "MODIS" or "Landsat"
## Outputs: 
#### a csv of the yearly and total average residuals in one dataframe 
##########################################################################################################################

# function to grab the actual residuals by cell for at least a moderate drought 
avgResid <- function(resid_files, n, sensor){
  resid_dt <- fread(resid_files[1]) # bring in the first residual file 
  for(i in 2:length(resid_files)){
    new_dt <- fread(resid_files[i], nThread = n)
    resid_dt <- merge(resid_dt, new_dt, by = "cellnum", all.x=T, all.y=T) # keep merging the residual files in 
  }
  
  resid_dt$meanResid <- as.numeric(rep(NA, nrow(resid_dt))) # initialize a column for the mean 
  sel_cols <- resid_dt[,names(resid_dt)[seq(2, ncol(resid_dt)-1)]] #get the columns to use in mean calc
  resid_dt[,meanResid := rowMeans(.SD, na.rm=T), .SDcols = sel_cols ] # calculate row means 
  fwrite(resid_dt, paste0(home, "/Analysis/outputs/", sensor, "/meanResiduals.csv"), nThread = n) # write out the csv
  
}

# same function but for the average ET anomaly
avgET <- function(resid_files, n, sensor){
  resid_dt <- fread(resid_files[1])
  for(i in 2:length(resid_files)){
    new_dt <- fread(resid_files[i], nThread = n)
    resid_dt <- merge(resid_dt, new_dt, by = "cellnum", all.x=T, all.y=T)
  }
  
  resid_dt$meanET <- as.numeric(rep(NA, nrow(resid_dt)))
  sel_cols <- resid_dt[,names(resid_dt)[seq(2, ncol(resid_dt)-1)]]
  resid_dt[,meanET := rowMeans(.SD, na.rm=T), .SDcols = sel_cols ]
  fwrite(resid_dt, paste0(home, "/Analysis/outputs/", sensor, "/meanET.csv"), nThread = n)
  
}

# same but for average SPI 
avgSPI <- function(resid_files, n, sensor){
  resid_dt <- fread(resid_files[1])
  for(i in 2:length(resid_files)){
    new_dt <- fread(resid_files[i], nThread = n)
    resid_dt <- merge(resid_dt, new_dt, by = "cellnum", all.x=T, all.y=T)
  }
  
  resid_dt$meanSPI <- as.numeric(rep(NA, nrow(resid_dt)))
  sel_cols <- resid_dt[,names(resid_dt)[seq(2, ncol(resid_dt)-1)]]
  resid_dt[,meanSPI := rowMeans(.SD, na.rm=T), .SDcols = sel_cols ]
  fwrite(resid_dt, paste0(home, "/Analysis/outputs/", sensor, "/meanSPI.csv"), nThread = n)
  
}


##############################################################################################################
## Function to convert a data.table with a column for cellnum and the variable of interest back into a raster
## Inputs: 
#### dtIN: data.table with column for cellnum and column of interest to convert to raster 
#### ExampleRaster: a raster with desired crs, extent, and resolution
#### fout: filename to write out the raster 
#### open: A true or false, do you want to also open the raster rn? 
## Outputs:
##### A raster 
###############################################################################################################

DTtoRast <- function(dtIN, ExampleRaster, fout, open){
  proj <- crs(ExampleRaster)
  cellnum <- data.table(cellnum = seq(1, ncell(ExampleRaster)))
  
  colnames(dtIN) <- c("cellnum", "varInt")
  dtIN_merged <- merge(cellnum, dtIN, by = "cellnum", all.x=T)
  
  output_matrix <- matrix(dtIN_merged$varInt, nrow=nrow(ExampleRaster),ncol=ncol(ExampleRaster),byrow=T)
  new_output_raster<-raster(output_matrix,xmn=xmin(ExampleRaster),ymn=ymin(ExampleRaster),xmx=xmax(ExampleRaster),
                            ymx=ymax(ExampleRaster), crs=proj)
  
  writeRaster(new_output_raster, fout, overwrite = T, format = "GTiff")
  if(open == TRUE){return(new_output_raster)}else{return()}
}


##############################################################################################################
## Function to find the average of ETanom or residuals based on SPI severity
## Inputs: 
#### varIN: data.table ETanom or residuals
#### spiIN: data.table of spi values used to mask ET -- column names must match (cellnum - years)
#### lowerBound: SPI lower bound (not included)
#### upperBound: SPI upper bound (included)
## Outputs:
##### a data.table with a column for cellnum and the average by cell for the cells meeting the spi requriement
###############################################################################################################

avgByDroughtSeverity <- function(varIN, spiIN, lowerBound, upperBound){
  
  # use conditions on spi to filter the variable of interest 
  varIN <- copy(varIN)
  spiIN <- copy(spiIN)
  for(col in names(spiIN[,2:ncol(spiIN)])) set(varIN, i = which(spiIN[[col]] > upperBound | spiIN[[col]] <= lowerBound), j=col, value=NA)
  
  # create a new mean column using the upper and lower bounds 
  varIN[[paste0("mean_", lowerBound, "_", upperBound)]] <- as.numeric(rep(NA, nrow(varIN)))
  sel_cols <- varIN[,names(varIN)[seq(2, ncol(varIN)-1)]]
  varIN[,paste0("mean_", lowerBound, "_", upperBound) := rowMeans(.SD, na.rm=T), .SDcols = sel_cols]
  
  keep <- c("cellnum", paste0("mean_", lowerBound, "_", upperBound))
  keepDT <- varIN[, ..keep]
  return(keepDT)
  
}



####################################################################################################
## Function to take the landsat tiles results at drought peak for each year and convert them to one full result
## Inputs: 
#### inFiles <- list of files of all tiles and all years for residuals, anoms, or spi
#### tile_template <- templates of each individual tile 
#### template <- template to regrid rasters onto 
#### output_dir <- file path to output directory 
#### n = number of cores to use 
## Outputs:
#### rasters of each year 
#### a dataframe of each year that matches the MODIS df 
#####################################################################################################

retiling <- function(inFiles, tile_template, template, output_dir, n){
  
  library(foreach)
  library(doParallel)
  UseCores <- n
  cl <- makeCluster(UseCores)
  registerDoParallel(cl)
  
  # convert list to a data table with a column for file, date, tile
  dt <- data.table(file = inFiles, 
                   tile = substr(inFiles, nchar(inFiles)-15, nchar(inFiles)-4),
                   year = substr(inFiles, nchar(inFiles)-20, nchar(inFiles)-17))
  # create another column that subs .csv to .tif 
  dt$tif_file <- gsub(".csv", ".tif", dt$file)
  
  # create another column with the file name of the landcover that corresponds to that tile 
  tile_template <- data.table(tile_file = tile_template, 
                              tile = substr(tile_template, nchar(tile_template) - 15, nchar(tile_template)-4))
  dt <- merge(dt, tile_template, by = "tile", all.x=T)
  print(dim(dt))
  # loop through unique date 
  ### Convert each file to a raster with the same name 
  ### Regrid all of the new rasters to the template 
  ### Convert the new full raster to a data.table with columns cellnum and the year 
  
  foreach(i = 1:nrow(dt))%dopar%{
    library(data.table)
    library(raster)
    home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
    #home <- "G:/My Drive/Chapter1_ET_Project"
    source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
    
    
    dtIN <- fread(dt$file[i])
    ExampleRaster <- raster(dt$tile_file[i])
    fout <- dt$tif_file[i]
    DTtoRast(dtIN, ExampleRaster, fout, open=F)
    #print(i)
  }
  
  years <- unique(dt$year)
  print(length(years))
  foreach(i = 1:length(years))%dopar%{
    library(data.table)
    library(raster)
    library(gdalUtils)
    home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
    #home <- "G:/My Drive/Chapter1_ET_Project"
    source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
    
    
    rast_list <- dt[year == years[i]]$tif_file
    out <- paste0(output_dir, years[i], ".tif")
    align_rasters(rast_list, template, dstfile = out, overwrite = T)
    
    outcsv <- gsub(".tif", ".csv", out)
    year_dt <- as.data.frame(raster(out))
    cellnum <- as.data.frame(seq(1, nrow(year_dt)))
    year_dt <- cbind(cellnum, year_dt)
    colnames(year_dt) <- c("cellnum", years[i])
    fwrite(year_dt, outcsv)
    #print(i)
  }
  
  stopCluster(cl)
}



####################################################################################################
## Function: Calculate the average ET anomaly and ET residual during months classified as nomral wetness 
## using SPI on a pixel basis. Normal wetness is spi <= 0.49 & spi >= -0.49
## Inputs:
#### drought_files: monthly raster SPI files 
#### ET_files: monthly raster ET files 
#### betas_file: the csv with linear model coefficients for each cell 
#### fmask: the forest mask 
#### cl: the cluster for parallel computing 
#### sensor: which sensor? (MODIS or Landsat)
#### tl: for Landsat use the tile ID; for MODIS put NA
## Outputs:
#### a csv with columns: cell number, average ET anomaly, and average ET residual 
#####################################################################################################

####################################################################################################
## Function: Calculate the average ET anomaly and ET residual during months classified as nomral wetness 
## This is used within the avgNA_tl function 
## Inputs:
#### x: a vector with SPI and ET values corresponding to one pixel. The B0 and B1 coeffficents are between the SPI and ET values in the vector. 
###### I put these all in one vector to make it easier to use the values with Rmpi functions 
## Outputs:
#### a vector with the average ET anomaly and average ET residual during normal wetness periods 
#####################################################################################################

avg_fun <- function(x){
  library(data.table)
  thresh <- (length(x)-2)/2
  spis <- x[1:thresh]
  ets <- x[(thresh+3):length(x)]
  B0 <- x[thresh+1]
  B1 <- x[thresh+2]
  sub <- data.table(spi=spis, et=ets)
  sub <- sub[complete.cases(sub),]
  sub <- sub[spi <= 0.49 & spi >= -0.49,]
  sub$residual <- (sub$et*B1) + B0
  
  if(dim(sub)[1] == 0){
    c(NA, NA)
  }else{
    etAnom <- mean(sub$et, na.rm=T)
    etRes <- mean(sub$residual)
    c(etAnom, etRes)
  }
  
}

avgNA_tl <- function(drought_files, ET_files, betas_file, fmask, cl, sensor, tl){
  
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
  spi_dt <- spi_dt[!is.na(f_dt$fmask),]
  et_dt <- et_dt[!is.na(f_dt$fmask),]
  print(dim(spi_dt))
  print(dim(et_dt))
  
  # join the betas file to the ET file 
  betas <- fread(betas_file)
  betas <- betas[complete.cases(betas),]
  et_dt <- merge(betas, et_dt, by='cellnum', all.y = T)
  
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
    beta_dt1 <- parApply(cl, combine_dt1, FUN=avg_fun, MARGIN = 1)
    beta_dt2 <- parApply(cl, combine_dt2, FUN=avg_fun, MARGIN=1)
    beta_dt3 <- parApply(cl, combine_dt3, FUN=avg_fun, MARGIN = 1)
    beta_dt4 <- parApply(cl, combine_dt4, FUN=avg_fun, MARGIN=1)
    beta_dt1 <- t(beta_dt1)
    beta_dt2 <- t(beta_dt2)
    beta_dt3 <- t(beta_dt3)
    beta_dt4 <- t(beta_dt4)
    
    beta_dt <- rbind(beta_dt1, beta_dt2, beta_dt3, beta_dt4)
  }
  print("finished calculating avg ET!")
  
  beta_dt <- cbind(spi_dt$cellnum, beta_dt)
  colnames(beta_dt) <- c("cellnum", "etAnomAvg", "etResAvg")
  rownames(beta_dt) <- c()
  
  if(is.na(tl) == T){
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/avgNonDroughtETResponse.csv"))
  }else{
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/avgNonDroughtETResponse_", tl,".csv"))
  }
  
}


####################################################################################################
## Function: Calculate the rolling and overall coupling of monthly SPI and monthly ET anomalies 
## Only calculate coupling if at least half of observations for the time period are available for the rolling coupling 
## Only calculate overall coupling at least 15 of observations were available 
## This is used within the rollingCoupling function 
## Inputs:
#### x: a vector with SPI and ET values corresponding to one pixel.  
###### I put these all in one vector to make it easier to use the values with Rmpi functions 
## Outputs:
#### a vector with the significant rolling coupling coefficients at each time step and the overall coupling last. 
###### Only correlations significant at 0.05 were included. All other values masked to NA
#####################################################################################################
coupling_fun <- function(x){
  library(data.table)
  thresh <- length(x)/2
  spis <- x[1:thresh]
  ets <- x[(thresh+1):length(x)]
  sub <- data.table(spi=spis, et=ets)

  # initialize a vector for rho and p-value
  rho <- c()
  pval <- c()
  
  # since we are looking at a 5 year rolling coupling and only looking at growing season months (6 months of year),look at blocks of 30 rows at a time
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
  rho[pval > 0.05] <- NA
  final_rho <- rho
  
  return(final_rho)
}

####################################################################################################
## Function: Calculate Sen's slope and Man Kendal trend test for each row where at least 25% of the correlation coefficeints are significant
## This is used within the rollingCoupling function 
## Inputs:
#### x: a vector of 5-year rolling coupling and overall coupling. each vector corresponds to a cell. 
## Outputs:
#### a vector with the Sen's slope and the Man Kendall significance level (p-value)
#####################################################################################################
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


####################################################################################################
## Function: Calculate the rolling and overall coupling of monthly SPI and monthly ET anomalies and the 
## Sen's slope and Man Kendall trend test of the rolling coupling for each pixel
## Inputs:
#### drought_files: monthly spi raster files 
#### ET_files: monthly ET raster files 
#### dates: A vector of dates corresponding to the monthly drought and ET files in the form YYYYmmdd
#### fmask: forest mask 
#### cl: intialized cluster 
#### sensor: MODIS or Landsat 
#### tl: tile ID or NA for MODIS 
## Outputs:
#### a .csv file where each row corresponds to a pixel, the first column is cell number,
###### and each column is the rolling coupling at each consective time step and the last column is the overall coupling 
#####################################################################################################
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
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/rollingCoupling.csv"))
  }else{
    fwrite(beta_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/rollingCoupling_", tl,".csv"))
  }
  
  
  sens_dt1 <- t(parApply(cl, beta_dt1, FUN=sens_row, MARGIN = 1))
  sens_dt2 <- t(parApply(cl, beta_dt2, FUN=sens_row, MARGIN = 1))
  sens_dt3 <- t(parApply(cl, beta_dt3, FUN=sens_row, MARGIN = 1))
  sens_dt4 <- t(parApply(cl, beta_dt4, FUN=sens_row, MARGIN = 1))
  
  sens_dt <- rbind(sens_dt1, sens_dt2, sens_dt3, sens_dt4)
  sens_dt <- cbind(spi_dt$cellnum, sens_dt)
  colnames(sens_dt) <- c("cellnum", "senSlope", "pvalue")
  
  if(is.na(tl) == T){
    fwrite(sens_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/coupling_sensslope.csv"))
  }else{
    fwrite(sens_dt, paste0(home, "/Analysis/outputs/", sensor, "/rollingCoupling/coupling_sensslope_", tl,".csv"))
  }
}
  

