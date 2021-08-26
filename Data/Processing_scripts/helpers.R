##########################################################################
##########################################################################
## Random functions that are helpful for geoprocessing or analysis 
##########################################################################
##########################################################################


## resample a raster that's already read in using gdalWarp 
## provide another raster for the paramters you want to resample to 
warp <- function(fin, template, fout, roi, nodata){
  
  res <- res(template)
  t1 <- c(xmin(template), ymin(template), 
          xmax(template), ymax(template))
  res_out <- crs(template)
  
  gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r='near', 
           cutline = roi, 
           crop_to_cutline = T, 
           dstnodata = nodata,
           overwrite = T)
}


warpM <- function(fin, template, fout, roi, nodata, re){
  
  res <- res(template)
  t1 <- c(xmin(template), ymin(template), 
          xmax(template), ymax(template))
  res_out <- crs(template)
  
  gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r=re, 
           cutline = roi, 
           crop_to_cutline = T, 
           dstnodata = nodata,
           overwrite = T)
}

warpMn <- function(fin, template, fout,re){
  
  res <- res(template)
  t1 <- c(xmin(template), ymin(template), 
          xmax(template), ymax(template))
  res_out <- crs(template)
  
  gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r=re, overwrite = T)
}





## function to calcualte the long term mean and sd of a bunch of rasters 
## input is the tile ID 
LTs <- function(TD){
  tile_dt <- file_dt[tile == TD]
  
  
  for(i in 1:length(unique_months)){
    sub <- tile_dt[month == unique_months[i]]$file
    stk <- do.call('stack', lapply(sub, raster))
    print("finished stacking")
    print(nlayers(stk))
    avg_stk <- clusterR(stk, fun = calc, args =list(fun = mean, na.rm = T))
    sd_stk <- clusterR(stk, fun = calc, args = list(fun=sd, na.rm=T))
    
    writeRaster(sd_stk, paste0(home, "/Landsat_ET/ltsd_month/", TD, "_",sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    
    writeRaster(avg_stk, paste0(home, "/Landsat_ET/ltm_month/et_",TD, "_", sprintf("%02d", unique_months[i]), ".tif"), 
                overwrite = T, format = "GTiff")
    print(i)
  }
  
  
}



