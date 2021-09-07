##########################################################################
##########################################################################
## Functions used in data preparation 
##########################################################################
##########################################################################

###########################################################################
## Resample a raster to match another rasters resolution, extent, and crs
## Input: 
#### fin: the filepath of the raster to be resampled 
#### template: the raster (read in) to serve as the template 
#### fout: the filepath to write the new resampled raster to 
#### re: resample method to be used. Use the same names as rgdal
## Output: 
#### a resampled raster that matches the template 
#############################################################################
warpMn <- function(fin, template, fout,re){
  
  res <- res(template)
  t1 <- c(xmin(template), ymin(template), 
          xmax(template), ymax(template))
  res_out <- crs(template)
  
  gdalwarp(fin, fout, t_srs=res_out, tr=res, te=t1, r=re, overwrite = T)
}


###########################################################################
## Calculate the monthly long term mean and standard deviation on a specified raster tile
## Input: 
#### TD: the tile ID 
## Output: 
##### Write out the long term mean raster for each month for the given tile
##### Write out the long term standard deviation raster for each month for the given tile 
#############################################################################
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



