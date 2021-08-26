library(data.table)
library(trend)
library(tidyr)
library(ggplot2)
library(doParallel)
library(parallel)
library(foreach)

home <- "G:/My Drive/Chapter1_ET_Project"
home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
# how many cores do i want 

UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)
print(cl)

# bring in residual,and SPI data tables
res_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanResiduals.csv"), nThread = UseCores)
res_dt <- res_dt[!Reduce(`&`, lapply(res_dt[,2:ncol(res_dt)], is.na))]

spi_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanSPI.csv"), nThread = UseCores)
spi_dt <- spi_dt[!Reduce(`&`, lapply(spi_dt[,2:ncol(spi_dt)], is.na))]
spi_dt <- spi_dt[cellnum %in% res_dt$cellnum,]


res_dt_long <- setDT(gather(res_dt[,1:(ncol(res_dt)-1)], Year, ET_Residual, `1985`:`2019`))
spi_dt_long <- setDT(gather(spi_dt[,1:(ncol(spi_dt)-1)], Year, SPI, `1985`:`2019`))
dim(res_dt_long)

rm(res_dt)
rm(spi_dt)

cols <- c("cellnum", "Year")
setkeyv(res_dt_long, cols)
setkeyv(spi_dt_long, cols)
dt_long <- res_dt_long[spi_dt_long, nomatch = 0]

rm(spi_dt_long)
rm(res_dt_long)

dt_long <- dt_long[complete.cases(dt_long),]
dt_long$Year <- as.numeric(dt_long$Year)
dim(dt_long)

## get the data tables of summarized info for all droughts for landsat
sens_fun_landsat(copy(dt_long), "all")

dt_long_mod <- copy(dt_long[SPI <= -1.3 & SPI > -1.6, ])
dim(dt_long_mod)
sens_fun_landsat(dt_long_mod, "moderate")

dt_long_severe <- copy(dt_long[SPI <= -1.6 & SPI > -2, ])
dim(dt_long_severe)
sens_fun_landsat(dt_long_severe, "severe")


dt_long_extreme <- copy(dt_long[SPI <= -2,])
dim(dt_long_extreme)
sens_fun_landsat(dt_long_extreme, "extreme")


# make figures on personal computer 
home <- "G:/My Drive/Chapter1_ET_Project"
sens_all <- sens_fun_landsat_figs(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/annual_all.csv"),
                       paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling3_all.csv"), 
                       paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling5_all.csv"))
 
sens_mod <- sens_fun_landsat_figs(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/annual_moderate.csv"),
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling3_moderate.csv"), 
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling5_moderate.csv"))
 
 
sens_severe <- sens_fun_landsat_figs(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/annual_severe.csv"),
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling3_severe.csv"), 
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling5_severe.csv"))
 
sens_extreme <- sens_fun_landsat_figs(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/annual_extreme.csv"),
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling3_extreme.csv"), 
                                   paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/rolling5_extreme.csv"))



######################################################################################
######################################################################################
# Now get the pixel-wise trends in drought ET residuals 
# only look at pixels that have at least 5 droughts to calculate a trend from 
dim(dt_long)
countDroughtInCell <- dt_long[!is.na(SPI), .N, .(cellnum)]
cellsOfInterest <- countDroughtInCell[N >=5]$cellnum
dt_long_sub <- dt_long[cellnum %in% cellsOfInterest, .(cellnum, Year, ET_Residual)]
dim(dt_long_sub)
rm(dt_long)

# put it back the wide way 
dt_wide_sub <- spread(dt_long_sub, Year, ET_Residual)
rm(dt_long_sub)
# make sure there are no rows that are all NA except for the cellnum 
idx <- apply(dt_wide_sub[,2:ncol(dt_wide_sub)], 1, FUN = function(x) { sum(is.na(x))})
dt_wide_sub <- dt_wide_sub[!idx== (ncol(dt_wide_sub) - 1),]


lm_fun <- function(x){
  library(trend)
  cellnum <- x[1]
  x <- as.numeric(x[2:length(x)])
  x <- x[!is.na(x)]
  s <- sens.slope(x)
  ss <- c(cellnum, s$estimates[1], s$p.value[1])
  
}

cellWiseSense <- parApply(cl, dt_wide_sub, MARGIN = 1, FUN =lm_fun)
cellWiseSense <- t(cellWiseSense)
fwrite(cellWiseSense, paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/pixelWiseETanomSlopes.csv"), nThread = UseCores)


slopes <- fread(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/pixelWiseETanomSlopes.csv"), nThread = UseCores)
slopes <- slopes[slopes$V3 <= 0.05,]
slopes <- slopes[,c(1,2)]
template <- raster(paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
out <- paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/pixelWiseETanomSlopes_raster.tif")
DTtoRast(slopes, template, out, open = T)
slopes_rast <- raster(out)

hist(values(slopes_rast), xlab = "Sen's slope (p<=0.05) for Et residuals at drought peak ")







