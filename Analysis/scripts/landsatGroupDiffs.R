##############################################################################
## Assess differences in composite ET drought response according to 
## forest composition, elevation, and height above nearest drainage for Landsat
##############################################################################
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(raster)

home <- "/share/klmarti3/kmcquil/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n = 10 # how many cores do i want 

# start a cluster 
UseCores <- detectCores()-1
cl <- makeCluster(UseCores)
registerDoParallel(cl)

# bring in residual, ETanom, and SPI data tables
# drop rows that are all NA except for cell number 
res_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanResiduals.csv"), nThread = n)
res_dt <- res_dt[!Reduce(`&`, lapply(res_dt[,2:ncol(res_dt)], is.na))]

et_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanET.csv"), nThread = n)
et_dt <- et_dt[!Reduce(`&`, lapply(et_dt[,2:ncol(et_dt)], is.na))]

spi_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanSPI.csv"), nThread = n)
spi_dt <- spi_dt[!Reduce(`&`, lapply(spi_dt[,2:ncol(spi_dt)], is.na))]
spi_dt <- spi_dt[cellnum %in% et_dt$cellnum,]

# first put the mean residduals into a raster so I can visualize 
ExampleRaster <- raster(paste0(home, "/Data/Landsat_ET/landsat_template/landsat_template_ROI.tif"))
sel_cols <- c("cellnum", "meanResid")
meanResid_dt <- res_dt[, ..sel_cols]
meanResid <- DTtoRast(meanResid_dt, ExampleRaster, paste0(home, "/Analysis/outputs/Landsat/compositeRasters/composite_residuals.tif"), open = T)

# now put the mean ET anomalies into a raster for comparison 
sel_cols <- c("cellnum", "meanET")
meanET_dt <- et_dt[, ..sel_cols]
meanET <- DTtoRast(meanET_dt, ExampleRaster, paste0(home, "/Analysis/outputs/Landsat/compositeRasters/composite_ETanom.tif"), open = T)


########################################################################################################################
## Create a DT and raster composite for ET anomalies at different drought severity levels (moderate, severe, extreme)
#########################################################################################################################
etIN <- et_dt[,1:(ncol(et_dt)-1)]
spiINN <- spi_dt[,1:(ncol(spi_dt)-1)]

lowers <- c(-1.6, -2, -5)
uppers <- c(-1.3, -1.6, -2)
et_comps <- stack()

for(i in 1:length(lowers)){
  out <- avgByDroughtSeverity(etIN, spiINN, lowers[i], uppers[i])
  et_dt <- merge(et_dt, out, by = "cellnum", all.x=T, all.y = T)
  
  et_comps <- stack(et_comps, DTtoRast(out, ExampleRaster, paste0(home, "/Analysis/outputs/Landsat/compositeRasters/ETanom_", names(out)[2], ".tif"), open=T))
}


########################################################################################################################
## Create a DT and raster composite for ET residuals at different drought severity levels (moderate, severe, extreme)
#########################################################################################################################
resIN <- res_dt[,1:(ncol(res_dt)-1)]
res_comps <- stack()

for(i in 1:length(lowers)){
  out <- avgByDroughtSeverity(resIN, spiINN, lowers[i], uppers[i])
  res_dt <- merge(res_dt, out, by = "cellnum", all.x=T, all.y = T)
  
  res_comps <- stack(res_comps, DTtoRast(out, ExampleRaster, paste0(home, "/Analysis/outputs/Landsat/compositeRasters/residual_", names(out)[2], ".tif"), open=T))
}


fwrite(et_dt, paste0(home, "/Analysis/outputs/Landsat/meanET.csv"), nThread = n)
fwrite(res_dt, paste0(home, "/Analysis/outputs/Landsat/meanResiduals.csv"), nThread = n)


###################################################################################################################################### 
## Create a datatable of the ET and residual drought composites + FC (Wilson and Riley) + Diffuse porous + Elevation + HAND + cellnum
######################################################################################################################################
et_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanET.csv"), nThread = n)
res_dt <- fread(paste0(home, "/Analysis/outputs/Landsat/meanResiduals.csv"), nThread = n)


# grab the first and last four columns from ET dt and rename 
colnum <- c(1,(ncol(et_dt) - 3):ncol(et_dt))
ET_cols <- et_dt[, ..colnum]
colnames(ET_cols ) <- c("cellnum","ET_all", "ET_moderate", "ET_severe", "ET_extreme")

# grab the last four columns from residuals dt and rename 
Res_cols <- res_dt[, ..colnum]
colnames(Res_cols) <- c("cellnum","Res_all", "Res_moderate", "Res_severe", "Res_extreme")

# merge based on cellnum
data <- merge(ET_cols, Res_cols, by = "cellnum", all.x=T, all.y=T)

# get rid of these large objects that are no longer needed 
rm(et_dt)
rm(res_dt)
rm(ET_cols)
rm(Res_cols)

# grab the forest composition products and merge them in 
FC_wilson <- as.data.frame(raster(paste0(home, "/Data/forest_composition/wilson_landsat_classified.tif")))
setDT(FC_wilson)
cell_col <- data.table(cellnum = seq(1, nrow(FC_wilson)))
FC_wilson <- cbind(cell_col, FC_wilson)
colnames(FC_wilson) <- c("cellnum","FCw")

FC_riley <- as.data.frame(raster(paste0(home, "/Data/forest_composition/riley_landsat_classified.tif")))
setDT(FC_riley)
FC_riley <- cbind(cell_col, FC_riley)
colnames(FC_riley) <- c("cellnum","FCr")

data <- merge(data, FC_wilson, by = "cellnum", all.x=T)
data <- merge(data, FC_riley, by = "cellnum", all.x=T)


# grab the elevation and merge it in 
elevation <- as.data.frame(raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation30m.tif")))
setDT(elevation)
elevation <- cbind(cell_col, elevation)
colnames(elevation) <- c("cellnum", "elevation")
data <- merge(data, elevation, by = "cellnum", all.x=T)

# grab the HAND and merge it in 
hand <- as.data.frame(raster(paste0(home, "/Data/Topography/HeightAboveNearestDrainage/hand_landsat.tif")))
setDT(hand)
hand <- cbind(cell_col, hand)
colnames(hand) <- c("cellnum", "HAND")
data <- merge(data, hand, by = "cellnum", all.x=T)

# Create labels for FC, elevation, and HAND 
#FC first : 1 = diffuse, 2 = ring, 3 = tracheid 
data$FCw <- ifelse(data$FCw == 1, "Diffuse", ifelse(data$FCw == 2, "Ring", ifelse(data$FCw == 3, "Tracheid", NA)))
data$FCr <- ifelse(data$FCr == 1, "Diffuse", ifelse(data$FCr == 2, "Ring", ifelse(data$FCr == 3, "Tracheid", NA)))

# elevation break into high (>= 1200), medium [800 - 1200), and low [0 - 800]
data$Elev_class <- ifelse(data$elevation >=0 & data$elevation < 800, "Low", 
                          ifelse(data$elevation >= 800 & data$elevation < 1200, "Med", 
                                 ifelse(data$elevation >= 1200, "High", NA)))

# for height above nearest drainage, riparian [0-50), slope [50 - 100), ridge [100 - ifn)
data$HAND_class <- ifelse(data$HAND >=0 & data$HAND < 50, "Riparian", 
                          ifelse(data$HAND >= 50 & data$HAND < 100, "Slope", 
                                 ifelse(data$HAND >= 100, "Ridge", NA)))

# now add in percent diffuse for Riley and Wilson 
riley_percent_diffuse <- as.data.frame(raster(paste0(home, "/Data/forest_composition/riley_landsat_diffuse_percent.tif")))
setDT(riley_percent_diffuse)
riley_percent_diffuse <- cbind(cell_col, riley_percent_diffuse)
colnames(riley_percent_diffuse) <- c("cellnum", "pctDiff_R")
data <- merge(data, riley_percent_diffuse, by = "cellnum", all.x=T)

wilson_percent_diffuse <- as.data.frame(raster(paste0(home, "/Data/forest_composition/wilson_landsat_diffuse_percent.tif")))
setDT(wilson_percent_diffuse)
wilson_percent_diffuse <- cbind(cell_col, wilson_percent_diffuse)
colnames(wilson_percent_diffuse) <- c("cellnum", "pctDiff_W")
data <- merge(data, wilson_percent_diffuse, by = "cellnum", all.x=T)


# add in sens slope of coupling 
sens_coupling <- as.data.frame(raster(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/sig_sens_coupling_final.tif")))
setDT(sens_coupling)
sens_coupling <- cbind(cell_col, sens_coupling)
colnames(sens_coupling) <- c("cellnum", "sensCoupling")
data <- merge(data, sens_coupling, by = "cellnum", all.x=T)

# add in overall coupling 
all_coupling <- as.data.frame(raster(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/sig_overall_coupling_final.tif")))
setDT(all_coupling)
all_coupling <- cbind(cell_col, all_coupling)
colnames(all_coupling) <- c("cellnum", "allCoupling")
data <- merge(data, all_coupling, by = "cellnum", all.x=T)

# add in sens slope of ET anomalies 
sens_anoms <- fread(paste0(home, "/Analysis/outputs/Landsat/aggDroughtPeaks/pixelWiseETanomSlopes.csv"))
sens_anoms <- sens_anoms[sens_anoms$V3 <= 0.05,]
sens_anoms <- sens_anoms[,c(1,2)]
colnames(sens_anoms) <- c("cellnum", "sensAnoms")
data <- merge(data, sens_anoms, by = "cellnum", all.x=T)


# add in the forest mask just to make sure 
forest_mask <- as.data.frame(raster(paste0(home, "/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")))
forest_mask <- cbind(cell_col, forest_mask)
colnames(forest_mask) <- c("cellnum", "fmask")
data <- merge(data, fmask, by = "cellnum", all.x=T)

fwrite(data, paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"), nThread = UseCores)


data <- fread(paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"), nThread = UseCores)

# bring in non-drought response 
noDrought <- fread(paste0(home, "/Analysis/outputs/Landsat/avgNonDroughtETResponse.csv"), nThread = UseCores)
colnames(noDrought) <- c("cellnum", 'avgNoDroughtETanom', 'avgNoDroughtETres')
data <- merge(data, noDrought, by = "cellnum", all.x=T)
fwrite(data, paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"), nThread = UseCores)


## Since we decided to use TWI instead of HAND at the last minute after all of the code for the steps after this had been written
## I am just going to put TWI into the HAND column so I don't have to rewrite a bunch of other scripts to match the new variable name 
cores <- detectCores() - 1
data <- fread(paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"), na.strings = (""), nThread = cores)
twi=values(raster("/share/klmarti3/kmcquil/Chapter1_ET_Project/Data/Topography/TWI/TWI_landsat.tif"))
twi <- data.table(twi=twi, 
                  cellnum = seq(1, length(twi), 1))
data[,HAND:=as.numeric(rep(NA, nrow(data)))]
setkey(data, cellnum)
setkey(twi, cellnum)
data[twi, HAND:=twi]

fwrite(data, paste0(home, "/Analysis/outputs/Landsat/final_drought_attribute_dt.csv"))



