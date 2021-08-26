##############################################################################
## Assess differences in composite ET drought response according to 
## forest composition, elevation, and height above nearest drainage for MODIS
##############################################################################

library(data.table)
library(tidyr)
library(ggplot2)
library(viridis)

home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Analysis/scripts/analysis_funcs.R"))
n = 20 # how many cores do i want 

# bring in residual, ETanom, and SPI data tables
res_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanResiduals.csv"))
et_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanET.csv"))
spi_dt <- fread(paste0(home, "/Analysis/outputs/MODIS/meanSPI.csv"))

# first put the mean residduals into a raster so I can visualize 
ExampleRaster <- raster(paste0(home, "/Data/MODIS_ET/clean/et_20000101.tif"))
sel_cols <- c("cellnum", "meanResid")
meanResid_dt <- res_dt[, ..sel_cols]
meanResid <- DTtoRast(meanResid_dt, ExampleRaster, paste0(home, "/Analysis/outputs/MODIS/compositeRasters/composite_residuals.tif"), open = T)

# now put the mean ET anomalies into a raster for comparison 
sel_cols <- c("cellnum", "meanET")
meanET_dt <- et_dt[, ..sel_cols]
meanET <- DTtoRast(meanET_dt, ExampleRaster, paste0(home, "/Analysis/outputs/MODIS/compositeRasters/composite_ETanom.tif"), open = T)


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
  
  et_comps <- stack(et_comps, DTtoRast(out, ExampleRaster, paste0(home, "/Analysis/outputs/MODIS/compositeRasters/ETanom_", names(out)[2], ".tif"), open=T))
}

par(mfrow = c(1,3))
hist(et_dt$`mean_-1.6_-1.3`)
hist(et_dt$`mean_-2_-1.6`)
hist(et_dt$`mean_-5_-2`)


########################################################################################################################
## Create a DT and raster composite for ET residuals at different drought severity levels (moderate, severe, extreme)
#########################################################################################################################
resIN <- res_dt[,1:(ncol(res_dt)-1)]
spiINN <- spi_dt[,1:(ncol(spi_dt)-1)]

lowers <- c(-1.6, -2, -5)
uppers <- c(-1.3, -1.6, -2)
res_comps <- stack()

for(i in 1:length(lowers)){
  out <- avgByDroughtSeverity(resIN, spiINN, lowers[i], uppers[i])
  res_dt <- merge(res_dt, out, by = "cellnum", all.x=T, all.y = T)
  
  res_comps <- stack(res_comps, DTtoRast(out, ExampleRaster, paste0(home, "/Analysis/outputs/MODIS/compositeRasters/residual_", names(out)[2], ".tif"), open=T))
}


####################################################################################################################
## Create a datatable of the ET and residual drought composites + FC (Wilson and Riley) + Elevation + HAND + cellnum
####################################################################################################################

# grab the first and last four columns from ET dt and rename 
colnum <- c(1,(ncol(et_dt) - 3):ncol(et_dt))
ET_cols <- et_dt[, ..colnum]
colnames(ET_cols ) <- c("cellnum","ET_all", "ET_moderate", "ET_severe", "ET_extreme")

# grab the last four columns from residuals dt and rename 
Res_cols <- res_dt[, ..colnum]
colnames(Res_cols) <- c("cellnum","Res_all", "Res_moderate", "Res_severe", "Res_extreme")

# merge based on cellnum
data <- merge(ET_cols, Res_cols, by = "cellnum", all.x=T, all.y=T)

# grab the forest composition products and merge them in 
FC_wilson <- as.data.frame(raster(paste0(home, "/Data/forest_composition/wilson_modis_classified.tif")))
setDT(FC_wilson)
cell_col <- data.table(cellnum = seq(1, nrow(FC_wilson)))
FC_wilson <- cbind(cell_col, FC_wilson)
colnames(FC_wilson) <- c("cellnum","FCw")

FC_riley <- as.data.frame(raster(paste0(home, "/Data/forest_composition/riley_modis_classified.tif")))
setDT(FC_riley)
FC_riley <- cbind(cell_col, FC_riley)
colnames(FC_riley) <- c("cellnum","FCr")

data <- merge(data, FC_wilson, by = "cellnum", all.x=T)
data <- merge(data, FC_riley, by = "cellnum", all.x=T)


# grab the elevation and merge it in 
elevation <- as.data.frame(raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation500m.tif")))
setDT(elevation)
elevation <- cbind(cell_col, elevation)
colnames(elevation) <- c("cellnum", "elevation")
data <- merge(data, elevation, by = "cellnum", all.x=T)

# grab the HAND and merge it in 
hand <- as.data.frame(raster(paste0(home, "/Data/Topography/HeightAboveNearestDrainage/hand_modis.tif")))
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
riley_percent_diffuse <- as.data.frame(raster(paste0(home, "/Data/forest_composition/riley_MODIS_diffuse_percent.tif")))
setDT(riley_percent_diffuse)
riley_percent_diffuse <- cbind(cell_col, riley_percent_diffuse)
colnames(riley_percent_diffuse) <- c("cellnum", "pctDiff_R")
data <- merge(data, riley_percent_diffuse, by = "cellnum", all.x=T)

wilson_percent_diffuse <- as.data.frame(raster(paste0(home, "/Data/forest_composition/wilson_MODIS_diffuse_percent.tif")))
setDT(wilson_percent_diffuse)
wilson_percent_diffuse <- cbind(cell_col, wilson_percent_diffuse)
colnames(wilson_percent_diffuse) <- c("cellnum", "pctDiff_W")
data <- merge(data, wilson_percent_diffuse, by = "cellnum", all.x=T)



# now add in the average ET anomaly and residual during non drought periods 
noDrought <- fread(paste0(home, "/Analysis/outputs/MODIS/avgNonDroughtETResponse.csv"))
colnames(noDrought) <- c("cellnum", 'avgNoDroughtETanom', 'avgNoDroughtETres')
data <- merge(data, noDrought, by = "cellnum", all.x=T)

## this is the data.table to work with from here on out to subset and look at significnat differences between groups! 
fwrite(data, paste0(home, "/Analysis/outputs/MODIS/final_drought_attribute_dt.csv"))


################################################################################################################
## Assess differences between groups: 
#################################################################################################################

# bring in the datatable with drought composite info and fc and topo attributes
data <- fread(paste0(home, "/Analysis/outputs/MODIS/final_drought_attribute_dt.csv"), na.strings = (""))

# factorize the Fc, hand classifications, and elevation classifications 
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)             
data$HAND_class <- as.factor(data$HAND_class)  
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))



## Group by just forest composition 
FC_diff <- kwTestbyDroughtGroup(copy(data), "FCr", "FC")

## Are there differences between FC at high elevation, medium elevation, or low elevation? 
FC_high_diff <- kwTestbyDroughtGroup(copy(data[Elev_class == "High"]), "FCr", "FC")
FC_med_diff <- kwTestbyDroughtGroup(copy(data[Elev_class == "Med"]), "FCr", "FC")
FC_low_diff <- kwTestbyDroughtGroup(copy(data[Elev_class == "Low"]), "FCr", "FC")


## Group by just elevation
Elev_diff <- kwTestbyDroughtGroup(copy(data), "Elev_class", "Elevation")

## Are there differrences between one forest composition type across elevation categories? 
Elev_Diffuse_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Diffuse"]), "Elev_class", "Elevation")
Elev_Ring_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Ring"]), "Elev_class", "Elevation")
Elev_Tracheid_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Tracheid"]), "Elev_class", "Elevation")


## group by just HAND
Hand_diff <- kwTestbyDroughtGroup(copy(data), "HAND_class", "HAND")

## Are there differences between FC at riparian, slope, or ridge (HAND)?
FC_riparian_diff <- kwTestbyDroughtGroup(copy(data[HAND_class == "Riparian"]), "FCr", "FC")
FC_slope_diff <- kwTestbyDroughtGroup(copy(data[HAND_class == "Slope"]), "FCr", "FC")
FC_ridge_diff <- kwTestbyDroughtGroup(copy(data[HAND_class == "Ridge"]), "FCr", "FC")


## Are there differences between one forest composition type across HAND categories? 
Hand_Diffuse_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Diffuse"]), "HAND_class", "HAND")
Hand_Ring_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Ring"]), "HAND_class", "HAND")
Hand_Tracheid_diff <- kwTestbyDroughtGroup(copy(data[FCr == "Tracheid"]), "HAND_class", "HAND")





















