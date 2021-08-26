#############################################################################################
## This is a scratch pad to make figures and do small calculations of MODIS/Landsat 
##############################################################################################
library(rasterVis)
library(raster)
library(pals)
library(RColorBrewer)
library(trend)
library(data.table)
library(viridis)
library(ggplot2)
home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))


## Calculate the correlation between elevation, HAND, and % diffuse porous for MODIS and Landsat 
data <- fread(paste0(home, "/Analysis/outputs/MODIS/final_drought_attribute_dt.csv"), na.strings = (""))
# factorize the Fc, hand classifications, and elevation classifications 
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)             
data$HAND_class <- as.factor(data$HAND_class)  
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))

preds <- data[,.(elevation, HAND, pctDiff_R)]
modis_corr <- cor(preds, use="pairwise.complete.obs")




# create a figure of elevation, HAND, and % diffuse porous maps next to each other 
elevation <- raster(paste0(home, "/Data/Topography/usgsNED_elevation/elevation30m.tif"))
hand <- raster(paste0(home, "/Data/Topography/HeightAboveNearestDrainage/hand_landsat.tif"))
pctDiffR <- raster(paste0(home, "/Data/forest_composition/riley_landsat_diffuse_percent.tif"))
pctDiffR <- pctDiffR*100


par(mfrow=c(1, 1), mai=c(0.5, 6, 0.5, 0),omi=c(0,0.8,0,0), new=FALSE)
plot(anom_stk_m[[1]], legend.only=TRUE, legend.shrink=1, legend.width=1, 
     zlim=c(quants[[1]], quants[[2]]),
     col=col5(n=color_levels), breaks=color_sequence,
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), font = 2,labels=pretty(quants[[1]]:quants[[2]])),
     legend.args=list(text='', side=4, font=2, line=2.3))


par(mfrow = c(1,3))
par(mfrow = c(1,3), mai = c(0, 0.25, 0, 0), omi=c(0,0.1,0,0.25))
plot(elevation, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))
plot(hand, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))
plot(pctDiffR, col = viridis(100), axes=F, box = F, legend.width = 1.5, axis.args = list(font = 2))


# create plots combining MODIS and Landsat looking at the percent of anomalies that are greater than 0 
# get the modis files 
modis_files <- list.files(paste0(home, "/Visualizing/MODIS/"), full.names=T, pattern = ".csv$")
landsat_files <- list.files(paste0(home, "/Visualizing/Landsat/"), full.names=T, pattern = ".csv$")
landsat_files <- landsat_files[grep("pct_greater0", landsat_files)]

# bring together the first file for landsat and modis to plot the % of anomalies greater than 0 
pg0_l <- fread(landsat_files[1])
pg0_l$sensor <- rep("Landsat", 5)
pg0_m <- fread(modis_files[1])
pg0_m$sensor <- rep("MODIS", 5)
pg0 <- setDT(rbind(pg0_l, pg0_m))
pg0$DroughtSeverity <- factor(pg0$DroughtSeverity, ordered=TRUE, 
                                 levels = c("Moderate", "Severe", "Extreme", "All", "None"))
pg0_sub <- pg0[DroughtSeverity == "Moderate" |DroughtSeverity == "Severe" | DroughtSeverity == "Extreme", ]
pg0_sub$pg0 <- pg0_sub$pg0*100
ggplot(pg0_sub, aes(x = DroughtSeverity, y = pg0, fill = sensor)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic()+
  xlab("Drought Severity") + 
  ylab(expression(bold(`%`~ET["dp"]~`>`~0))) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=12))+
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65))
  
# bring together the elevation files 
elev_l <- fread(landsat_files[2])
elev_l$sensor <- rep("Landsat", nrow(elev_l))
elev_m <- fread(modis_files[2])
elev_m$sensor <- rep("MODIS", nrow(elev_m))

elev <- setDT(rbind(elev_l, elev_m))
elev$DroughtSeverity <- factor(elev$DroughtSeverity, ordered=TRUE, 
                               levels = c("Moderate", "Severe", "Extreme", "All", "None"))

elev <- elev[DroughtSeverity == "Moderate" |DroughtSeverity == "Severe" | DroughtSeverity == "Extreme", ]
elev_modis <- elev[sensor == "MODIS",]

elev_modis_plot <- ggplot(elev_modis, aes(x = elevation_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned Elevation") + 
  ylab("% Anomalies > 0") + 
  ggtitle("MODIS")+
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))


 
elev_landsat <- elev[sensor == "Landsat",]
elev_landsat_plot <- ggplot(elev_landsat, aes(x = elevation_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned Elevation") + 
  ylab("% Anomalies > 0") + 
  ggtitle("Landsat")+
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

  
# bring together the hand files 
hand_l <- fread(landsat_files[3])
hand_l$sensor <- rep("Landsat", nrow(hand_l))
hand_m <- fread(modis_files[3])
hand_m$sensor <- rep("MODIS", nrow(hand_m))

hand <- setDT(rbind(hand_l, hand_m))
hand$DroughtSeverity <- factor(hand$DroughtSeverity, ordered=TRUE, 
                               levels = c("Moderate", "Severe", "Extreme", "All", "None"))

hand <- hand[DroughtSeverity == "Moderate" |DroughtSeverity == "Severe" | DroughtSeverity == "Extreme", ]
hand_modis <- hand[sensor == "MODIS",]

hand_modis_plot <- ggplot(hand_modis, aes(x = HAND_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned HAND") + 
  ylab("% Anomalies > 0") + 
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

hand_landsat <- hand[sensor == "Landsat",]
hand_landsat_plot <- ggplot(hand_landsat, aes(x = HAND_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned HAND") + 
  ylab("% Anomalies > 0") + 
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))



# bring together the percent diffuse porous files 
diff_l <- fread(landsat_files[4])
diff_l$sensor <- rep("Landsat", nrow(diff_l))
diff_m <- fread(modis_files[4])
diff_m$sensor <- rep("MODIS", nrow(diff_m))

diff <- setDT(rbind(diff_l, diff_m))
diff$DroughtSeverity <- factor(diff$DroughtSeverity, ordered=TRUE, 
                               levels = c("Moderate", "Severe", "Extreme", "All", "None"))

diff <- diff[DroughtSeverity == "Moderate" |DroughtSeverity == "Severe" | DroughtSeverity == "Extreme", ]
diff_modis <- diff[sensor == "MODIS",]

diff_modis_plot <-ggplot(diff_modis, aes(x = pctDiff_R_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned % Diffuse Porous") + 
  ylab("% Anomalies > 0") + 
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

diff_landsat <- diff[sensor == "Landsat",]
diff_landsat_plot <- ggplot(diff_landsat, aes(x = pctDiff_R_bin, y = pg0, fill = DroughtSeverity)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_color_manual(values = viridis(3))+
  scale_fill_manual(values = viridis(3))+
  theme_bw()+
  xlab("Binned % Diffuse Porous") + 
  ylab("% Anomalies > 0") + 
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

  
  
  
  
library(cowplot) 
plot_grid(elev_modis_plot, elev_landsat_plot,
          hand_modis_plot, hand_landsat_plot, 
          diff_modis_plot, diff_landsat_plot, labels = "AUTO", nrow = 3)
  
  





###################################################################################
## create a bivariate map of coupling and sen's slope 
####################################################################################
library(classInt)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
correlation_modis <- raster(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_overallCoupling.tif"))
trend_modis <- raster(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/BIV_slope.tif"))
correlation_landsat <- raster(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/BIV_sig_overall_coupling_final.tif"))
trend_landsat <- raster(paste0(home, "/Analysis/outputs/Landsat/rollingCoupling/BIV_sig_sens_coupling_final.tif"))

nquantiles = 10
cor_mod_quants <- quantile(values(correlation_modis), na.rm=TRUE, probs = c(seq(0,1,1/nquantiles)))
cor_landsat_quants <- quantile(values(correlation_landsat), na.rm=TRUE, probs = c(seq(0,1,1/nquantiles)))
# use the landsat corerlation quantiles 

trend_mod_quants <- quantile(values(trend_modis), na.rm=TRUE, probs = c(seq(0,1,1/nquantiles)))
trend_landsat_quants <- quantile(values(trend_landsat), na.rm=TRUE, probs = c(seq(0,1,1/nquantiles)))
# use the landsat trend quantiles 


colmat.custom<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label", x_quants, y_quants){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3,  xaxt = "n", yaxt="n")
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  axis(1, at=seq(0.0, 1, 0.1), labels = round(x_quants, 2))
  axis(2, at=seq(0.0, 1, 0.1), labels = round(y_quants, 3))
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]
}


bivariate.map.flex<-function(rasterx, rastery, colormatrix=col.matrix, brks1, brks2){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks1, labels = 2:length(brks1),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks2, labels = 2:length(brks2),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)
}



col.matrix<-colmat.custom(nquantiles=10, xlab = "R", ylab = "Trend in R",bottomleft="grey", x_quants = cor_landsat_quants, y_quants = trend_landsat_quants)
bivmap_modis<-bivariate.map.flex(correlation_modis, trend_modis, colormatrix=col.matrix, brks1 = cor_landsat_quants, brks2 = trend_landsat_quants)
bivmap_landsat<-bivariate.map.flex(correlation_landsat, trend_landsat, colormatrix=col.matrix, brks1 = cor_landsat_quants, brks2 = trend_landsat_quants)

sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")
par(mfrow= c(1,2), omi=c(0.1,0,0,0))
plot(bivmap_modis,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
plot(sbr, add = T)
plot(bivmap_landsat,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
plot(sbr, add = T)









##########################################################################################
## make plots of the correlation and the trend that combine the MODIS and Landsat 
##########################################################################################
dir <- "G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

# bring in the datatable with drought composite info and fc and topo attributes
data <- fread(paste0(home, "/Analysis/outputs/MODIS/final_drought_attribute_dt.csv"), na.strings = (""))
data$FCw <- as.factor(data$FCw)
data$FCr <- as.factor(data$FCr)             
data$HAND_class <- as.factor(data$HAND_class)  
data$HAND_class <- ordered(data$HAND_class, levels = c("Riparian", "Slope", "Ridge"))
data$Elev_class <- as.factor(data$Elev_class)
data$Elev_class <- ordered(data$Elev_class, levels = c("Low", "Med", "High"))
pd_bins <- seq(0.05, 0.85, 0.05)
pd_names <- pd_bins[1:length(pd_bins)-1]*100
data$pctDiff_R_bin <- cut(data$pctDiff_R, breaks = pd_bins, labels = pd_names)
data$pctDiff_W_bin <- cut(data$pctDiff_W, breaks = pd_bins, labels = pd_names)
e_bins <- seq(200, 1950, 50)
e_names <- e_bins[1:length(e_bins)-1]
data$elevation_bin <- cut(data$elevation, breaks = e_bins, labels = e_names)
h_bins<- seq(0, 270, 10)
h_names <- h_bins[1:length(h_bins)-1]
data$HAND_bin <- cut(data$HAND, breaks = h_bins, labels = h_names)
forest_mask <- as.data.frame(raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/MODIS_FOREST/modis_permanent_forest_resampled.tif"))
forest_mask <- cbind(seq(1, nrow(forest_mask)), forest_mask)
colnames(forest_mask) <- c("cellnum", "fmask")
forest_mask1 <- forest_mask[forest_mask$fmask== 1 & !is.na(forest_mask$fmask),]
coupling_slope <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/coupling_sensslope.csv"))
coupling_slope <- coupling_slope[coupling_slope$cellnum %in% forest_mask1$cellnum,]
coupling_slope_sig <- coupling_slope[coupling_slope$pvalue < 0.05,]
rollingCoupling <- fread(paste0(home, "/Analysis/outputs/MODIS/rollingCoupling/rollingCoupling.csv"))
rollingCoupling <- rollingCoupling[rollingCoupling$cellnum %in% forest_mask1$cellnum,]
data <- merge(data, coupling_slope_sig[,c("cellnum", "senSlope")], by = "cellnum", all.x=T)
data <- merge(data, rollingCoupling[, c("cellnum", "all")], by = "cellnum", all.x=T)

agg_data <- function(DF, X, Y){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  agg <- agg[order(agg[[1]]),]
  # anom_sen <- sens.slope(agg$meanAnom)
  # results <- data.table(VAR = X, 
  #                       senslope = anom_sen$estimates, 
  #                       pval = anom_sen$p.value)
  return(agg)
}


############################################################################################3
# Start with overall coupling 
### first do elevation
m_elev <- agg_data(data, "elevation_bin", "all")
m_elev$Sensor <- rep("MODIS", nrow(m_elev))
l_elev <- fread(paste0(dir, "/elevation_allCoupling.csv"))
l_elev$Sensor <- rep("Landsat", nrow(l_elev))
combo_elev <- setDT(rbind(m_elev, l_elev))
combo_elev$elevation_bin <- as.numeric(as.character(combo_elev$elevation_bin))

m_elev$elevation_bin <- as.numeric(as.character(m_elev$elevation_bin))
m_elev <- m_elev[order(m_elev$elevation_bin),]
m_elev_sens <- sens.slope(m_elev$meanAnom)
m_elev_label <- paste0("MODIS slope = ", round(m_elev_sens$estimates, 4), ifelse(m_elev_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_elev_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_elev_sens$p.value <= 0.05, "*", ""))))
l_elev$elevation_bin <- as.numeric(as.character(l_elev$elevation_bin))
l_elev <- l_elev[order(l_elev$elevation_bin),]
l_elev_sens <- sens.slope(l_elev$meanAnom)
l_elev_label <- paste0("Landsat slope = ", round(l_elev_sens$estimates, 4), ifelse(l_elev_sens$p.value <= 0.001, "***",
                                                                                        ifelse(l_elev_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(l_elev_sens$p.value <= 0.05, "*", ""))))

elev_correlation <- ggplot(data = combo_elev, aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("Elevation (m)") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 725, y = c(-0.3, -0.45), label = c(m_elev_label, l_elev_label))
  

### now do HAND 
m_hand <- agg_data(data, "HAND_bin", "all")
m_hand$Sensor <- rep("MODIS", nrow(m_hand))
l_hand <- fread(paste0(dir, "/HAND_allCoupling.csv"))
l_hand$Sensor <- rep("Landsat", nrow(l_hand))
combo_hand <- setDT(rbind(m_hand, l_hand))
combo_hand$HAND_bin <- as.numeric(as.character(combo_hand$HAND_bin))

m_hand$HAND_bin <- as.numeric(as.character(m_hand$HAND_bin))
m_hand <- m_hand[order(m_hand$HAND_bin),]
m_hand_sens <- sens.slope(m_hand$meanAnom)
m_hand_label <- paste0("MODIS Slope = ", round(m_hand_sens$estimates, 4), ifelse(m_hand_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_hand_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_hand_sens$p.value <= 0.05, "*", ""))))
l_hand$HAND_bin <- as.numeric(as.character(l_hand$HAND_bin))
l_hand <- l_hand[order(l_hand$HAND_bin),]
l_hand_sens <- sens.slope(l_hand$meanAnom)
l_hand_label <- paste0("Landsat Slope = ", round(l_hand_sens$estimates, 4), ifelse(l_hand_sens$p.value <= 0.001, "***",
                                                                                        ifelse(l_hand_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(l_hand_sens$p.value <= 0.05, "*", ""))))

hand_correlation <- ggplot(data = combo_hand, aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("HAND (m)") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 125, y = c(-0.17, -0.27), label = c(m_hand_label, l_hand_label))

### now do percent diffuse porous 
m_diff <- agg_data(data, "pctDiff_R_bin", "all")
m_diff$Sensor <- rep("MODIS", nrow(m_diff))
l_diff <- fread(paste0(dir, "/pctDiffR_allCoupling.csv"))
l_diff$Sensor <- rep("Landsat", nrow(l_diff))
combo_diff <- setDT(rbind(m_diff, l_diff))
combo_diff$pctDiff_R_bin<- as.numeric(as.character(combo_diff$pctDiff_R_bin))

m_diff$pctDiff_R_bin <- as.numeric(as.character(m_diff$pctDiff_R_bin))
m_diff <- m_diff[order(m_diff$pctDiff_R_bin),]
m_diff_sens <- sens.slope(m_diff$meanAnom)
m_diff_label <- paste0("MODIS Slope = ", round(m_diff_sens$estimates, 4), ifelse(m_diff_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_diff_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_diff_sens$p.value <= 0.05, "*", ""))))

l_diff$pctDiff_R_bin <- as.numeric(as.character(l_diff$pctDiff_R_bin))
l_diff <- l_diff[order(l_diff$pctDiff_R_bin),]
l_diff_sens <- sens.slope(l_diff$meanAnom)
l_diff_label <- paste0("Landsat Slope = ", round(l_diff_sens$estimates, 4), ifelse(l_diff_sens$p.value <= 0.001, "***",
                                                                                        ifelse(l_diff_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(l_diff_sens$p.value <= 0.05, "*", ""))))

diff_correlation <- ggplot(data = combo_diff, aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("% Diffuse Porous") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.05, -0.12), label = c(m_diff_label, l_diff_label))


############################################################################################3
# Now do the sen's slope 
### first do elevation
m_elev <- agg_data(data, "elevation_bin", "senSlope")
m_elev$Sensor <- rep("MODIS", nrow(m_elev))
l_elev <- fread(paste0(dir, "/elevation_sensCoupling.csv"))
l_elev$Sensor <- rep("Landsat", nrow(l_elev))
combo_elev <- setDT(rbind(m_elev, l_elev))
combo_elev$elevation_bin <- as.numeric(as.character(combo_elev$elevation_bin))

m_elev$elevation_bin <- as.numeric(as.character(m_elev$elevation_bin))
m_elev <- m_elev[order(m_elev$elevation_bin),]
m_elev_sens <- sens.slope(m_elev$meanAnom)
m_elev_label <- paste0("MODIS Slope = ", round(m_elev_sens$estimates, 4), ifelse(m_elev_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_elev_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_elev_sens$p.value <= 0.05, "*", ""))))
l_elev$elevation_bin <- as.numeric(as.character(l_elev$elevation_bin))
l_elev <- l_elev[order(l_elev$elevation_bin),]
l_elev_sens <- sens.slope(l_elev$meanAnom)
l_elev_label <- paste0("Landsat Slope = ", round(l_elev_sens$estimates, 4), ifelse(l_elev_sens$p.value <= 0.001, "***",
                                                                                           ifelse(l_elev_sens$p.value <= 0.01, "**", 
                                                                                                  ifelse(l_elev_sens$p.value <= 0.05, "*", ""))))

elev_trend <- ggplot(data = combo_elev, aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("Elevation (m)") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 600, y = c(-0.02, -0.028), label = c(m_elev_label, l_elev_label))



### next do HAND 
m_hand <- agg_data(data, "HAND_bin", "senSlope")
m_hand$Sensor <- rep("MODIS", nrow(m_hand))
l_hand <- fread(paste0(dir, "/HAND_sensCoupling.csv"))
l_hand$Sensor <- rep("Landsat", nrow(l_hand))
combo_hand <- setDT(rbind(m_hand, l_hand))
combo_hand$HAND_bin <- as.numeric(as.character(combo_hand$HAND_bin))

m_hand$HAND_bin <- as.numeric(as.character(m_hand$HAND_bin))
m_hand <- m_hand[order(m_hand$HAND_bin),]
m_hand_sens <- sens.slope(m_hand$meanAnom)
m_hand_label <- paste0("MODIS Slope = ", round(m_hand_sens$estimates, 4), ifelse(m_hand_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_hand_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_hand_sens$p.value <= 0.05, "*", ""))))
l_hand$HAND_bin <- as.numeric(as.character(l_hand$HAND_bin))
l_hand <- l_hand[order(l_hand$HAND_bin),]
l_hand_sens <- sens.slope(l_hand$meanAnom)
l_hand_label <- paste0("Landsat Slope = ", round(l_hand_sens$estimates, 4), ifelse(l_hand_sens$p.value <= 0.001, "***",
                                                                                        ifelse(l_hand_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(l_hand_sens$p.value <= 0.05, "*", ""))))

hand_trend <- ggplot(data = combo_hand, aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("HAND (m)") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 100, y = c(-0.022, -0.028), label = c(m_hand_label, l_hand_label))


### next do % diffuse porous 
m_diff <- agg_data(data, "pctDiff_R_bin", "senSlope")
m_diff$Sensor <- rep("MODIS", nrow(m_diff))
l_diff <- fread(paste0(dir, "/pctDiffR_sensCoupling.csv"))
l_diff$Sensor <- rep("Landsat", nrow(l_diff))
combo_diff <- setDT(rbind(m_diff, l_diff))
combo_diff$pctDiff_R_bin<- as.numeric(as.character(combo_diff$pctDiff_R_bin))

m_diff$pctDiff_R_bin <- as.numeric(as.character(m_diff$pctDiff_R_bin))
m_diff <- m_diff[order(m_diff$pctDiff_R_bin),]
m_diff_sens <- sens.slope(m_diff$meanAnom)
m_diff_label <- paste0("MODIS Slope = ", round(m_diff_sens$estimates, 4), ifelse(m_diff_sens$p.value <= 0.001, "***",
                                                                                        ifelse(m_diff_sens$p.value <= 0.01, "**", 
                                                                                               ifelse(m_diff_sens$p.value <= 0.05, "*", ""))))

l_diff$pctDiff_R_bin <- as.numeric(as.character(l_diff$pctDiff_R_bin))
l_diff <- l_diff[order(l_diff$pctDiff_R_bin),]
l_diff_sens <- sens.slope(l_diff$meanAnom)
l_diff_label <- paste0("Landsat Slope = ", round(l_diff_sens$estimates, 4), ifelse(l_diff_sens$p.value <= 0.001, "***",
                                                                                          ifelse(l_diff_sens$p.value <= 0.01, "**", 
                                                                                                 ifelse(l_diff_sens$p.value <= 0.05, "*", ""))))

diff_trend <- ggplot(data = combo_diff, aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("%Diffuse Porous") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.023, -0.029), label = c(m_diff_label, l_diff_label))

library(cowplot)
plot_grid(elev_correlation, elev_trend, 
          hand_correlation, hand_trend, 
          diff_correlation, diff_trend, 
          labels = "AUTO", nrow = 3)








##################################################################################################3
## plotting the mapped anomalies for MODIS andn Landsat together with one colorbar 
#### make maps of the rasters for ET anomalies and ET residuals at each level of drought 
et_tifs_m <- list.files(paste0(home, "/Analysis/outputs/MODIS/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs_m <- et_tifs_m[grep("ETanom", et_tifs_m)]
anom_stk_m <- do.call("stack", lapply(anom_tifs_m, raster))

et_tifs_l <- list.files(paste0(home, "/Analysis/outputs/Landsat/compositeRasters"), full.names = T, pattern = ".tif$")
anom_tifs_l <- et_tifs_l[grep("ETanom", et_tifs_l)]
anom_stk_l <- do.call("stack", lapply(anom_tifs_l, raster))
forest_mask_landsat <- raster("G:/My Drive/Chapter1_ET_Project/Data/landcover/LANDSAT_FOREST/landsat_permanent_forest_resampled.tif")
forest_mask_landsat[forest_mask_landsat == 0] <- NA
anom_stk_l <- mask(anom_stk_l, forest_mask_landsat)

library(rgdal)
sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")

########################################################################################
## diverging color bar for all plots 
quants <- quantile(c(values(anom_stk_m), values(anom_stk_l)), c(0.02, 0.98), na.rm=T)

#col5 <- colorRampPalette(c('red', 'gray96', 'blue'))  #create color ramp starting from blue to red
col5 <- colorRampPalette(c("#d73027", 'gray96', "#313695"))  #create color ramp starting from blue to red
color_levels=50 #the number of colors to use
max_absolute_value=max(abs(c(quants[[1]], quants[[2]]))) #what is the maximum absolute value of raster?
color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)

par(mfrow = c(2,4), mai = c(0, 0, 0.2, 0), omi=c(0,0.1,0.2,1))
plot(anom_stk_m[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "All", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Moderate", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Severe", cex.main = 2)
plot(sbr, add = T)
plot(anom_stk_m[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F, main = "Extreme", cex.main = 2)
plot(sbr, add = T)

plot(anom_stk_l[[1]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[2]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[3]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)
plot(anom_stk_l[[4]], col=col5(n=color_levels), breaks=color_sequence, zlim = c(quants[[1]], quants[[2]]), legend=FALSE, axes=F, box = F)
plot(sbr, add = T)

par(mfrow=c(1, 1), mai=c(0.5, 6, 0.5, 0),omi=c(0,0.8,0,0), new=FALSE)
plot(anom_stk_m[[1]], legend.only=TRUE, legend.shrink=1, legend.width=1, 
     zlim=c(quants[[1]], quants[[2]]),
     col=col5(n=color_levels), breaks=color_sequence,
     axis.args=list(at=pretty(quants[[1]]:quants[[2]]), font = 2,labels=pretty(quants[[1]]:quants[[2]])),
     legend.args=list(text='', side=4, font=2, line=2.3))







############################################################################################
## Make plots of the ETdp across gradients broken down by drought severity 
############################################################################################
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

# start with landsat get the .csv and make the plots individually 
dir <- "G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat"
l_e <- contPlots_1(paste0(dir, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
l_d <- contPlots_1(paste0(dir, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
l_h <- contPlots_1(paste0(dir, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "HAND (m)")


# get the MODIS data_long_sub dataframe initialized in modis_visualize.R and then run this 
m_e <- contPlots(data_long_sub, "elevation_bin", "DroughtSeverity", "Elevation (m)")
m_d <- contPlots(data_long_sub, "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
m_h <- contPlots(data_long_sub, "HAND_bin", "DroughtSeverity", "HAND (m)")


# Add the sen's slope of M, S, and E drought with a * for < 0.05 and ** < 0.01, and *** < 0.001
droughts <- c("Moderate", "Severe", "Extreme")

# do elevation first for MODIS and Landsat 
results <- l_e[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
LE <- l_e[[1]] + annotate("text", x = 1500, y = c(-0.61, -0.85, -1.07), label = labels)

results <- m_e[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
ME <- m_e[[1]] + annotate("text", x = 1250, y = c(-1, -1.2, -1.4), label = labels)


# Next do Hand 
results <- l_h[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
LH <- l_h[[1]] + annotate("text", x = 350, y = c(-0.33, -0.48, -0.63), label = labels)

results <- m_h[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
MH <- m_h[[1]] + annotate("text", x = 110, y = c(-0.88, -1.03, -1.2), label = labels)



# Finally do % Diff porous BA
results <- l_d[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
LD <- l_d[[1]] + annotate("text", x = 70, y = c(-0.3, -0.45, -0.60), label = labels)


results <- m_d[[2]]
labels <- c()
for(i in 1:length(droughts)){
  s <- results[DroughtSeverity == droughts[[i]], 2]
  p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
              ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                     ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
  labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
  
}
MD <- m_d[[1]] + annotate("text", x = 55, y = c(-0.9, -1.1, -1.3), label = labels)


library(cowplot)
plot_grid(ME, LE, 
         MH, LH, 
          MD, LD,
          labels = "AUTO", nrow = 3)


