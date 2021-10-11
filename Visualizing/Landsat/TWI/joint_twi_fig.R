library(rasterVis)
library(raster)
library(pals)
library(RColorBrewer)
library(trend)
library(data.table)
library(viridis)
library(ggplot2)
library(rgdal)
library(cowplot)
library(ggpubr)
home <- "G:/My Drive/Chapter1_ET_Project"
source(paste0(home, "/Visualizing/MODIS/helpers.R"))

dir_l <- "G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/TWI"
dir_m <- "G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/TWI"

# Start with overall coupling 
### first do elevation
elev_corr_in <- prep(paste0(dir_m, "/elevation_allCoupling.csv"), paste0(dir_l, "/elevation_allCoupling.csv"), "elevation_bin")
elev_correlation <- ggplot(data =elev_corr_in[[1]], aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
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
  annotate("text", x = 725, y = c(-0.3, -0.45), label = c(elev_corr_in[[2]], elev_corr_in[[3]]))


### now do HAND 
hand_corr_in <- prep(paste0(dir_m, "/hand_allCoupling.csv"), paste0(dir_l, "/hand_allCoupling.csv"), "HAND_bin")
hand_correlation <- ggplot(data = hand_corr_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("TWI") + 
  ylab("R(SPI~ET)") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 15, y = c(0, -0.15), label = c(hand_corr_in[[2]], hand_corr_in[[3]]))

### now do percent diffuse porous 
diff_corr_in <- prep(paste0(dir_m, "/pctDiffR_allCoupling.csv"), paste0(dir_l, "/pctDiffR_allCoupling.csv"), "pctDiff_R_bin")
diff_correlation <- ggplot(data = diff_corr_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
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
  annotate("text", x = 25, y = c(-0.05, -0.12), label = c(diff_corr_in[[2]], diff_corr_in[[3]]))


# start the sen's slope plots 
elev_trend_in <- prep(paste0(dir_m, "/elevation_sensCoupling.csv"), paste0(dir_l, "/elevation_sensCoupling.csv"), "elevation_bin")
elev_trend <- ggplot(data = elev_trend_in[[1]], aes(x=elevation_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
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
  annotate("text", x = 600, y = c(-0.02, -0.028), label = c(elev_trend_in[[2]], elev_trend_in[[3]]))



### next do HAND 
hand_trend_in <- prep(paste0(dir_m, "/HAND_sensCoupling.csv"), paste0(dir_l, "/HAND_sensCoupling.csv"), "HAND_bin")
hand_trend <- ggplot(data = hand_trend_in[[1]], aes(x=HAND_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("TWI") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1))+ 
  annotate("text", x = 15, y = c(-0.022, -0.028), label = c(hand_trend_in[[2]], hand_trend_in[[3]]))


### next do % diffuse porous 
diff_trend_in <- prep(paste0(dir_m, "/pctDiffR_sensCoupling.csv"), paste0(dir_l, "/pctDiffR_sensCoupling.csv"), "pctDiff_R_bin")
diff_trend <- ggplot(data = diff_trend_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
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
  annotate("text", x = 25, y = c(-0.01, -0.017), label = c(diff_trend_in[[2]], diff_trend_in[[3]]))


## make a plot with just the legend -- no plot 
legend <- ggplot(data = diff_trend_in[[1]], aes(x=pctDiff_R_bin, y = meanAnom, group = Sensor)) + 
  geom_line(aes(color = Sensor),size = 1.5) + 
  #geom_point(aes(color = Sensor), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA, fill = Sensor), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  xlab("%Diffuse Porous") + 
  ylab("Trend") +
  theme_classic()+
  theme(legend.title = element_blank()) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 25, y = c(-0.023, -0.029), label = c(diff_trend_in[[2]], diff_trend_in[[3]]))
legend <- ggpubr::get_legend(legend)

#tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_7.tiff", units="in", width=8, height=8, res=800)
ggarrange(elev_correlation, elev_trend, 
          hand_correlation, hand_trend, 
          diff_correlation, diff_trend, 
          labels = "AUTO", nrow = 3, ncol = 2, 
          common.legend = TRUE, 
          legend.grob = legend,
          legend = "bottom")
#dev.off()



############################################################################################
## Make plots of the ETdp across gradients broken down by drought severity 
############################################################################################
# start with landsat get the .csv and make the plots individually 
l_e <- contPlots_1(paste0(dir_l, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
l_d <- contPlots_1(paste0(dir_l, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
l_h <- contPlots_1(paste0(dir_l, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "TWI")

# then do modis 
m_e <- contPlots_1(paste0(dir_m, "/elevation_DS.csv"), "elevation_bin", "DroughtSeverity", "Elevation (m)")
m_d <- contPlots_1(paste0(dir_m, "/pctDiffR_DS.csv"), "pctDiff_R_bin", "DroughtSeverity", "% Diffuse Porous")
m_h <- contPlots_1(paste0(dir_m, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "TWI")

# Add the sen's slope of M, S, and E drought with a * for < 0.05 and ** < 0.01, and *** < 0.001
# do elevation first for MODIS and Landsat 
le_labels <- makeLabels(l_e[[2]])
LE <- l_e[[1]] + annotate("text", x = 1500, y = c(-0.61, -0.85, -1.07), label = le_labels)
me_labels <- makeLabels(m_e[[2]])
ME <- m_e[[1]] + annotate("text", x = 1250, y = c(-1, -1.2, -1.4), label = me_labels)

# Next do Hand 
lh_labels <- makeLabels(l_h[[2]])
LH <- l_h[[1]] + annotate("text", x = 15, y = c(-0.33, -0.48, -0.63), label = lh_labels)
mh_labels <- makeLabels(m_h[[2]])
MH <- m_h[[1]] + annotate("text", x = 10, y = c(-0.88, -1.03, -1.2), label = mh_labels)


# Finally do % Diff porous BA
ld_labels <- makeLabels(l_d[[2]])
LD <- l_d[[1]] + annotate("text", x = 70, y = c(-0.3, -0.45, -0.60), label = ld_labels)
md_labels <- makeLabels(m_d[[2]])
MD <- m_d[[1]] + annotate("text", x = 55, y = c(-0.9, -1.1, -1.3), label = md_labels)

# Make a legend 
m_h_legend <- contPlots_legend(paste0(dir_m, "/HAND_DS.csv"), "HAND_bin", "DroughtSeverity", "HAND (m)")
m_h_legend <- ggpubr::get_legend(m_h_legend)

#tiff("G:/My Drive/Chapter1_ET_Project/Figures/figure_5.tiff", units="in", width=9, height=7, res=800)
ggarrange(ME, LE, 
          MH, LH, 
          MD, LD,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend.grob = m_h_legend,
          legend="bottom", 
          labels = "AUTO")
#dev.off()





#######################################################################33
## Look at variation in pct diff BA across hillslope bins 



contPlots_11 <- function(DF, X, Z, xlab){
  # this is different from regular contPlots because it additionally groups by elevation bin 
  library(trend)
  agg <- fread(DF)
  agg$DroughtSeverity <- ordered(agg$DroughtSeverity, levels = c("None", "Moderate", "Severe", "Extreme", "All"))
  agg$Grouping <- paste0(agg$DroughtSeverity, " ", agg$HAND_bin)
  groups <-  unique(agg$Grouping)
  results <- data.table(Group = groups, senslope_anom = as.numeric(rep(NA, length(groups))), pval_anom = as.numeric(rep(NA, length(groups))), 
                        senslope_res = as.numeric(rep(NA, length(groups))), pval_res = as.numeric(rep(NA, length(groups))))
  for(i in 1:nrow(results)){
    sub <- agg[agg$Grouping == results[[1]][i],]
    sub <- sub[order(sub[[1]]),]
    print(nrow(sub))
    if(nrow(sub) < 5){
      results$senslope_anom[i] <- NA
      results$pval_anom[i] <- NA
      results$senslope_res[i] <- NA
      results$pval_res[i] <- NA
    }else{
      anom_sen <- sens.slope(sub$meanAnom)
      results$senslope_anom[i] <- anom_sen$estimates
      results$pval_anom[i] <- anom_sen$p.value
      
      res_sen <- sens.slope(sub$meanRes)
      results$senslope_res[i] <- res_sen$estimates
      results$pval_res[i] <- res_sen$p.value
    }
    
  }
  # upslope
  ap1 <- ggplot(data = agg[agg$HAND_bin == "5-10",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 2) + 
    geom_point(size = 3) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(expression(bold(ET["dp"]))) +
    theme_classic()+
    theme(legend.title = element_blank()) + 
    #theme(legend.position = "bottom") +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12, color = 'black'),
          axis.title=element_text(size=12,face="bold")) 
  
  # mid slope
  ap2 <- ggplot(data = agg[agg$HAND_bin == "10-15",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 2) + 
    geom_point(size = 3) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(expression(bold(ET["dp"]))) +
    theme_classic()+
    theme(legend.title = element_blank()) + 
    #theme(legend.position = "bottom") +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12, color = 'black'),
          axis.title=element_text(size=12,face="bold")) 
  
  # downslope
  ap3 <- ggplot(data = agg[agg$HAND_bin == "15-30",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 2) + 
    geom_point(size = 3) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(expression(bold(ET["dp"]))) +
    theme_classic()+
    theme(legend.title = element_blank()) + 
    #theme(legend.position = "bottom") +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12, color = 'black'),
          axis.title=element_text(size=12,face="bold")) 
  
  
  results_list <- list(ap1, ap2, ap3, results)
  return(results_list)
}

modis_ETdp_diff_plots <- contPlots_11("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/TWI_bin/pctDiffR_DS1.csv", 
                                      "pctDiff_R_bin", 
                                      "DroughtSeverity",
                                      "% Diffuse Porous")

modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("\\d+", "", modis_ETdp_diff_plots[[4]]$Group)
modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("-", "", modis_ETdp_diff_plots[[4]]$DroughtSeverity)
modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub(" ", "", modis_ETdp_diff_plots[[4]]$DroughtSeverity)

m_upslope <- makeLabels(modis_ETdp_diff_plots[[4]][grep("5-10", Group),])
m_midslope <- makeLabels(modis_ETdp_diff_plots[[4]][grep("10-15", Group),])
m_downslope <- makeLabels(modis_ETdp_diff_plots[[4]][grep("15-30", Group),])
M1 <- modis_ETdp_diff_plots[[1]] + annotate("text", x = 35, y=c(-0.25, -0.45, -0.65), label=m_upslope) + ggtitle('Landsat: Upslope')
M2 <- modis_ETdp_diff_plots[[2]] + annotate("text", x = 35, y=c(-0.25, -0.45, -0.65), label=m_midslope) + ggtitle('Landsat: Midslope')
M3 <- modis_ETdp_diff_plots[[3]] + annotate("text", x = 35, y=c(-0.25, -0.45, -0.65), label=m_downslope)+ ggtitle('Landsat: Downslope')

ggarrange(M1, M2, M3,
          nrow = 1, ncol = 3,
          common.legend = TRUE,
          legend="bottom", 
          labels = "AUTO")






## Now make the plots for R(SPI-ET) and dR(SPI-ET) 

## function to take the 
sens_labels <- function(DF, VAR){
  hand_groups <- unique(DF$hand)
  labels <- c()
  for(i in 1:length(hand_groups)){
    sub <- DF[hand == hand_groups[i],]
    sub <- sub[order(get(VAR)),]
    ss <- sens.slope(sub$meanAnom)
    s <- ss$estimates[[1]]
    p <- ifelse(ss$p.value[[1]] <= 0.001, "***", 
                ifelse(ss$p.value[[1]] <= 0.01, "**", 
                       ifelse(ss$p.value[[1]] <= 0.05, "*", "")))
    labels <- c(labels, paste0(hand_groups[[i]], ": Slope = ", round(s, 4), p))
    
  }
  return(labels)
}


landsat_R <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/TWI_bin/pctDiffR_allCoupling.csv")
landsat_R$hand <- ifelse(landsat_R$HAND_bin == "5-10", "Upslope", 
                              ifelse(landsat_R$HAND_bin == "10-15", 'Midslope', 
                                     ifelse(landsat_R$HAND_bin == "15-30", 'Downslope', NA)))
L_R_labs <- sens_labels(landsat_R, "pctDiff_R_bin")
L_R <- ggplot(data = landsat_R, aes(x=pctDiff_R_bin, y = meanAnom, color = hand, fill = hand)) + 
  geom_line(aes(color = hand),size = 1.5) + 
  geom_point(aes(color = hand), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("% Diffuse Porous") + 
  ylab("R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 40, y = c(-0.25, -0.3, -0.35), label = L_R_labs)+ 
  ggtitle("Landsat")



landsat_dR <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/TWI_bin/pctDiffR_sensCoupling.csv")
landsat_dR$hand <- ifelse(landsat_dR$HAND_bin == "5-10", "Upslope", 
                         ifelse(landsat_dR$HAND_bin == "10-15", 'Midslope', 
                                ifelse(landsat_dR$HAND_bin == "15-30", 'Downslope', NA)))
L_dR_labs <- sens_labels(landsat_dR, "pctDiff_R_bin")
L_dR <- ggplot(data = landsat_dR, aes(x=pctDiff_R_bin, y = meanAnom, color = hand, fill = hand)) + 
  geom_line(aes(color = hand),size = 1.5) + 
  geom_point(aes(color = hand), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("% Diffuse Porous") + 
  ylab("Trend in R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 40, y = c(-0.005, -0.0075, -0.01), label = L_dR_labs)

ggarrange(L_R, L_dR,
          nrow = 1, ncol = 2, 
          common.legend = TRUE, 
          legend = 'bottom', 
          labels = "AUTO")
