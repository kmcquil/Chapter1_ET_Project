##########################################################################################
## now make the plots of HAND and % diffuse porous broken down over elevation (high medium low )
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


##########################################################################################
# Start looking at ETdp 
contPlots_11 <- function(DF, X, Z, xlab){
  # this is different from regular contPlots because it additionally groups by elevation bin 
  library(trend)
  agg <- fread(DF)
  agg$DroughtSeverity <- ordered(agg$DroughtSeverity, levels = c("None", "Moderate", "Severe", "Extreme", "All"))
  agg$Grouping <- paste0(agg$DroughtSeverity, " ", agg$elevation_bin)
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
  
  ap1 <- ggplot(data = agg[agg$elevation_bin == "200-800",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
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
  ap2 <- ggplot(data = agg[agg$elevation_bin == "800-1400",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
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
  ap3 <- ggplot(data = agg[agg$elevation_bin == "1400-2000",], aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
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

# MODIS 
# create the plots 
modis_ETdp_diff_plots <- contPlots_11("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/pctDiffR_DS1.csv", 
            "pctDiff_R_bin", 
            "DroughtSeverity",
            "% Diffuse Porous")

# make the labels 
modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("\\d+", "", modis_ETdp_diff_plots[[4]]$Group)
modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("-", "", modis_ETdp_diff_plots[[4]]$DroughtSeverity)
modis_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub(" ", "", modis_ETdp_diff_plots[[4]]$DroughtSeverity)

m_low <- makeLabels(modis_ETdp_diff_plots[[4]][grep("200-800", Group),])
m_med <- makeLabels(modis_ETdp_diff_plots[[4]][grep("800-1400", Group),])
m_high <- makeLabels(modis_ETdp_diff_plots[[4]][grep("1400-2000", Group),])
M1 <- modis_ETdp_diff_plots[[1]] + annotate("text", x = 25, y=c(-0.85, -1.05, -1.2), label=m_low) + ggtitle('MODIS: Low Elevation')
M2 <- modis_ETdp_diff_plots[[2]] + annotate("text", x = 25, y=c(-0.85, -1.05, -1.2), label=m_med) + ggtitle('Mid Elevation')
M3 <- modis_ETdp_diff_plots[[3]] + annotate("text", x = 35, y=c(-0.75, -0.95, -1.15), label=m_high)+ ggtitle('High Elevation')

# Now do Landsat
landsat_ETdp_diff_plots <- contPlots_11("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/pctDiffR_DS1.csv", 
                                      "pctDiff_R_bin", 
                                      "DroughtSeverity",
                                      "% Diffuse Porous")
landsat_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("\\d+", "", landsat_ETdp_diff_plots[[4]]$Group)
landsat_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub("-", "", landsat_ETdp_diff_plots[[4]]$DroughtSeverity)
landsat_ETdp_diff_plots[[4]]$DroughtSeverity <- gsub(" ", "", landsat_ETdp_diff_plots[[4]]$DroughtSeverity)

l_low <- makeLabels(landsat_ETdp_diff_plots[[4]][grep("200-800", Group),])
l_med <- makeLabels(landsat_ETdp_diff_plots[[4]][grep("800-1400", Group),])
l_high <- makeLabels(landsat_ETdp_diff_plots[[4]][grep("1400-2000", Group),])
L1 <- landsat_ETdp_diff_plots[[1]] + annotate("text", x = 30, y=c(-0.6, -0.8, -1), label=l_low) + ggtitle('Landsat: Low Elevation')
L2 <- landsat_ETdp_diff_plots[[2]] + annotate("text", x = 25, y=c(-0.30, -0.5, -0.7), label=l_med) + ggtitle('Mid Elevation')
L3 <- landsat_ETdp_diff_plots[[3]] + annotate("text", x = 35, y=c(-0.25, -0.45, -0.65), label=l_high) + ggtitle('High Elevation')

ggarrange(M1, L1, 
          M2, L2, 
          M3, L3,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend.grob = m_h_legend,
          legend="bottom", 
          labels = "AUTO")



################################################3
## Now do for HAND 
modis_ETdp_hand_plots <- contPlots_11("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/HAND_DS1.csv", 
                                      "HAND_bin", 
                                      "DroughtSeverity",
                                      "HAND (m)")
# make the labels 
modis_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub("\\d+", "", modis_ETdp_hand_plots[[4]]$Group)
modis_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub("-", "", modis_ETdp_hand_plots[[4]]$DroughtSeverity)
modis_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub(" ", "", modis_ETdp_hand_plots[[4]]$DroughtSeverity)

m_low <- makeLabels(modis_ETdp_hand_plots[[4]][grep("200-800", Group),])
m_med <- makeLabels(modis_ETdp_hand_plots[[4]][grep("800-1400", Group),])
m_high <- makeLabels(modis_ETdp_hand_plots[[4]][grep("1400-2000", Group),])
M1 <- modis_ETdp_hand_plots[[1]] + annotate("text", x = 40, y=c(-0.85, -1.05, -1.2), label=m_low) + ggtitle('MODIS: Low Elevation')
M2 <- modis_ETdp_hand_plots[[2]] + annotate("text", x = 60, y=c(-0.85, -1.05, -1.2), label=m_med) + ggtitle('Mid Elevation')
M3 <- modis_ETdp_hand_plots[[3]] + annotate("text", x = 90, y=c(-0.75, -0.95, -1.15), label=m_high)+ ggtitle('High Elevation')


landsat_ETdp_hand_plots <- contPlots_11("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/HAND_DS1.csv", 
                                        "HAND_bin", 
                                        "DroughtSeverity",
                                        "HAND (m)")

landsat_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub("\\d+", "", landsat_ETdp_hand_plots[[4]]$Group)
landsat_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub("-", "", landsat_ETdp_hand_plots[[4]]$DroughtSeverity)
landsat_ETdp_hand_plots[[4]]$DroughtSeverity <- gsub(" ", "", landsat_ETdp_hand_plots[[4]]$DroughtSeverity)

l_low <- makeLabels(landsat_ETdp_hand_plots[[4]][grep("200-800", Group),])
l_med <- makeLabels(landsat_ETdp_hand_plots[[4]][grep("800-1400", Group),])
l_high <- makeLabels(landsat_ETdp_hand_plots[[4]][grep("1400-2000", Group),])
L1 <- landsat_ETdp_hand_plots[[1]] + annotate("text", x = 150, y=c(-0.6, -0.8, -1), label=l_low) + ggtitle('Landsat: Low Elevation')
L2 <- landsat_ETdp_hand_plots[[2]] + annotate("text", x = 150, y=c(-0.30, -0.5, -0.7), label=l_med) + ggtitle('Mid Elevation')
L3 <- landsat_ETdp_hand_plots[[3]] + annotate("text", x = 155, y=c(-0.25, -0.45, -0.65), label=l_high) + ggtitle('High Elevation')

ggarrange(M1, L1, 
          M2, L2, 
          M3, L3,
          nrow = 3, ncol = 2,
          common.legend = TRUE, 
          legend.grob = m_h_legend,
          legend="bottom", 
          labels = "AUTO")




#################################################################################################
## Now make the plots for R(SPI-ET) and dR(SPI-ET) 

## function to take the 
sens_labels <- function(DF, VAR){
  elevation_groups <- unique(DF$elevation)
  labels <- c()
  for(i in 1:length(elevation_groups)){
    sub <- DF[elevation == elevation_groups[i],]
    sub <- sub[order(get(VAR)),]
    ss <- sens.slope(sub$meanAnom)
    s <- ss$estimates[[1]]
    p <- ifelse(ss$p.value[[1]] <= 0.001, "***", 
                ifelse(ss$p.value[[1]] <= 0.01, "**", 
                       ifelse(ss$p.value[[1]] <= 0.05, "*", "")))
    labels <- c(labels, paste0(elevation_groups[[i]], ": Slope = ", round(s, 4), p))
  }
  return(labels)
}


modis_R <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/pctDiffR_allCoupling1.csv")
modis_R$elevation <- ifelse(modis_R$elevation_bin == "200-800", "Low", 
                            ifelse(modis_R$elevation_bin == "800-1400", 'Mid', 
                                   ifelse(modis_R$elevation_bin == "1400-2000", 'High', NA)))
M_R_labs <- sens_labels(modis_R, "pctDiff_R_bin")
M_R <- ggplot(data = modis_R, aes(x=pctDiff_R_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
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
  annotate("text", x = 40, y = c(0, 0.05, 0.1), label = M_R_labs) + 
  ggtitle("MODIS")


modis_dR <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/pctDiffR_sensCoupling1.csv")
modis_dR$elevation <- ifelse(modis_dR$elevation_bin == "200-800", "Low", 
                            ifelse(modis_dR$elevation_bin == "800-1400", 'Mid', 
                                   ifelse(modis_dR$elevation_bin == "1400-2000", 'High', NA)))
M_dR_labs <- sens_labels(modis_dR, "pctDiff_R_bin")
M_dR <- ggplot(data = modis_dR, aes(x=pctDiff_R_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
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
  annotate("text", x = 40, y = c(-0.02, -0.025, -0.03), label = M_dR_labs)+ 
  ggtitle("MODIS")



# Now Landsat 
landsat_R <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/pctDiffR_allCoupling.csv")
landsat_R$elevation <- ifelse(landsat_R$elevation_bin == "200-800", "Low", 
                            ifelse(landsat_R$elevation_bin == "800-1400", 'Mid', 
                                   ifelse(landsat_R$elevation_bin == "1400-2000", 'High', NA)))
L_R_labs <- sens_labels(landsat_R, "pctDiff_R_bin")
L_R <- ggplot(data = landsat_R, aes(x=pctDiff_R_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
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

landsat_dR <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/pctDiffR_sensCoupling.csv")
landsat_dR$elevation <- ifelse(landsat_dR$elevation_bin == "200-800", "Low", 
                             ifelse(landsat_dR$elevation_bin == "800-1400", 'Mid', 
                                    ifelse(landsat_dR$elevation_bin == "1400-2000", 'High', NA)))
L_dR_labs <- sens_labels(landsat_dR, "pctDiff_R_bin")
L_dR <- ggplot(data = landsat_dR, aes(x=pctDiff_R_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
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
  annotate("text", x = 40, y = c(-0.01, -0.015, -0.02), label = L_dR_labs)+ 
  ggtitle("Landsat")

ggarrange(M_R, L_R, 
          M_dR, L_dR, 
          nrow = 2, ncol = 2, 
          common.legend = TRUE, 
          legend = 'bottom', 
          labels = "AUTO")
















############################################33
# Now do HAND
modis_R <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/HAND_allCoupling1.csv")
modis_R$elevation <- ifelse(modis_R$elevation_bin == "200-800", "Low", 
                            ifelse(modis_R$elevation_bin == "800-1400", 'Mid', 
                                   ifelse(modis_R$elevation_bin == "1400-2000", 'High', NA)))
M_R_labs <- sens_labels(modis_R, "HAND_bin")
M_R <- ggplot(data = modis_R, aes(x=HAND_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("HAND (m)") + 
  ylab("R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 50, y = c(0.2, 0.18, 0.15), label = M_R_labs) + 
  ggtitle("MODIS")


modis_dR <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/MODIS/trial/HAND_sensCoupling1.csv")
modis_dR$elevation <- ifelse(modis_dR$elevation_bin == "200-800", "Low", 
                             ifelse(modis_dR$elevation_bin == "800-1400", 'Mid', 
                                    ifelse(modis_dR$elevation_bin == "1400-2000", 'High', NA)))
M_dR_labs <- sens_labels(modis_dR, "HAND_bin")
M_dR <- ggplot(data = modis_dR, aes(x=HAND_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("HAND (m)") + 
  ylab("Trend in R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 50, y = c(-0.02, -0.025, -0.03), label = M_dR_labs)+ 
  ggtitle("MODIS")


# Now Landsat 
landsat_R <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/HAND_allCoupling.csv")
landsat_R$elevation <- ifelse(landsat_R$elevation_bin == "200-800", "Low", 
                              ifelse(landsat_R$elevation_bin == "800-1400", 'Mid', 
                                     ifelse(landsat_R$elevation_bin == "1400-2000", 'High', NA)))
L_R_labs <- sens_labels(landsat_R, "HAND_bin")
L_R <- ggplot(data = landsat_R, aes(x=HAND_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("HAND (m)") + 
  ylab("R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 150, y = c(-0.25, -0.3, -0.35), label = L_R_labs)+ 
  ggtitle("Landsat")

landsat_dR <- fread("G:/My Drive/Chapter1_ET_Project/Visualizing/Landsat/trial/HAND_sensCoupling.csv")
landsat_dR$elevation <- ifelse(landsat_dR$elevation_bin == "200-800", "Low", 
                               ifelse(landsat_dR$elevation_bin == "800-1400", 'Mid', 
                                      ifelse(landsat_dR$elevation_bin == "1400-2000", 'High', NA)))
L_dR_labs <- sens_labels(landsat_dR, "HAND_bin")
L_dR <- ggplot(data = landsat_dR, aes(x=HAND_bin, y = meanAnom, color = elevation, fill = elevation)) + 
  geom_line(aes(color = elevation),size = 1.5) + 
  geom_point(aes(color = elevation), size = 4) + 
  geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd, color = NA), alpha = 0.2) + 
  scale_color_manual(values = viridis(4))+
  scale_fill_manual(values = viridis(4))+
  theme_classic() + 
  xlab("HAND (m)") + 
  ylab("Trend in R(SPI~ET)") +
  theme(legend.title = element_blank()) + 
  theme(axis.text=element_text(size=12, color = 'black'),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1)) + 
  annotate("text", x = 140, y = c(-0.01, -0.015, -0.02), label = L_dR_labs)+ 
  ggtitle("Landsat")

ggarrange(M_R, L_R, 
          M_dR, L_dR, 
          nrow = 2, ncol = 2, 
          common.legend = TRUE, 
          legend = 'bottom', 
          labels = "AUTO")










####################################################################################3
## Looking at TWI 
# MODIS 
home <- "G:/My Drive/Chapter1_ET_Project/Data/Topography/"
twi <- raster(paste0(home, "usgsNED_elevation/MODIS_TWI/TWI_SBR.tif"))
elevation <- raster(paste0(home, "usgsNED_elevation/elevation500m.tif"))
hand <- raster(paste0(home, "HeightAboveNearestDrainage/HAND_modis.tif"))
diffBA <- raster("G:/My Drive/Chapter1_ET_Project/Data/forest_composition/riley_modis_diffuse_percent.tif")


# Landsat 
# grab the TWI and crop and resample it to match elevation, hand, and % diff BA 
landsat_elevation <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/elevation30m.tif")
landsat_hand <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/HeightAboveNearestDrainage/hand_landsat.tif")
landsat_diffBA <- raster("G:/My Drive/Chapter1_ET_Project/Data/forest_composition/riley_landsat_diffuse_percent.tif")
sbr <- readOGR("G:/My Drive/Chapter1_ET_Project/Data/NA_CEC_Eco_Level3/blue_ridge.shp")
landsat_twi <- raster("G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/Landsat_TWI/TWI.tif")
landsat_twi <- mask(crop(landsat_twi, sbr), sbr)

# resample the landsat TWI to MODIS resolution 
modis_twi <- resample(landsat_twi, diffBA)
writeRaster(modis_twi, "G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/MODIS_TWI/TWI_resampled.tif",
            format = "GTiff", overwrite = T)
preds <- data.table(twi=values(modis_twi), 
                    elevation=values(elevation), 
                    hand=values(hand), 
                    diffBA=values(diffBA))
modis_correlations <- cor(preds[complete.cases(preds),])


# resampe landat twi to match landsat data 
landsat_twi <- resample(landsat_twi, landsat_diffBA)
writeRaster(landsat_twi, "G:/My Drive/Chapter1_ET_Project/Data/Topography/usgsNED_elevation/Landsat_TWI/TWI_resampled.tif", 
             format="GTiff", overwrite = T)


landsat_preds <- data.table(twi=values(landsat_twi), 
                    elevation=values(landsat_elevation), 
                    hand=values(landsat_hand), 
                    diffBA=values(landsat_diffBA))
landsat_correlations <- cor(landsat_preds[complete.cases(landsat_preds),])


