########################################################################################################
########################################################################################################
## Function: Calculate the mean and +/- 1 sd of the Y variable(ET anomaly and residual) by the binned 
#### independent variable and grouping variable (drought severity) and save as a .csv for plotting later. 
#### agg_DF is for Landsat and agg_DF_m is for MODIS. They are the same except for file path I'm writing out to. 
## Inputs 
##### DF= dataframe with FC, elevation, Hand and ET anomalies averaged by all, extreme, severe, and moderate drought
##### X = independent variable that has been binned (elevation, HAND, or % diffuse porous)
##### Z = grouping variable (Drought severity) 
## Ouputs 
#### a .csv with columns of the mean ET anomaly and residual, +/- 1 sd for each, the X bin, and the drought severity
########################################################################################################
########################################################################################################
agg_DF <- function(DF, X, Z, out){
  agg <- DF[, .(meanAnom = mean(ET_Anomaly, na.rm=T), 
                meanRes = mean(ET_Residual, na.rm=T), 
                sdAnom = sd(ET_Anomaly, na.rm=T), 
                sdRes = sd(ET_Residual, na.rm=T), 
                count = .N), by=c(X, Z)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,] # only keep bins with at least 100 observations 
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  agg$res_plus_sd <- agg$meanRes + agg$sdRes
  agg$res_minus_sd <- agg$meanRes - agg$sdRes
  
  
  fwrite(agg, paste0(home, "/Visualizing/Landsat", out))
}


agg_DF_m <- function(DF, X, Z, out){
  agg <- DF[, .(meanAnom = mean(ET_Anomaly, na.rm=T), 
                meanRes = mean(ET_Residual, na.rm=T), 
                sdAnom = sd(ET_Anomaly, na.rm=T), 
                sdRes = sd(ET_Residual, na.rm=T), 
                count = .N), by=c(X, Z)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  agg$res_plus_sd <- agg$meanRes + agg$sdRes
  agg$res_minus_sd <- agg$meanRes - agg$sdRes
  
  
  fwrite(agg, paste0(home, "/Visualizing/MODIS", out))
}


########################################################################################################
########################################################################################################
## Function: Calculate the mean and +/- 1 sd of the Y variable by the binned X variable
#### agg_DF_XY is for Landsat and agg_DF_XY_m is for MODIS. They are the same except for file path I'm writing out to. 
## Inputs 
##### DF= dataframe with FC, elevation, Hand and ET anomalies averaged by all, extreme, severe, and moderate drought
##### X = independent variable that has been binned (elevation, HAND, or % diffuse porous)
##### Y = dependent variable to be summarized 
## Ouputs 
#### a .csv with columns of the mean Y variable, +/- 1 sd for the Y variable and the X bin 
##########################################################################################################
##########################################################################################################
agg_DF_XY <- function(DF, X, Y, out){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,] # don't keep groups with less than 100 observations
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  agg <- agg[order(agg[[1]]),]
  fwrite(agg, paste0(home, "/Visualizing/Landsat", out))
}

agg_DF_XY_m <- function(DF, X, Y, out){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  agg <- agg[order(agg[[1]]),]
  fwrite(agg, paste0(home, "/Visualizing/MODIS", out))
}

########################################################################################################
########################################################################################################
## Function: Calculate the sen's slope and  pvalue of the ET at drought peak across the binned gradient of interest 
#### (elevation, hand, or % diff) and broken down by drought severity. Also plot the ET at drought peak across the 
#### binned graadient broken down by drought severity 
## Inputs 
#### DF: the dataframe with the binned mean and +/- 1 sd ET at drought peak, binned variable of interest, and drought severity created from agg_XY
#### X: the variable of interest column name in the DF 
#### Z: the groupling variable of inteest column name in the DF 
#### xlab: the name for the x label in the plot 
## Ouputs 
#### a ggplot 
#### a dataframe of the sens slope and pvalue of the trend of ET anomalies arcoss the var of interest broke down by drought severity 
########################################################################################################################
contPlots_1 <- function(DF, X, Z, xlab){
  library(trend)
  agg <- fread(DF)
  agg$DroughtSeverity <- ordered(agg$DroughtSeverity, levels = c("None", "Moderate", "Severe", "Extreme", "All"))
  groups <- unique(agg[,..Z])
  results <- data.table(groups[,1], senslope_anom = as.numeric(rep(NA, nrow(groups))), pval_anom = as.numeric(rep(NA, nrow(groups))), 
                        senslope_res = as.numeric(rep(NA, nrow(groups))), pval_res = as.numeric(rep(NA, nrow(groups))))
  for(i in 1:nrow(results)){
    sub <- agg[agg[[Z]] == results[[1]][i],]
    sub <- sub[order(sub[[1]]),]
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
  
  ap <- ggplot(data = agg, aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(expression(bold(ET["dp"]))) +
    theme_classic()+
    theme(legend.title = element_blank()) + 
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12, color = 'black'),
          axis.title=element_text(size=12,face="bold"))
  
  results_list <- list(ap, results)
  return(results_list)
}

########################################################################################################
########################################################################################################
## Function: Make labels for each group of the sens slope and the significance. This is used for the plots produced by contPlots_1
## Inputs 
#### DF: The DF produced in the contPlots_1 function
## Ouputs: 
#### text that serves as labels on the plots 
########################################################################################################################
makeLabels <- function(DF){
  droughts <- c("Moderate", "Severe", "Extreme")
  results <- DF
  labels <- c()
  for(i in 1:length(droughts)){
    s <- results[DroughtSeverity == droughts[[i]], 2]
    p <- ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.001, "***", 
                ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.01, "**", 
                       ifelse(results[DroughtSeverity == droughts[[i]], 3] <= 0.05, "*", "")))
    labels <- c(labels, paste0(droughts[[i]], " Slope = ", round(s, 4), p))
    
  }
  
  return(labels)
}

##################################################################################################################
## Function to calculate the % of anomalies greater than 0 at each drought level 
#### and break down the percentage of anomalies greater than 0 across elevation, HAND, and pctDiff_R gradients 
## Input: 
#### DT: the overall DT with columns of the average ET anomaly by pixel at each level of drought severity 
#### home: filepath to save .csvs
#### sensor: Landsat or MODIS 
## Output: 
#### a csv file of the percent of anomalies greater than 0 
##################################################################################################################
pct_greater_0 <- function(DT, home, sensor){
  
  pctGR0 <- DT[, .(pg0=sum(ET_Anomaly >0)/.N),.(DroughtSeverity)]
  
  pctGR0_elev <- DT[, .(pg0=sum(ET_Anomaly >0)/.N),.(DroughtSeverity, elevation_bin)]
  pctGR0_hand <- DT[, .(pg0=sum(ET_Anomaly >0)/.N),.(DroughtSeverity, HAND_bin)] 
  pctGR0_pctDiffR <- DT[, .(pg0=sum(ET_Anomaly >0)/.N),.(DroughtSeverity, pctDiff_R_bin)]
  
  fwrite(pctGR0, paste0(home, "/Visualizing/", sensor, "/pct_greater0.csv"))
  fwrite(pctGR0_elev, paste0(home, "/Visualizing/", sensor, "/pct_greater0_elevation.csv"))
  fwrite(pctGR0_hand, paste0(home, "/Visualizing/", sensor, "/pct_greater0_HAND.csv"))
  fwrite(pctGR0_pctDiffR, paste0(home, "/Visualizing/", sensor, "/pct_greater0_pctDiffR.csv"))
  
}



########################################################################################################
########################################################################################################
## Function: Calculate the sen's slope and  pvalue of some variable across a gradient for MODIS and Landsat 
## Inputs 
#### modis_path: the path to the modis DF
#### landsat_path: the path to the landsat DF
#### VAR: the column name of the gradient the trend is calculated across 
## Ouputs 
#### A dataframe of the combined landsat and modis data 
#### The MODIS label that tells the sens slope and significance 
#### The Landsat label that tells the sens slope and signifiance
########################################################################################################################
prep <- function(modis_path, landsat_path, VAR){
  m_elev <- fread(modis_path)
  m_elev$Sensor <- rep("MODIS", nrow(m_elev))
  l_elev <- fread(landsat_path)
  l_elev$Sensor <- rep("Landsat", nrow(l_elev))
  
  combo_elev <- setDT(rbind(m_elev, l_elev))
  combo_elev[,(VAR):= as.numeric(as.character(combo_elev[[VAR]]))]
 
  m_elev[,(VAR) := as.numeric(as.character(m_elev[[VAR]]))]
  m_elev <- m_elev[order(m_elev[[VAR]]),]
  m_elev_sens <- sens.slope(m_elev$meanAnom)
  m_elev_label <- paste0("MODIS slope = ", round(m_elev_sens$estimates, 4), ifelse(m_elev_sens$p.value <= 0.001, "***",
                                                                                   ifelse(m_elev_sens$p.value <= 0.01, "**", 
                                                                                          ifelse(m_elev_sens$p.value <= 0.05, "*", ""))))
  
  
  l_elev[,(VAR) := as.numeric(as.character(l_elev[[VAR]]))]
  l_elev <- l_elev[order(l_elev[[VAR]]),]
  l_elev_sens <- sens.slope(l_elev$meanAnom)
  l_elev_label <- paste0("Landsat slope = ", round(l_elev_sens$estimates, 4), ifelse(l_elev_sens$p.value <= 0.001, "***",
                                                                                     ifelse(l_elev_sens$p.value <= 0.01, "**", 
                                                                                            ifelse(l_elev_sens$p.value <= 0.05, "*", ""))))
  
  results <- list(combo_elev, m_elev_label, l_elev_label)
  return(results)
}



