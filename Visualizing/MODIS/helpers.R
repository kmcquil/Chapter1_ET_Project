#####################################################################################################
## Function to make the color bar diverge at 0 
#####################################################################################################
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)    
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=101)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(100+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(100)[zlim[-length(zlim)]]
  p
}


########################################################################################################
## Function to plot relationships between binned conntinuous variables, grouping variable, and the ET anom and res
## Inputs 
##### DF= dataframe (this way i can subset before putting in )
##### X = independent variable 
##### Z = grouping variable 
## Ouputs 
##### plot of the anomalies 
##### plot of the residuals 
########################################################################################################

contPlots <- function(DF, X, Z, xlab){
  library(trend)
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
  
  groups <- unique(agg[,..Z])
  results <- data.table(groups[,1], senslope_anom = as.numeric(rep(NA, nrow(groups))), pval_anom = as.numeric(rep(NA, nrow(groups))), 
                        senslope_res = as.numeric(rep(NA, nrow(groups))), pval_res = as.numeric(rep(NA, nrow(groups))))
  for(i in 1:nrow(results)){
    sub <- agg[agg[[Z]] == results[[1]][i],]
    sub <- sub[order(sub[[1]]),]
    anom_sen <- sens.slope(sub$meanAnom)
    results$senslope_anom[i] <- anom_sen$estimates
    results$pval_anom[i] <- anom_sen$p.value
    
    res_sen <- sens.slope(sub$meanRes)
    results$senslope_res[i] <- res_sen$estimates
    results$pval_res[i] <- res_sen$p.value
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
  

  rp <- ggplot(data = agg, aes(x = !!ensym(X), y = meanRes, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    geom_ribbon(aes(ymin = res_minus_sd, ymax = res_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab("ET Residual") +
    theme_bw()+
    theme(legend.title = element_blank())
  

  #library(cowplot)
  #results_list <- list(plot_grid(ap, rp, labels = "AUTO"), results)
  results_list <- list(ap, results)
  return(results_list)
  
  
}




# DF <- data
# X <- "elevation_bin"
# Z <- "DroughtSeverity"
# Y <- "senSlope"


contPlots_single <- function(DF, X,Y, xlab, ylab){
  library(trend)
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom

  agg <- agg[order(agg[[1]]),]
  anom_sen <- sens.slope(agg$meanAnom)
  results <- data.table(VAR = X, 
                        senslope = anom_sen$estimates, 
                        pval = anom_sen$p.value)
  
  ap <- ggplot(data = agg, aes(x = !!ensym(X), y = meanAnom, group = 1)) + 
    geom_line(size = 1.5, color = viridis(1)) + 
    geom_point(size = 4, color = viridis(1), fill = viridis(1)) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), fill = viridis(1), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(ylab) +
    theme_bw()+
    theme(legend.title = element_blank()) + 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16,face="bold"))
  
  
  results_list <- list(ap, results)
  return(results_list)
  
}

# DF <- data
# X <- "pctDiff_R_bin"
# Y <- "all"
# Z <- "Elev_class"
# xlab = "Binned % Diffuse Porous"
# ylab = "R (ET ~ SPI)"


contPlots_rolling <- function(DF, X, Y, Z, xlab, ylab){
  library(trend)
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T),
                sdAnom = sd(get(Y), na.rm=T),
                count = .N), by=c(X, Z)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  
  groups <- unique(agg[,..Z])
  results <- data.table(groups[,1], senslope_anom = as.numeric(rep(NA, nrow(groups))), pval_anom = as.numeric(rep(NA, nrow(groups))))
  for(i in 1:nrow(results)){
    sub <- agg[agg[[Z]] == results[[1]][i],]
    sub <- sub[order(sub[[1]]),]
    anom_sen <- sens.slope(sub$meanAnom)
    results$senslope_anom[i] <- anom_sen$estimates
    results$pval_anom[i] <- anom_sen$p.value
  }
  
  ap <- ggplot(data = agg, aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 1.5) + 
    geom_point(size = 4) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(ylab) +
    theme_bw()+
    theme(legend.title = element_blank()) 
  

  results_list <- list(ap, results)
  return(results_list)
  
  
}


agg_DF <- function(DF, X, Z, out){
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
  
  
  fwrite(agg, paste0(home, "/Visualizing/Landsat", out))
}



agg_DF_XY <- function(DF, X, Y, out){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T), 
                sdAnom = sd(get(Y), na.rm=T), 
                count = .N), by=c(X)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  
  agg <- agg[order(agg[[1]]),]
  fwrite(agg, paste0(home, "/Visualizing/Landsat", out))
}


agg_DF_XYZ <- function(DF, X, Y, Z, out){
  agg <- DF[, .(meanAnom = mean(get(Y), na.rm=T),
                sdAnom = sd(get(Y), na.rm=T),
                count = .N), by=c(X, Z)]
  agg <- agg[complete.cases(agg),]
  agg <- agg[!count < 100,]
  
  agg$anom_plus_sd <- agg$meanAnom + agg$sdAnom
  agg$anom_minus_sd <- agg$meanAnom - agg$sdAnom
  fwrite(agg, paste0(home, "/Visualizing/Landsat", out))
}






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
  
  
  
  
  rp <- ggplot(data = agg, aes(x = !!ensym(X), y = meanRes, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 1.5) + 
    geom_point(size = 4) + 
    geom_ribbon(aes(ymin = res_minus_sd, ymax = res_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab("ET Residual") +
    theme_bw()+
    theme(legend.title = element_blank())
  
  
  #library(cowplot)
  #results_list <- list(plot_grid(ap, rp, labels = "AUTO"), results)
  results_list <- list(ap, results)
  return(results_list)
  
  
}





contPlots_single_1 <- function(DF, X,Y, xlab, ylab){
  library(trend)
  agg <- fread(DF)
  anom_sen <- sens.slope(agg$meanAnom)
  results <- data.table(VAR = X, 
                        senslope = anom_sen$estimates, 
                        pval = anom_sen$p.value)
  
  ap <- ggplot(data = agg, aes(x = !!ensym(X), y = meanAnom, group = 1)) + 
    geom_line(size = 1.5, color = viridis(1)) + 
    geom_point(size = 4, color = viridis(1), fill = viridis(1)) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), fill = viridis(1), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(ylab) +
    theme_bw()+
    theme(legend.title = element_blank()) + 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16,face="bold"))
  
  
  results_list <- list(ap, results)
  return(results_list)
  
}



contPlots_rolling_1 <- function(DF, X, Y, Z, xlab, ylab){
  library(trend)
  agg <- fread(DF)
  
  groups <- unique(agg[,..Z])
  results <- data.table(groups[,1], senslope_anom = as.numeric(rep(NA, nrow(groups))), pval_anom = as.numeric(rep(NA, nrow(groups))))
  for(i in 1:nrow(results)){
    sub <- agg[agg[[Z]] == results[[1]][i],]
    sub <- sub[order(sub[[1]]),]
    anom_sen <- sens.slope(sub$meanAnom)
    results$senslope_anom[i] <- anom_sen$estimates
    results$pval_anom[i] <- anom_sen$p.value
  }
  
  ap <- ggplot(data = agg, aes(x = !!ensym(X), y = meanAnom, color = !!ensym(Z), fill = !!ensym(Z))) + 
    geom_line(size = 1.5) + 
    geom_point(size = 4) + 
    geom_ribbon(aes(ymin = anom_minus_sd, ymax = anom_plus_sd), alpha = 0.20, color = NA) + 
    scale_color_manual(values = viridis(4))+
    scale_fill_manual(values = viridis(4))+
    xlab(xlab) + 
    ylab(ylab) +
    theme_bw()+
    theme(legend.title = element_blank()) 
  
  
  results_list <- list(ap, results)
  return(results_list)
  
  
}









##################################################################################################################
## Function to calculate the % of anomalies greater than 0 at each drought level 
## and break down the percentage of anomalies greater than 0 across elevation, HAND, and pctDiff_R gradients 
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







