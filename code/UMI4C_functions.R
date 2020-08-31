## Functions for UMI4C analysis

smoothedTrend <- function(p4c_obj,
                          min_win_mols=NULL,
                          min_res=1,
                          sd=2) {
  if(is.null(min_win_mols)) min_win_mols <- getMinWinMols(sum(p4c_obj$dgram[,2]))
  
  ### Start code from umi4cpackage ------------
  dgram <- p4c_obj$dgram
  dgram[is.na(dgram)] <- 0
  bait_x <- p4c_obj$bait$start
  vec1 <- rep(NA, nrow(dgram))
  scale <- rep(NA, nrow(dgram))
  coord <- rep(NA, nrow(dgram))
  base_coord <- dgram[, 1]
  mean_coord <- numeric()
  base_scales <- p4c_obj$dgram_params$dgram_scales
  for (i in (2 + min_res):ncol(dgram)) {
    cur_scale <- base_scales[i - 2]
    f <- is.na(vec1) & (dgram[, i] * cur_scale * 2 > min_win_mols)
    vec1[f] <- dgram[f, i]
    mean_coords <- umi4cPackage:::.p4cWinGeoMeanCoordinate(base_coord, 
                                            cur_scale * 2, bait_x)
    coord[f] <- mean_coords[f]
    scale[f] <- cur_scale
  }
  vec1[is.na(vec1)] <- NA
  coord[is.na(coord)] <- base_coord[is.na(coord)]
  p4c_obj$smoothedTrend <- list(start = sort(coord), trend = vec1, 
                                scale = scale)
  ## End code from umi4cpackage
  
  ## Convert to df for ggplot2
  trends.df <- data.frame(start=p4c_obj$smoothedTrend$start,
                          trend=p4c_obj$smoothedTrend$trend,
                          scale=p4c_obj$smoothedTrend$scale)
  dev <- sd * sqrt(trends.df$trend/(trends.df$scale*2))
  trends.df$devP <- trends.df$trend + dev
  trends.df$devM <- trends.df$trend - dev
  
  ## Assign group to allow plotting lines separately
  trends.df$group <- "upstream"
  trends.df$group[trends.df$start > p4c_obj$bait$start] <- "downstream"
  return(trends.df)
}

getMinWinMols <- function(sumMols) 50*ceiling(sumMols/2000) # Min mols for smoothing

getTrend <- function(treat,
                     ctrl,
                     min_win_mols=NULL,
                     sd=2) {
  if(is.null(min_win_mols)) min_win_mols <- getMinWinMols(sum(ctrl$dgram[,2]))
  
  treat <- umi4cPackage:::.p4cNormDgram.p4cProfile(treat,
                                                   ctrl)
  comp <- umi4cPackage:::.p4cSmoothedTrendComp.p4cProfile(treat,
                                                          ctrl,
                                                          min_win_cov=min_win_mols)
  trends.df <- data.frame(comp$trend_mat)
  colnames(trends.df)[2:3] <- c("treat", "ctrl")
  trends.df <- reshape2::melt(trends.df,
                              id.vars=c(1,4),
                              value.vars=2:3,
                              variable.name="sample",
                              value.name="trend")
  
  trends.df <- trends.df[!is.na(trends.df$scale),]
  dev <- sd * sqrt(trends.df$trend/(trends.df$scale*2))
  trends.df$devP <- trends.df$trend + dev
  trends.df$devM <- trends.df$trend - dev
  
  ## Assign group to allow plotting lines separately
  trends.df$group <- "upstream"
  trends.df$group[trends.df$start > ctrl$bait$start] <- "downstream"
  
  return(trends.df)
}

plotUMI4C <- function(res,
                      diff,
                      xlim,
                      ymax
                      ) {
  
  
  umi <-
    ggplot(res$norm_trend,
           aes(start, trend)) +
    geom_ribbon(aes(ymin=devM, ymax=devP, group=interaction(group, sample)),
                color=NA, fill="light grey") +
    geom_line(aes(color=sample, group=interaction(group, sample)),
              lwd=0.7) +
    scale_color_manual(values=c(ctrl="#1f78b4", treat="#d95f02")) +
    scale_y_continuous(name="# UMI-4C contacts",
                       limits=c(0, NA),
                       breaks=scales::pretty_breaks()) +
    theme(legend.position="none")
}
