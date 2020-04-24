binRegions <- function(regions,
                       scope,
                       bin) {
  regions.ext <- resize(regions, width=scope*2, fix="center")
  regions.ext$center <- start(regions.ext) + scope
  
  # Bin regions
  regions.bin <- tile(regions.ext, width=bin)
  n <- unique(sapply(regions.bin, length))
  regions.bin <- unlist(regions.bin)
  regions.bin$GeneID <- rep(regions.ext$GeneID, each=n) ## add identifier
  regions.bin$center <- rep(regions.ext$center, each=n) ## add center
  regions.bin$pos <- start(regions.bin) - regions.bin$center ## add pos relative to center
  return(regions.bin)
}

getTSSenrichment <- function(bam_files,
                             saf,
                             suffix=".offset.bam") {
  counts <- Rsubread::featureCounts(files,
                                    annot.ext=saf,
                                    allowMultiOverlap=TRUE,
                                    nthreads=8)
  
  names <- pipelineNGS::getNameFromPath(files, suffix=suffix)
  colnames(counts$counts) <- names
  colnames(counts$stat) <- c("Status", names)
  
  anno <- cbind(saf, counts$counts)
  anno.m <- reshape2::melt(anno,
                           id.vars=1:7,
                           value.vars=8:ncol(anno),
                           variable.name="sample",
                           value.name="reads")
  anno.m$Position[anno.m$Strand=="-"] <- - anno.m$Position[anno.m$Strand=="-"]
  
  tss <- anno.m %>%
    group_by(sample, Position) %>%
    summarise(mean=mean(reads),
              sd=sd(reads),
              median=median(reads),
              mad=mad(reads))
  
  return(tss)
}

comparisons <- function(files,
                        title,
                        xlab="# of samples",
                        suffix="_peaks.narrowPeak") {
  # df <- lapply(files, read.delim, header=FALSE)
  if (is.character(files)) {
    names <- pipelineNGS::getNameFromPath(files, suffix=suffix)
    
    gr <- lapply(files, regioneR::toGRanges)
    names(gr) <- names
    gr <- GRangesList(gr)
  } else {
    gr <- GRangesList(files)
  }
  
  
  comb <- sapply(1:length(gr),
                 function(x) utils::combn(names(gr), m=x))
  
  comb.df.all <- data.frame()
  for (l in 1:length(comb)) {
    for (c in 1:ncol(comb[[l]])) {
      samples <- comb[[l]][,c]
      
      # Merge regions
      sel <- gr[names(gr) %in% samples]
      sel <- unlist(sel)
      sel <- regioneR::joinRegions(sel)
      
      # DF
      comb.df <- data.frame(num_islets=as.factor(l),
                            id_islets=paste0(samples, collapse=" "),
                            num_peaks=length(sel))
      comb.df.all <- rbind(comb.df.all, comb.df)
    }
  }
  
  
  stats <- comb.df.all %>%
    group_by(num_islets) %>%
    summarise(mean=mean(num_peaks),
              sd=sd(num_peaks),
              median=median(num_peaks),
              mad=mad(num_peaks))
  
  res <- list(df=comb.df.all,
              stats=stats)
  
  ggplot() +
    geom_line(data=res$stats,
              aes(num_islets, mean, group=1), size=0.7, linetype=2) +
    geom_jitter(data=res$df,
                aes(num_islets, num_peaks,
                    fill=num_islets), width=0.1, pch=21,
                color="black", size=2) +
    geom_pointrange(data=res$stats,
                    aes(x=num_islets, y=mean,
                        ymin=mean-sd, ymax=mean+sd), size=1) +
    xlab(xlab) +
    scale_y_continuous(label=comma, 
                       breaks=scales::pretty_breaks(),
                       name="# of peaks") +
    theme(legend.position="none") +
    ggtitle(title)
  
}

