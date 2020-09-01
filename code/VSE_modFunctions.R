variantSetEnrichment.mod <- function(avs, 
                                     mrvs, 
                                     regions,
                                     maxgap=-1L) 
{
  no_of_tags <- length(avs)
  no_of_beds <- length(regions$SampleID)
  intersect_matrix <- matrix(NA, nrow = no_of_beds, ncol = 1)
  row.names(intersect_matrix) <- regions$SampleID
  overlap <- matrix(NA, nrow = no_of_beds, ncol = no_of_tags)
  beds_list <- list()
  message("Loading regions")
  for (i in 1:no_of_beds) {
    bed_path <- as.character(regions$Peaks[i])
    bed.gr <- bedToGRanges(bed_path)
    beds_list <- c(beds_list, bed.gr)
    overlap[i, ] <- ifelse(countOverlaps(avs, bed.gr, maxgap=maxgap, 
                                         ignore.strand = TRUE) > 
                             0, 1, 0)
    message(paste0(sum(overlap[i, ]), " LD blocks intersect with ", 
                   regions$SampleID[i]))
  }
  intersect_matrix[, 1] <- rowSums(overlap)
  null_intersects <- matrix(NA, nrow = no_of_beds, ncol = length(mrvs))
  row.names(null_intersects) <- regions$SampleID
  message("Tallying MRVS ", "\r", appendLF = FALSE)
  for (i in 1:length(mrvs)) {
    flush.console()
    imod <- i%%length(mrvs)
    message(ifelse(imod%%10 == 0, i, "."), "\r", appendLF = FALSE)
    overlap <- matrix(NA, nrow = no_of_beds, ncol = length(mrvs[[i]]))
    for (k in 1:length(beds_list)) {
      overlap[k, ] <- ifelse(countOverlaps(mrvs[[i]], 
                                           beds_list[[k]], 
                                           maxgap=maxgap,
                                           ignore.strand = TRUE) > 0, 1, 
                             0)
    }
    null_intersects[, i] <- rowSums(overlap)
  }
  intersect_matrix <- cbind(intersect_matrix, null_intersects)
  normalityCutoff = 1
  matrix_norm <- t(apply(intersect_matrix, 1, function(x) x + 
                           runif(length(x), 0, 1)))
  ksvalues <- apply(matrix_norm, 1, function(x) {
    kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x), 
                   exact = TRUE)
    kst$p.value
  })
  message("\nNormalizing null distribution")
  for (i in 1:length(ksvalues)) {
    if (ksvalues[i] < normalityCutoff) {
      ks_lambda <- data.frame(lambda = seq(-2, 2, 0.1), 
                              ksp = rep(0, length(seq(-2, 2, 0.1))))
      for (j in 1:length(ks_lambda$lambda)) {
        x <- car::bcPower(matrix_norm[i, ], ks_lambda$lambda[j])
        kst <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x), 
                       exact = TRUE)
        ks_lambda$ksp[j] <- kst$p.value
      }
      ks_lambda <- ks_lambda[order(ks_lambda$ksp, decreasing = TRUE), 
                             ]
      matrix_norm[i, ] <- car::bcPower(matrix_norm[i, 
                                                   ], ks_lambda[1, 1])
      ksvalues[i] <- ks_lambda[1, 2]
    }
  }
  matrix_scaled <- t(apply(matrix_norm, 1, function(x) x <- (x - 
                                                               median(x))/sd(x)))
  vse_matrices <- list(intersect_matrix, matrix_norm, matrix_scaled)
  return(vse_matrices)
}


VSESummary.mod <- function (data,
                            method="bonferroni") 
{
  p_values <- c()
  intersect_matrix <- data[[1]]
  matrix_scaled <- data[[3]]
  for (i in 1:nrow(matrix_scaled)) {
    avs <- matrix_scaled[i, 1]
    null <- matrix_scaled[i, c(2:ncol(matrix_scaled))]
    p <- 1 - pnorm(matrix_scaled[i, ], mean(matrix_scaled[i, 
                                                          ]), sd(matrix_scaled[i, ]))
    p_values[i] <- p[1]
  }
  padjust_values <- p.adjust(p_values, method = method)
  summary_df <- data.frame(region = row.names(matrix_scaled), 
                           avs = intersect_matrix[, 1], enrichment = matrix_scaled[, 
                                                                                   1], p.value = p_values, p.adjusted = padjust_values)
  return(summary_df)
}

VSEplot.mod <- function (data, 
                         padj = 0.01, 
                         method="bonferroni", 
                         ...) 
{
  matrix_norm <- data[[2]]
  matrix_scaled <- data[[3]]
  p_values <- c()
  for (i in 1:nrow(matrix_scaled)) {
    avs <- matrix_scaled[i, 1]
    null <- matrix_scaled[i, c(2:ncol(matrix_scaled))]
    p <- 1 - pnorm(matrix_scaled[i, ], mean(matrix_scaled[i, 
                                                          ]), sd(matrix_scaled[i, ]))
    p_values[i] <- p[1]
  }
  padjust_values <- p.adjust(p_values, method = method)
  boxplot(t(matrix_scaled), ylim = c(min(matrix_scaled), max(matrix_scaled)), 
          ylab = "Enrichment Score", outline = FALSE, notch = TRUE, 
          ...)
  for (i in 1:length(padjust_values)) {
    points(i, matrix_scaled[i, 1], pch = 19, col = ifelse(padjust_values[i] <= 
                                                            padj, "red", "black"))
  }
}

makeMRVS.mod <- function (avs, bgSize = 100, mc.cores = 6) 
{
  if (!exists("nullblocks.08")) {
    nullblocks.08 <- list()
    tmpfile <- tempfile(fileext = ".rda")
    url <- "http://www.hansenhelab.org/null0.8.rda"
    download.file(url, destfile = tmpfile, method = "curl")
    load(tmpfile)
  }
  if (!exists("nullblocks.08")) {
    stop("Downloading null failed. Are you connected to internet?")
  }
  no_of_tags <- length(avs)
  tally <- nullblocks.08$tally
  null_id_list <- matrix(NA, nrow = no_of_tags, ncol = bgSize + 
                           1)
  for (i in 1:no_of_tags) {
    null_id_list[i, 1] <- as.character(elementMetadata(avs[[i]])[1, 
                                                                 2])
    ld_tally <- length(avs[[i]])
    tlist <- tally[tally$X0.8 == ld_tally, 1]
    if (length(tlist)==0) {
      ld_tally <- unique(tally$X0.8[which(abs(tally$X0.8-ld_tally)==min(abs(tally$X0.8-ld_tally)))])
      tlist <- tally[tally$X0.8 == ld_tally, 1]
      null_id_list[i, c(2:ncol(null_id_list))] <- as.character(tlist[sample(length(tlist), 
                                                                            bgSize, 
                                                                            replace = ifelse(length(tlist) < bgSize, 
                                                                                                     TRUE, FALSE))])
    } else {
      null_id_list[i, c(2:ncol(null_id_list))] <- as.character(tlist[sample(length(tlist), 
                                                                            bgSize, 
                                                                            replace = ifelse(length(tlist) < bgSize, 
                                                                                                     TRUE, FALSE))])
    }
  }
  
  message(paste0("Using ", mc.cores, " cores"))
  mrvs <- list()
  for (i in 2:ncol(null_id_list)) {
    message(paste0("Computing MRVS no. ", i - 1))
    null_glist <- GRangesList(parallel::mclapply(null_id_list[,i], 
                                                 function(x) {
                                                                chr <- strsplit(x, ":")[[1]][1]
                                                                pos <- strsplit(x, ":")[[1]][2]
                                                                var <- paste0("nullblocks.08", "$", chr)
                                                                chr <- as.factor(unlist(lapply(chr, function(z) gsub("chr", 
                                                                                                                     "", z))))
                                                                lds <- eval(parse(text = var))[eval(parse(text = var))$BP_A %in% 
                                                                                                 pos, 1]
                                                                gr <- GRanges(seqnames = chr, ranges = IRanges::IRanges(start = lds, 
                                                                                                                        width = 1))
                                                              }, 
                                                 mc.cores = mc.cores))
    mrvs <- c(mrvs, null_glist)
  }
  return(mrvs)
}

