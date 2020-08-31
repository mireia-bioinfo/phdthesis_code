library(umi4cPackage)
library(GenomicRanges)
devtools::load_all("~/tools/umi4cCatcheR/")

w <- 3e3

# Load IREs -------------------------
load("../../data/REs/REs_hi_fc1_padj0.05_granges_subgroup.rda")
re <- re[!grepl("Other", re$type),]

# Load UMIs -------------------------
prom <- c("SOCS1", "DOCK9", "IFIH1", "GBP1", "CIITA", "TNFSF10", "ATF3", "IRF1",
          "CMPK2", "LAMP3", "RSAD2", "CXCL11", "IFIT1", "DEXI")

df <- read.delim("data/UMI4C_promoters_views.tsv", stringsAsFactors=F, header=T)
df <- df[df$bait %in% prom,]

conf <- "../../data/UMI4C/merged/conf/"
umi4cPackage::p4cLoadConfFiles(conf)
gtracks <- gtrack.ls()
gtracks <- gtracks[unlist(sapply(prom, grep, gtracks))]
gtracks <- gtracks[grepl("_m_", gtracks)]

for (i in 1:nrow(df)) {
  load(paste0("data/UMI4C_norm_results_", df$bait[i], ".rda"))
  win <- res$view
  
  ## Create p4cObjects
  gt.sel <- gtracks[grepl(df$bait[i], gtracks)]
  ctrl <- suppressMessages(p4cNewProfile(gt.sel[grep("ctrl", gt.sel)],
                                         scope_5=5e5, scope_3=5e5))
  
  cyt <- suppressMessages(p4cNewProfile(gt.sel[grep("cyt", gt.sel)],
                                        scope_5=5e5, scope_3=5e5))
  
  ## Load REs
  re.sel <- subsetByOverlaps(re, win)
  re.sel <- resize(re.sel, width=w, fix="center")
  re.sel <- re.sel[order(re.sel)]
  
  ## Perform tests
  test <- suppressMessages(mapply(p4cIntervalsMean,
                                  start=start(re.sel),
                                  end=end(re.sel),
                                  MoreArgs=list(p4c_obj=cyt, 
                                                ref_p4c_obj=ctrl,
                                                min_win_cov=df$min_win_mols[i]),
                                  SIMPLIFY=FALSE))
  test <- do.call(rbind, test)
  colnames(test)[3:4] <- c("cyt", "ctrl")
  test <- cbind(data.frame(re.sel)[,-c(4:5)],
                test[,-c(1:2)])
  
  save(test,
       file=paste0("data/UMI4C_test_promoters_", df$bait[i], ".rda"))
}