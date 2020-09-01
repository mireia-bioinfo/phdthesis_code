load("T1D_risk_loci_granges.rda")
t1d <- loci
load("T2D_risk_loci_granges.rda")
t2d <- loci

ols <- findOverlaps(t1d, t2d)
t1d.cmn <- unique(t1d$topSNP[queryHits(ols)])
t2d.cmn <- unique(t2d$topSNP[subjectHits(ols)])

loci <- t1d[!(t1d$topSNP %in% t1d.cmn),]
save(loci,
     file="T1Drm_risk_loci_granges.rda")

loci <- t2d[!(t2d$topSNP %in% t2d.cmn),]
save(loci,
     file="T2Drm_risk_loci_granges.rda")

## Proxies ---------------------------------------
library(VSE)

# T1D
load("T1D_risk_proxies_granges.rda")
proxies <- proxies[!(proxies$locus %in% t1d.cmn),]
save(proxies,
     file="T1Drm_risk_proxies_granges.rda")

proxies <- unique(proxies[,c(1,8)])
colnames(mcols(proxies)) <- c("idLd", "idTag")
t1d.avs <- makeAVS(proxies)
save(t1d.avs, file="T1Drm_VSE_AVS.rda")

t1d.mrvs.500 <- makeMRVS(t1d.avs, bgSize=500, mc.cores = 8)
save(t1d.mrvs.500, file="T1Drm_VSE_MRVS_500.rda")

# T2D
source("../../../code/VSE_modFunctions.R")

load("T2D_risk_proxies_granges.rda")
proxies <- proxies[!(proxies$locus %in% t2d.cmn),]
save(proxies,
     file="T2Drm_risk_proxies_granges.rda")

proxies <- proxies[,c(1,8)]
colnames(mcols(proxies)) <- c("idLd", "idTag")
t2d.avs <- makeAVS(proxies)
save(t2d.avs, file="T2Drm_VSE_AVS.rda")

t2d.mrvs.500 <- makeMRVS.mod(t2d.avs, bgSize=500, mc.cores = 8)
save(t2d.mrvs.500, file="T2Drm_VSE_MRVS_500.rda")
