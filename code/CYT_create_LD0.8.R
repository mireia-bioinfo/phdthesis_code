library(GenomicRanges)
library(VSE)

######################################################################
## T1D LD 0.8 --------------------------------------------------------
######################################################################
t1d <- read.delim("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/raw/T1D_gwas-association-downloaded_2019-04-18-EFO_0001359-withChildTraits.tsv",
                 header=T, stringsAsFactors=T)

t1d <- unique(t1d[t1d$DISEASE.TRAIT=="Type 1 diabetes", c(22,15)])
t1d <- rbind(t1d,
             data.frame(SNPS="rs78037977",
                        MAPPED_GENE="FASLG"))

## Get genomic coordinates --------------------------------------------------
library(biomaRt)
snp37 <- useMart(biomart="ENSEMBL_MART_SNP", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_snp")

anno <- getBM(attributes = c("chr_name", "chrom_start", "chrom_end", "refsnp_id"),
              filters = "snp_filter",
              values = t1d$SNPS,
              mart = snp37)

t1d <- dplyr::left_join(anno, t1d,
                        by=c("refsnp_id"="SNPS"))
t1d <- t1d[t1d$chr_name %in% c(1:22, "X", "Y"),]

##-------------------------------------------------------------------------------------
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1D_info_proxies_all.rda")

## Create loci with LD 0.8
loci <- lapply(info.all,
               function(x) x[["proxies"]])

loci <- lapply(loci,
               function(x) x[x$R.squared>=0.8,])

loci <- lapply(loci,
               function(x) GRanges(paste0("chr", x$CHROM, ":",
                                   min(x$POS, na.rm=T), "-",
                                   max(x$POS, na.rm=T))))
loci <- unlist(GRangesList(loci))
loci <- regioneR::joinRegions(loci)
loci <- loci[seqnames(loci)!="chrNA",]

t1d.gr <- regioneR::toGRanges(t1d)
seqlevelsStyle(t1d.gr) <- "UCSC"

ols <- findOverlaps(loci, t1d.gr)
snps <- split(t1d.gr$refsnp_id, queryHits(ols))
snps <- sapply(snps, paste0, collapse=", ")
loci$topSNP <- snps

gene <- split(t1d.gr$MAPPED_GENE, queryHits(ols))
gene <- sapply(gene, 
               function(x) paste0(unique(x), collapse=", "))
loci$mappedGene <- gene
save(loci, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1D_LD0.8_risk_loci_granges.rda")

# Proxies
proxies <- do.call(rbind,
                   lapply(info.all,
                          function(x) x[["proxies"]]))
proxies <- proxies[!is.na(proxies$POS),]
proxies <- proxies[proxies$R.squared>=0.8,]
proxies <- regioneR::toGRanges(proxies[,c(1,2,2,3:10)])
seqlevelsStyle(proxies) <- "UCSC"

ols <- findOverlaps(proxies, loci)
proxies$locus <- loci$topSNP[subjectHits(ols)]
save(proxies, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1D_LD0.8_risk_proxies_granges.rda")

######################################################################
## T2D LD 0.8 --------------------------------------------------------
######################################################################
t2d <- read.delim("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/raw/T2D_gwas-association-downloaded_2019-04-18-EFO_0001360-withChildTraits.tsv",
                  header=T, stringsAsFactors=T)

t2d <- unique(t2d[t2d$DISEASE.TRAIT=="Type 2 diabetes", c(22,15)])


## Get genomic coordinates --------------------------------------------------
suppressPackageStartupMessages(library(biomaRt))
snp37 <- useMart(biomart="ENSEMBL_MART_SNP", 
                 host="grch37.ensembl.org", 
                 path="/biomart/martservice", 
                 dataset="hsapiens_snp")

anno <- getBM(attributes = c("chr_name", "chrom_start", "chrom_end", "refsnp_id"),
              filters = "snp_filter",
              values = t2d$SNPS,
              mart = snp37)

t2d <- dplyr::left_join(anno, t2d,
                        by=c("refsnp_id"="SNPS"))
t2d <- t2d[t2d$chr_name %in% c(1:22, "X", "Y"),]

## Create loci with LD 0.8
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2D_info_proxies_all.rda")
loci <- lapply(info.all,
               function(x) x[["proxies"]])

loci <- lapply(loci,
               function(x) x[x$R.squared>=0.8,])

loci <- lapply(loci,
               function(x) GRanges(paste0("chr", x$CHROM, ":",
                                          min(x$POS, na.rm=T), "-",
                                          max(x$POS, na.rm=T))))
loci <- unlist(GRangesList(loci))
loci <- regioneR::joinRegions(loci)
loci <- loci[seqnames(loci)!="chrNA",]

t2d.gr <- regioneR::toGRanges(t2d)
seqlevelsStyle(t2d.gr) <- "UCSC"

ols <- findOverlaps(loci, t2d.gr)
snps <- split(t2d.gr$refsnp_id[subjectHits(ols)], queryHits(ols))
snps <- sapply(snps, paste0, collapse=", ")
loci$topSNP[unique(queryHits(ols))] <- snps

gene <- split(t2d.gr$MAPPED_GENE[subjectHits(ols)], queryHits(ols))
gene <- sapply(gene, 
               function(x) paste0(unique(x), collapse=", "))
loci$mappedGene <- NA
loci$mappedGene[unique(queryHits(ols))] <- gene
save(loci, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2D_LD0.8_risk_loci_granges.rda")

# Proxies
proxies <- do.call(rbind,
                   lapply(info.all,
                          function(x) x[["proxies"]]))
proxies <- proxies[!is.na(proxies$POS),]
proxies <- proxies[proxies$R.squared>=0.8,]
proxies <- regioneR::toGRanges(proxies[,c(1,2,2,3:10)])
seqlevelsStyle(proxies) <- "UCSC"

ols <- findOverlaps(proxies, loci)
proxies$locus <- loci$topSNP[subjectHits(ols)]
save(proxies, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2D_LD0.8_risk_proxies_granges.rda")

######################################################################
## Remove common -----------------------------------------------------
######################################################################
## Remove common with T2D
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1D_LD0.8_risk_loci_granges.rda")
t1d <- loci
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2D_LD0.8_risk_loci_granges.rda")
t2d <- loci

ols <- findOverlaps(t1d, t2d)
t1d.cmn <- unique(t1d$topSNP[queryHits(ols)])
t2d.cmn <- unique(t2d$topSNP[subjectHits(ols)])

loci <- t1d[!(t1d$topSNP %in% t1d.cmn),]
save(loci,
     file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1Drm_LD0.8_risk_loci_granges.rda")

loci <- t2d[!(t2d$topSNP %in% t2d.cmn),]
save(loci,
     file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2Drm_LD0.8_risk_loci_granges.rda")

source("../../code/VSE_modFunctions.R")

# T1D
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1D_LD0.8_risk_proxies_granges.rda")
proxies <- proxies[!(proxies$locus %in% t1d.cmn),]
save(proxies,
     file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1Dr_LD0.8m_risk_proxies_granges.rda")

proxies <- unique(proxies[,c(1,8)])
colnames(mcols(proxies)) <- c("idLd", "idTag")
t1d.avs <- makeAVS(proxies)
save(t1d.avs, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1Drm_LD0.8_VSE_AVS.rda")

t1d.mrvs.500 <- makeMRVS(t1d.avs, bgSize=500, mc.cores = 8)
save(t1d.mrvs.500, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T1Drm_LD0.8_VSE_MRVS_500.rda")

### T2D
load("~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2D_LD0.8_risk_proxies_granges.rda")
proxies <- proxies[!(proxies$locus %in% t2d.cmn),]
save(proxies,
     file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2Drm_LD0.8_risk_proxies_granges.rda")

proxies <- proxies[,c(1,8)]
colnames(mcols(proxies)) <- c("idLd", "idTag")
t2d.avs <- makeAVS(proxies)
save(t2d.avs, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2Drm_LD0.8_VSE_AVS.rda")

t2d.mrvs.500 <- makeMRVS.mod(t2d.avs, bgSize=500, mc.cores = 8)
save(t2d.mrvs.500, file="~/Projects/CYT_hg19/data/T1D_loci_JT/T1D_public_loci/T2Drm_LD0.8_VSE_MRVS_500.rda")
