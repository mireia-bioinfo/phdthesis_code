## Load data from GWAS catalog  ---------------------------------------------
t1d <- read.delim("T1D_gwas-association-downloaded_2019-04-18-EFO_0001359-withChildTraits.tsv",
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

## Get proxies  --------------------------------------------------------------
getProxiesProcess <- function(chrom, 
                              pos, 
                              window_size=5e5, 
                              pop="EUR", 
                              rsquared=0.5) {
  message(paste0("> Query ", chrom, ":", pos))
  df <- proxysnps::get_proxies(chrom=chrom, 
                               pos=pos,
                                 window_size=window_size,
                                 pop=pop)
  ## Create DF proxies
  df$topSNP <- unique(df$ID[df$CHOSEN])[1]
  message(paste0("---  ", unique(df$ID[df$CHOSEN])[1]))
  df <- df[df$R.squared >= rsquared,]
  
  ## Create GRanges locus
  locus <- GRanges(paste0("chr", df$CHROM, ":",
                          min(df$POS, na.rm=T), "-",
                          max(df$POS, na.rm=T)))
  locus$topSNP <- unique(df$ID[df$CHOSEN])[1]
  
  info <- list(proxies=df,
               locus=locus)
  return(info)
}

info.all <- mapply(getProxiesProcess,
                   chrom=t1d$chr_name,
                   pos=t1d$chrom_start,
                   MoreArgs = list(window_size=2.5e5,
                                   rsquared=0.5),
                   SIMPLIFY = FALSE)
names(info.all) <- t1d$refsnp_id
save(info.all, file="T1D_info_proxies_all.rda")

## Create output data to work with --------------------------------------------
load("T1D_info_proxies_all.rda")

# T1D loci
loci <- unlist(GRangesList(lapply(info.all,
               function(x) x[["locus"]])))
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
save(loci, file="T1D_risk_loci_granges.rda")

# Proxies
proxies <- do.call(rbind,
                   lapply(info.all,
                          function(x) x[["proxies"]]))
proxies <- proxies[!is.na(proxies$POS),]
proxies <- regioneR::toGRanges(proxies[,c(1,2,2,3:10)])
seqlevelsStyle(proxies) <- "UCSC"

ols <- findOverlaps(proxies, loci)
proxies$locus <- loci$topSNP[subjectHits(ols)]
save(proxies, file="T1D_risk_proxies_granges.rda")

## Prepare data for VSE ---------------------------------
library(VSE)
load("T1D_risk_proxies_granges.rda")
proxies <- proxies[,c(1,8)]
colnames(mcols(proxies)) <- c("idLd", "idTag")

t1d.avs <- makeAVS(proxies)
save(t1d.avs, file="T1D_VSE_AVS.rda")

t1d.mrvs.500 <- makeMRVS(t1d.avs, bgSize=500, mc.cores = 8)
save(t1d.mrvs.500, file="T1D_VSE_MRVS_500.rda")

## Prepare data for VSE 2 --------------------------------
library(VSE)
load("T1D_risk_proxies_granges.rda")
proxies <- proxies[proxies$R.squared>=0.8,c(1,8)]
colnames(mcols(proxies)) <- c("idLd", "idTag")

t1d.avs <- makeAVS(proxies)
save(t1d.avs, file="T1D_VSE_AVS_ld0.8.rda")

t1d.mrvs.500 <- makeMRVS(t1d.avs, bgSize=500, mc.cores = 8)
save(t1d.mrvs.500, file="T1D_VSE_MRVS_500_ld0.8.rda")
