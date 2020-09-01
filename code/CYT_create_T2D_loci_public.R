suppressPackageStartupMessages(library(GenomicRanges))

## Load data from GWAS catalog  ---------------------------------------------
t2d <- read.delim("raw/T2D_gwas-association-downloaded_2019-04-18-EFO_0001360-withChildTraits.tsv",
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

## Get proxies  --------------------------------------------------------------
# getProxiesProcess <- function(chrom, 
#                               pos, 
#                               window_size=5e5, 
#                               pop="EUR", 
#                               rsquared=0.5) {
#   message(paste0("> Query ", chrom, ":", pos))
#   if (!file.exists(paste0("tmp_", chrom, ":", pos, ".rda"))) {
#     
#     df <- proxysnps::get_proxies(chrom=chrom, 
#                                  pos=pos,
#                                  window_size=window_size,
#                                  pop=pop)
#     ## Create DF proxies
#     df$topSNP <- unique(df$ID[df$CHOSEN])[1]
#     message(paste0("---  ", unique(df$ID[df$CHOSEN])[1]))
#     df <- df[df$R.squared >= rsquared,]
#     
#     ## Create GRanges locus
#     if (nrow(df)>0) {
#       locus <- GRanges(paste0("chr", df$CHROM, ":",
#                               min(df$POS, na.rm=T), "-",
#                               max(df$POS, na.rm=T)))
#       locus$topSNP <- unique(df$ID[df$CHOSEN])[1]
#       
#       info <- list(proxies=df,
#                    locus=locus)
#       save(info, file=paste0("tmp_", chrom, ":", pos, ".rda"))
#     }
#   } else {
#     message("--- File exists")
#   }
# }
# 
# info.all <- mapply(getProxiesProcess,
#                    chrom=t2d$chr_name,
#                    pos=t2d$chrom_start,
#                    MoreArgs = list(window_size=2.5e5,
#                                    rsquared=0.5),
#                    SIMPLIFY = FALSE)
# 
# files <- list.files(".",
#                     pattern="tmp_")
# info.all <- list()
# for (i in files) {
#   load(i)
#   info.all[[gsub("tmp_", "", gsub(".rda", "", i))]] <- info[!is.na(info)]
# }
# save(info.all, file="T2D_info_proxies_all.rda")

## Create output data to work with --------------------------------------------
load("T2D_info_proxies_all.rda")

# T2D loci
loci <- unlist(GRangesList(lapply(info.all,
               function(x) x[["locus"]])))
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
save(loci, file="T2D_risk_loci_granges.rda")

# Proxies
proxies <- do.call(rbind,
                   lapply(info.all,
                          function(x) x[["proxies"]]))
proxies <- proxies[!is.na(proxies$POS),]
proxies <- regioneR::toGRanges(proxies[,c(1,2,2,3:10)])
seqlevelsStyle(proxies) <- "UCSC"

ols <- findOverlaps(proxies, loci)
proxies$locus <- loci$topSNP[subjectHits(ols)]
save(proxies, file="T2D_risk_proxies_granges.rda")

## Prepare data for VSE ---------------------------------
library(VSE)
load("T2D_risk_proxies_granges.rda")
proxies <- proxies[,c(1,8)]
colnames(mcols(proxies)) <- c("idLd", "idTag")
proxies <- proxies[!is.na(proxies$idTag),]

t2d.avs <- makeAVS(proxies)
save(t2d.avs, file="T2D_VSE_AVS.rda")

source("../../../code/VSE_modFunctions.R")
t2d.mrvs.500 <- makeMRVS.mod(t2d.avs, bgSize=500, mc.cores = 8)
save(t2d.mrvs.500, file="T2D_VSE_MRVS_500.rda")
