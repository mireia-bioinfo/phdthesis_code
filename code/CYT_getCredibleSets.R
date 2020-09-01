uk <- list.files("/imppc/labs/lplab/mramos/GWAS_T1D_JT/version_ext/",
                 pattern=".txt", full.names=TRUE)
uk <- uk[!grepl("sardinia",uk)]
names <- pipelineNGS::getNameFromPath(uk, suffix="_lorenzo_resized.txt")

uk <- lapply(uk, read.delim, stringsAsFactors=FALSE)
uk <- lapply(uk, function(x) regioneR::toGRanges(x[,c(3,4,4,1:2,5:12)]))
names(uk) <- names

loci <- read.delim("~/Projects/EndoC_CYT_hg19/results_paper/data05/T1D_riskLociWithInducedRE.txt")
loci$seqnames <- gsub("chr", "", loci$seqnames)
loci <- regioneR::toGRanges(loci)

t <- lapply(uk, findOverlaps, loci)

lens.ol <- sapply(t, length)
lens.snps <- sapply(uk, length)

uk.df <- lapply(uk, data.frame)
uk.long <- do.call(rbind, uk.df)
uk.long$riskLocus <- unlist(mapply(rep, names(uk), each=lens.snps, SIMPLIFY=TRUE))
save(uk.long, file="data/ukDataset.rda")

## Perform credible set analysis with available loci
load("data/ukDataset.rda")

## filter maf > 0.01
uk <- uk.long[uk.long$controls_maf_affy>0.01 & uk.long$controls_maf_illu>0.01,]

## snps with pval < 5e-8
table(uk.long$riskLocus, uk.long$pmeta<5e-8)
table(uk.long$riskLocus, uk.long$pmeta<1e-4) # sec ond threshold

library(dplyr)
leading <- uk.long %>%
  group_by(riskLocus) %>%
  summarise(minPval=min(pmeta),
            leading=ID[minPval==pmeta])

uk.long <- dplyr::left_join(uk.long, leading[,c(1:3)])

## Get r2 values for pairs of leading-proxies
uk.long$r2 <- NA
uk.long <- uk.long[uk.long$riskLocus=="riskLocus_4",]

for (i in 2800:2900){#nrow(uk.long)) {
  cmd <- paste("curl -k -X GET",
               shQuote(paste0("https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldpair?var1=",uk.long$leading[i],
                              "&var2=", uk.long$ID[i],
                              "&pop=CEU%2BTSI%2BFIN%2BGBR%2BIBS&token=bdb73ad74ff8")))
  r2 <- system(cmd, intern = TRUE)
  r2 <- r2[grep("R2:", r2)]
  r2 <- gsub("          R2: ", "", r2)
  if (length(r2)==0) r2 <- NA
  uk.long$r2[i] <- as.numeric(r2)
}

save(uk.long, file="data/ukDataset_filtMAF_withR2.rda")

load("data/ukDataset_filtMAF_withR2.rda")

uk.filt <- uk.long[uk.long$r2>0.1,]
save(uk.filt, file="data/ukDataset_filtMAF_filtR2.rda")

load("data/ukDataset_filtMAF_filtR2.rda")
uk.filt$top <- uk.filt$leading == uk.filt$ID

credSet <- uk.filt %>%
  mutate(r=0.04/(seall*seall + 0.04),
         z=beta/seall,
         abf=sqrt(1-r)/exp(-r*z*z/2)) %>%
  group_by(riskLocus) %>%
  arrange(desc(abf), .by_group=TRUE) %>%
  group_by(riskLocus) %>%
  mutate(PP=abf/sum(abf),
         cumsum=cumsum(PP),
         inCredible=(cumsum<=0.99 | top))

save(credSet, file="data/uk_credibleSet.rda")

table(credSet$riskLocus, credSet$inCredible)

load("data/ukDataset_filtMAF_withR2.rda")
uk.long$r2 <- as.numeric(uk.long$r2)
uk.long$seqnames <- paste0("chr", uk.long$seqnames)

uk.long <- dplyr::left_join(uk.long, credSet)
uk.long$inCredible[is.na(uk.long$inCredible)] <- FALSE

uk.full <- uk.long
save(uk.full, file="data/uk_credibleSet_fullr2.rda")
