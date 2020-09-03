## Parameters -----------------------------------------
scope <- 2e3 # Scope for TSS enrichment
bin <- 10 # bin for TSS enrichment

## Libraries ------------------------------------------
setwd("data/CYT/ATAC/")

library(pipelineNGS)
library(dplyr)
library(GenomicRanges)
source("../../code/CYT_QC_ATAC_functions.R")

## Stats ----------------------------------------------
bams <- list.files("BAM",
                   pattern=".raw.bam$",
                   full.names=TRUE)
stats <- getStats(raw_bam=bams, path_logs = "Logs/")
stats$factor <- 1e7/stats$final_reads

save(stats, file="QC/ATAC_stats.rda")

## Enrichment at TSS -----------------------------------
# Get transcripts ---------
library(biomaRt)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position",
                            "strand", "ensembl_gene_id", "external_gene_name"),
               filters="biotype", values="protein_coding",
               mart=grch37)

genes$strand[genes$strand==-1] <- "-"
genes$strand[genes$strand==1] <- "+"

genes <- regioneR::toGRanges(genes)
strand(genes) <- genes$strand
mcols(genes) <- mcols(genes)[,-1]
genes <- genes[seqnames(genes) %in% c(1:22, "X", "Y"),]
seqlevelsStyle(genes) <- "UCSC"
genes$GeneID <- genes$ensembl_gene_id

# Extend regions to scope*2
regions <- unique(promoters(genes, upstream=0, downstream=1))
regions.bin <- binRegions(regions,
                          scope=scope,
                          bin=bin)

# Create saf for annotating regions
saf <- data.frame(regions.bin)[,c(1:3,5,6,8)]
colnames(saf) <- c("Chr", "Start", "End", "Strand", "txID", "Position")
saf$GeneID <- paste0(as.character(saf$txID), "_", saf$Position)
save(saf, file="QC/ATAC_tss_saf.rda")

# Create random saf
rnd <- regioneR::randomizeRegions(genes)
rnd$GeneID <- 1:length(rnd)

regions <- unique(promoters(rnd, upstream=0, downstream=1))
regions.bin <- binRegions(regions,
                          scope=scope,
                          bin=bin)

saf <- data.frame(regions.bin)[,c(1:3,5,6,8)]
colnames(saf) <- c("Chr", "Start", "End", "Strand", "txID", "Position")
saf$GeneID <- paste0(as.character(saf$txID), "_", saf$Position)
save(saf, file="QC/ATAC_tss_random_saf.rda")

## Annotate reads to regions ------------------------
files <- list.files("BAM",
                    pattern="*offset.bam$",
                    full.names=TRUE)
load("QC/ATAC_tss_saf.rda")

tss <- getTSSenrichment(bam=files,
                 saf=saf)

## Random regions
load("QC/ATAC_tss_random_saf.rda")

tss.rnd <- getTSSenrichment(bam=files,
                            saf=saf)

tss$dataset <- "TSS annotation"
tss.rnd$dataset <- "Random control"

tss <- rbind(tss, tss.rnd)
save(tss,
     file="QC/ATAC_tss_enrichment.rda")

## % of reads in peaks -----------------------------------
# For each bam load peaks and get number of assigned reads
files.bam <- list.files("BAM",
                          pattern=".offset.bam$",
                          full.names=TRUE)

files.peaks <- gsub(".offset.bam", "_peaks.narrowPeak", files.bam)
files.peaks <- gsub("BAM/", "Peaks/", files.peaks)

names <- pipelineNGS::getNameFromPath(files.bam,
                suffix=".offset.bam")

df <- data.frame()
for (i in 1:length(files.bam)) {
  ## Annotate as distal or proximal
  peaks <- regioneR::toGRanges(files.peaks[i])
  anno <- data.frame(ChIPseeker::annotatePeak(peaks,
                                              TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))
  anno$annotation <- "Promoter"
  anno$annotation[abs(anno$distanceToTSS)>2e3] <- "Distal"

  saf <- anno[,c(1:3,6,12)]
  colnames(saf) <- c("Chr", "Start", "End", "GeneID", "Annotation")
  saf$Strand <- "+"

  counts <- Rsubread::featureCounts(files.bam[i],
                                    annot.ext=saf,
                                    nthreads=8)

  anno <- cbind(saf, counts$counts)
  colnames(anno)[ncol(anno)] <- "reads"

  sum <- anno %>%
    group_by(Annotation) %>%
    summarise(reads=sum(reads))
  sum[3,1] <- "Unassigned"
  sum[3,2] <- sum(counts$stat[grepl("Unas", counts$stat$Status),2])

  sum$sample <- names[i]
  df <- rbind(df, sum)

  ## With randomized peaks
  peaks <- regioneR::randomizeRegions(peaks)
  peaks$id <- paste0("peak_", 1:length(peaks))
  anno <- data.frame(ChIPseeker::annotatePeak(peaks,
                                              TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))
  anno$annotation <- "Promoter"
  anno$annotation[abs(anno$distanceToTSS)>2e3] <- "Distal"

  saf <- anno[,c(1:3,6,7)]
  colnames(saf) <- c("Chr", "Start", "End", "GeneID", "Annotation")
  saf$Strand <- "+"

  counts <- Rsubread::featureCounts(files.bam[i],
                                    annot.ext=saf,
                                    nthreads=8)

  anno <- cbind(saf, counts$counts)
  colnames(anno)[ncol(anno)] <- "reads"

  sum <- anno %>%
    group_by(Annotation) %>%
    summarise(reads=sum(reads))
  sum[3,1] <- "Unassigned"
  sum[3,2] <- sum(counts$stat[grepl("Unas", counts$stat$Status),2])

  sum$sample <- paste0(names[i], "_rndm")
  df <- rbind(df, sum)
}

df$Annotation <- factor(df$Annotation,
                        levels=c("Unassigned", "Distal", "Promoter"))

df$group <- gsub("[[:digit:]]*_", "_", df$sample)
df$group[grep("rndm", df$group)] <- "Random"


## Convert to %
total <- df %>%
  group_by(sample) %>%
  summarise(total=sum(reads))

df <- dplyr::left_join(df, total)
df$percent <- df$reads/df$total*100

stats <- df %>%
  group_by(group, Annotation) %>%
  summarise(mean=mean(percent),
            sd=sd(percent))
save(stats, file = "QC/ATAC_noise.rda")

## PhantomPeak QC ----------------------------------------------
bam <- list.files("BAM/",
                  pattern="offset.bam$",
                  full.names=T)
bam <- bam[grep("hi44", bam)]
path_script <- "~/tools/phantompeakqualtools/run_spp.R"
dir.create("tmp/", F)

for (i in bam) {

  name <- pipelineNGS::getNameFromPath(i, suffix=".offset.bam")
  cmd <- paste("Rscript-3.5.1-bioc-3.8",
               path_script,
               paste0("-c='", i, "'"),
               "-p=5 -tmpdir='tmp/'",
               "-savp -rf",
               paste0("-out='QC/QC_phantompeakqualtools_", name, ".txt'"))

  system(cmd)
}

# Load files
txt <- list.files("QC/",
                  pattern="QC_phantompeakqualtools_",
                  full.names=T)

txt <- lapply(txt, read.delim, stringsAsFactors=F, header=F)
txt <- do.call(rbind, txt)

colnames(txt) <- c("sampleID", "numReads", "estFragLen", "corr_estFragLen", "phantomPeak",
                   "corr_phantomPeak", "argmin_corr", "min_corr", "NSC", "RSC", "QualityTag")

txt$sampleID <- gsub(".offset.bam", "", txt$sampleID)

save(txt, file="QC/QC_scores.rda")
