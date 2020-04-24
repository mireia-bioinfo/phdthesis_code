####################################
## Correlation 10Kb binned genome ##
####################################

## Parameters
args <-  commandArgs(trailingOnly=TRUE)
win <- 1e4
pathbam <- args[1]
pathQC <- args[1]

## Load libraries ------------------
library(DESeq2)
library("BiocParallel")

## Code -----------------------------
## Create tiled genome (5kb)
hg19 <- read.delim("~/data/hg19.len", header=F)
seqlengths <- hg19$V2
names(seqlengths) <- hg19$V1

hg19.tile <- unlist(tileGenome(seqlengths, tilewidth = win))
hg19.tile <- hg19.tile[seqnames(hg19.tile) %in% c(paste0("chr", 1:22), "chrX", "chrY"),]

saf <- data.frame(hg19.tile)[,1:3]
colnames(saf) <- c("Chr", "Start", "End")
saf$GeneID <- paste0("tile_", 1:nrow(saf))
saf$Strand <- "+"

# List bams
bam <- list.files(pathbam,
                  pattern="offset.bam$",
                  full.names=TRUE)

cov <- Rsubread::featureCounts(bam,
                               annot.ext=saf,
                               nthreads=10)

names <- pipelineNGS::getNameFromPath(bam, suffix=".offset.bam")
cov$counts <- cov$counts[rowSums(cov$counts)>0,]
mat <- cov$counts
colnames(mat) <- names

# Normalize with DESeq2
coldata <- data.frame(ids=names,
                      sample=unlist(lapply(strsplit(names, "_"),
                                           function(x) x[1])),
                      treatment=unlist(lapply(strsplit(names, "_"),
                                              function(x) x[2])))
coldata$batch <- gsub("hi", "", gsub("endoc", "", coldata$sample))
coldata$tissue <- gsub("[[:digit:]]*", "", coldata$sample)

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = coldata,
                              design= ~ tissue + treatment)


register(MulticoreParam(10))

dds <- DESeq(dds,
             parallel=TRUE,
             BPPARAM=MulticoreParam(10))

counts <- counts(dds, normalized=TRUE)
save(counts, file=paste0(pathQC, "COR_10kb_norm.rda"))
