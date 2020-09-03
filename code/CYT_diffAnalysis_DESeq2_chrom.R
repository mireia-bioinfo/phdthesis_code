#####################################
## Differential analysis Chromatin ##
#####################################

## Create options for script ---------------------------------------
library("optparse")

option_list = list(
  make_option(c("-f", "--foldChange"), type="numeric", default=1,
              help="log2FC of the RE dataset to use", metavar="numeric"),
  make_option(c("-q", "--padj"), type="numeric", default=0.05,
              help="Adjuste P-value cuttoff of the RE dataset to use", metavar="numeric"),
  make_option(c("-p", "--promoter"), type="numeric", default=2e3,
              help="Distance to TSS for classifying regions as promoters. [default= %default]", metavar="numeric"),
  make_option(c("-b", "--batch"), type="logical", default=T,
              help="Correct for batch effect in differential analysis [default= %default]",
              metavar="logical"),
  make_option(c("-s", "--sample"), type="character", default="endoc",
              help="Sample to use for this analysis. [default= %default]",
              metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=4,
              help="Number of cores to use for this analysis. [default= %default]",
              metavar="numeric"),
  make_option(c("-e", "--experiment"), type="character", default="ATAC",
              help="Name of the experiment",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Set parameters --------------------------------------
dir.create("diffAnalisis/", F)
if (opt$batch) b <- "Batch" else b <- ""

## Load libraries --------------------------------------
library(GenomicRanges)
library(DESeq2)
library(BiocParallel)

## Create integrated peak dataset ----------------------
message(paste(Sys.time(), ">> Create integrated peak dataset"))
peaks <- list.files("Peaks",
                    pattern="_peaks.broadPeak",
                    full.names=TRUE)
peaks <- peaks[!grepl("[[:digit:]]", peaks) &
                 grepl(opt$sample, peaks)]

peaks <- lapply(peaks, regioneR::toGRanges)
peaks <- GRangesList(peaks)
peaks <- unlist(peaks)
peaks <- regioneR::joinRegions(peaks)

saf <- data.frame(peaks)[,1:3]
colnames(saf) <- c("Chr", "Start", "End")
saf$GeneID <- paste0("peak_", 1:nrow(saf))
saf$Strand <- "+"
save(saf, file=paste0("diffAnalisis/", opt$experiment, "_", opt$sample, "_saf.rda"))

## Get matrix of counts --------------------------------
message(paste(Sys.time(), ">> Get matrix of counts"))
# List bams
bam <- list.files("BAM",
                  pattern=".bam$",
                  full.names=TRUE)
bam <- bam[!grepl("raw", bam)]
bam <- bam[grepl(opt$sample, bam)]

cov <- Rsubread::featureCounts(bam,
                               annot.ext=saf,
                               nthreads=10)

names <- pipelineNGS::getNameFromPath(bam, suffix=".bam")
mat <- cov$counts
colnames(mat) <- names

## Differential analysis with DESeq2 -------------------
message(paste(Sys.time(), ">> Differential analysis with DESeq2"))

coldata <- data.frame(ids=names,
                      sample=unlist(lapply(strsplit(names, "_"),
                                           function(x) x[1])),
                      treatment=unlist(lapply(strsplit(names, "_"),
                                              function(x) x[2])))
coldata$batch <- gsub("hi", "", gsub("endoc", "", coldata$sample))
coldata$tissue <- gsub("[[:digit:]]*", "", coldata$sample)

if(batch) form <- formula(~ batch + treatment) else form <- formula(~ treatment)

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = coldata,
                              design= form)

register(MulticoreParam(opt$cores))

dds <- DESeq(dds,
             parallel=TRUE)
save(dds,
     file=paste0("diffAnalisis/", opt$experiment, "_", opt$sample,
                 "_dds", b, ".rda"))

## Perform rlog transformation --------------------------
message(paste(Sys.time(), ">> Perform rlog transformation"))

rld <- rlog(dds)
save(rld,
     file=paste0("diffAnalisis/", opt$experiment, "_", opt$sample,
                 "_rld", b, ".rda"))

## Get normalized counts --------------------------------
message(paste(Sys.time(), ">> Get normalized counts"))

counts <- counts(dds, normalized=TRUE)
save(counts, file=paste0("diffAnalisis/", opt$experiment, "_",
                         opt$sample, "_normCounts", b, ".rda"))

## Annotate and create GRanges object -------------------
message(paste(Sys.time(), ">> Annotate and create GRanges object"))

load(paste0("diffAnalisis/", opt$experiment, "_", opt$sample, "_saf.rda"))
load(paste0("diffAnalisis/", opt$experiment, "_", opt$sample, "_dds", b, ".rda"))

res <- results(dds)
res$GeneID <- rownames(res)
res <- dplyr::left_join(data.frame(res), saf[,-5])
res.gr <- regioneR::toGRanges(res[,c(8:10,7,1:6)])

# Set type
res.gr$type <- "stable"
res.gr$type[res.gr$padj<opt$padj &
              res.gr$log2FoldChange>opt$foldChange] <- "gained"
res.gr$type[res.gr$padj<opt$padj &
              res.gr$log2FoldChange<(-opt$foldChange)] <- "lost"

# Get distance to TSS
## Get coding genes list
library(biomaRt)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position",
                            "ensembl_gene_id", "external_gene_name"),
               filters="biotype", values="protein_coding",
               mart=grch37)

genes <- regioneR::toGRanges(genes)
genes <- genes[seqnames(genes) %in% c(1:22, "X", "Y"),]
seqlevelsStyle(genes) <- "UCSC"

## Annotate ATAC-seq regions to nearest gene
anno <- data.frame(ChIPseeker::annotatePeak(res.gr, TxDb=genes))
anno <- anno[,c(6,21,19:20)]
mcols(res.gr) <- dplyr::left_join(as.data.frame(mcols(res.gr)), anno)
res.gr$annotation <- "Promoter"
res.gr$annotation[abs(res.gr$distanceToTSS)>opt$promoter] <- "Distal"

save(res.gr, file=paste0("diffAnalisis/", opt$experiment, "_", opt$sample, "_fc", opt$foldChange,
                         "_padj", opt$padj, "_GRanges", b, ".rda"))


