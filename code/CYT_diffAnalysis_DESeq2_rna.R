###################################
## Differential analysis RNA-seq ##
###################################

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
sample <- opt$sample # tissue to analyze
batch <- opt$batch # correct for batch
cores <- opt$cores # Number of cores to use for analysis
fc <- opt$foldChange # log2FC
padj <- opt$padj # FDR adjusted P

if (batch) b <- "Batch" else b <- ""

## Load libraries --------------------------------------
library(GenomicRanges)
library(DESeq2)
library(BiocParallel)

## Create data.frame with htseq-count info -------------
message(paste(Sys.time(), ">> Create integrated peak dataset"))
files <- list.files("htseq-count",
                    pattern=".txt",
                    full.names=FALSE)

names <- pipelineNGS::getNameFromPath(files,
                                      suffix=".htseq-count.txt")
condition <- unlist(lapply(strsplit(names, "_"),
                            function(x) x[2]))
batch.n <- unlist(lapply(strsplit(names, "_"),
                         function(x) x[1]))

sampleTable <- data.frame(sampleName = names,
                          fileName = files,
                          treatment = condition,
                          batch = batch.n)

## Differential analysis with DESeq2 -------------------
message(paste(Sys.time(), ">> Differential analysis with DESeq2"))

if(batch) form <- formula(~ batch + treatment) else form <- formula(~ treatment)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "htseq-count/",
                                  design = form)
register(MulticoreParam(cores))

dds <- DESeq(dds,
             parallel=TRUE)
save(dds,
     file=paste0("diffAnalisis/RNA_", sample,
                 "_dds", b, ".rda"))

## Perform rlog transformation --------------------------
message(paste(Sys.time(), ">> Perform rlog transformation"))

rld <- rlog(dds)
save(rld,
     file=paste0("diffAnalisis/RNA_", sample,
                 "_rld", b, ".rda"))

## Annotate and create GRanges object -------------------
message(paste(Sys.time(), ">> Annotate and create GRanges object"))

load(paste0("diffAnalisis/RNA_", sample, "_dds", b, ".rda"))

res <- results(dds)
res$GeneID <- rownames(res)

# Set type
res$type <- "stable"
res$type[res$padj<0.05 &
              res$log2FoldChange>fc] <- "gained"
res$type[res$padj<0.05 &
              res$log2FoldChange<(-fc)] <- "lost"

# Add coordinates for genes
library(biomaRt)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position",
                            "strand", "ensembl_gene_id", "external_gene_name", "gene_biotype"),
               mart=grch37)

res$GeneID <- gsub("\\.[[:digit:]]*", "", res$GeneID)
res <- dplyr::left_join(data.frame(res), genes,
                        by=c(GeneID="ensembl_gene_id"))
res$strand[res$strand==-1] <- "-"
res$strand[res$strand==1] <- "+"

res.gr <- regioneR::toGRanges(res[!is.na(res$chromosome_name),c(9:12,7,13,14,1:6,8)])
strand(res.gr) <- res.gr$strand
mcols(res.gr) <- mcols(res.gr)[,-1]
seqlevelsStyle(res.gr) <- "UCSC"

save(res.gr, file=paste0("diffAnalisis/RNA_", sample, "_GRanges", b, ".rda"))

## Get normalized counts --------------------------------
message(paste(Sys.time(), ">> Get normalized counts"))

counts <- counts(dds, normalized=TRUE)
save(counts, file=paste0("diffAnalisis/RNA_",
                         sample, "_normCounts", b, ".rda"))
