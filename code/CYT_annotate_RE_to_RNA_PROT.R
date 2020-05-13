################################
## Correlation RNA-seq and RE ##
################################

## Create options for script ---------------------------------------
library("optparse")

option_list = list(
  make_option(c("-f", "--foldChange"), type="numeric", default=1,
              help="log2FC of the RE dataset to use", metavar="numeric"),
  make_option(c("-q", "--padj"), type="numeric", default=1,
              help="Adjusted p-value of the RE dataset to use", metavar="numeric"),
  make_option(c("-w", "--winWidth"), type="numeric", default=4e4,
              help="Size of the window from TSS used for annotating REs to genes [default= %default]",
              metavar="numeric"),
  make_option(c("-s", "--sample"), type="character", default="endoc",
              help="Sample to use for this analysis. [default= %default]",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Load necessary data ---------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))

load(paste0("../../data/RNA/diffAnalisis/RNA_", opt$sample, "_GRangesBatch.rda")) # RNA annotation
load(paste0("../../data/REs/REs_", opt$sample, "_fc", opt$foldChange,
            "_padj", opt$padj, "_granges.rda")) # RE annotation
load(paste0("../../data/RNA/diffAnalisis/RNA_", opt$sample, "_normCountsBatch.rda")) # RNA counts
load("../../data/Proteomics/proteomics_data_type.rda") # Proteomics data

## Create windows and annotate REs ---------------------------------
message(">> Annotating windows")
rna <- res.gr[res.gr$gene_biotype=="protein_coding",] # Select protin coding
win <- promoters(rna) # Reduce to promoters
win <- resize(rna, width=opt$winWidth, fix="center") # Resize to window

ols <- findOverlaps(re, win) # find overlaps RE and window
anno <- unique(cbind(mcols(re)[queryHits(ols), c(1,16)],
                     mcols(win)[subjectHits(ols), c(1,2,5)]))

## Annotate with RNA-seq counts ------------------------------------
message(">> Adding RNA-seq counts")

# Converto to data.frame
counts <- data.frame(counts)
counts$GeneID <- rownames(counts)
counts$GeneID <- gsub("\\.[[:digit:]]*", "", counts$GeneID)

# Reshape to long format
counts <- reshape2::melt(counts,
                         id.vars=11,
                         value.vars=1:10,
                         value.name="counts",
                         variable.name="sample")
counts$treatment <- unlist(lapply(strsplit(as.character(counts$sample), "_"),
                                  function(x) x[2]))

# Calculate mean counts
counts.mean <- counts %>%
  group_by(GeneID, treatment) %>%
  summarise(counts=mean(counts))

# Merge with annotation
anno.counts <- left_join(data.frame(anno), counts.mean)
anno.counts <- anno.counts[!is.na(anno.counts$type),]

## Annotate with # of IREs ----------------------------------------
message(">> Grouping genes by # of IREs")

t <- data.frame(table(anno$GeneID, anno$type))
t <- t[t$Var2=="IRE",]

t$group <- "0"
t$group[t$Freq==1] <- "1"
t$group[t$Freq==2] <- "2"
t$group[t$Freq>=3] <- ">=3"
colnames(t)[1] <- "GeneID"

anno.group <- left_join(unique(data.frame(anno)[,c(3:5)]), t)
anno.group$group <- factor(anno.group$group,
                           levels=c("0", "1", "2", ">=3"))

## Save RNA data --------------------------------------------------
save(anno, anno.counts, anno.group, file=paste0("data/RNA_annotation_", opt$sample,
                                                "_fc", opt$foldChange,
                                                "_padj", opt$padj, "_win", opt$winWidth, ".rda"))

## Annotat to proteins --------------------------------------------
message(">> Annotating to proteins")

prot <- data[,c(9,4,17)]
colnames(prot)[2:3] <- paste0("prot.", colnames(prot)[2:3])

anno.prot <- dplyr::left_join(data.frame(anno),
                              prot)

anno.prot <- anno.prot[!is.na(anno.prot$prot.log2FoldChange),]
anno.prot <- anno.prot[!is.na(anno.prot$type),]
anno.prot$type <- factor(anno.prot$type, levels=c("SRE", "IRE"))

save(anno.prot, file=paste0("data/PROT_annotation_", opt$sample,
                            "_fc", opt$foldChange,
                            "_padj", opt$padj, "_win", opt$winWidth, ".rda"))
