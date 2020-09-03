########################
## FOOTPRINT ANALYSIS ##
########################

library(maRge)
out_dir <- "output/cyt_03_subclassification"

## Make tag directories ------------------------------------
dir.create("data/CYT/ATAC/HOMER_tags", F)
bam.files <- list.files("data/CYT/ATAC/BAM",
                        pattern="offset.bam$", full.names=TRUE)
bam.files <- bam.files[grep("endoc", bam.files)]

makeTagDir("data/CYT/ATAC/HOMER_tags/endoc_ctrl",
           bam.files=bam.files[grep("ctrl", bam.files)])

makeTagDir("/data/CYT/ATAC/HOMER_tags/endoc_cyt",
           bam.files=bam.files[grep("cyt", bam.files)])

## Get motif coverage ---------------------------------------
tag.dirs <- list.dirs("data/CYT/ATAC/HOMER_tags",
                      recursive=F)

cov <-
  coverageFootprintHOMER(regions.file="data/bedfiles/IREs_endoc_fc1_padj0.05_Neo_distal.bed",
                         motif.file=file.path(out_dir, "HOMER_IREs_endoc_fc1_padj0.05_opening_distal_mask/homerResults/motif1.motif"),
                         out.name=file.path(out_dir,"FOOT_ISRE_neo_distal"),
                         tags.dirs=tag.dirs)

save(cov, file=file.path(out_dir, "FOOT_ISRE_neo_distal_cov.rda"))

cov <-
  coverageFootprintHOMER(regions.file="data/bedfiles/IREs_endoc_fc1_padj0.05_Primed_distal.bed",
                       motif.file=file.path(out_dir, "HOMER_IREs_endoc_fc1_padj0.05_primed_distal_mask/homerResults/motif1.motif"),
                       out.name=file.path(out_dir,"FOOT_ISRE_primed_distal"),
                       tags.dirs=tag.dirs)

save(cov, file=file.path(out_dir, "FOOT_ISRE_primed_distal_cov.rda"))
