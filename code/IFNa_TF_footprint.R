#############################
## CALCULATE TF FOOTPRINTS ##
#############################

## Load packages ------------------------------
library(maRge)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)

## Load data ----------------------------------
out_dir <- "output/ifn_03_footprint"
load("data/IFNa/clusters_RNA/tf_clusters_ifna.rda")
clust.tf$cluster <- as.character(clust.tf$cluster)
clust.tf$TFs <- as.character(clust.tf$TFs)

matJaspar <- getMotifMatrixJASPAR(matrixtype="PFM")

##--------------------------------------------
## 2 hours
##--------------------------------------------
## Define bam files
bam.ctrl_2h <- list.files("data/IFNa/ATAC/BAMs", pattern="ctrl-2h_[[:digit:]].bam$",
                          full.names=TRUE)
bam.ifn_2h <- list.files("data/IFNa/ATAC/BAMs", pattern="ifn-2h_[[:digit:]].bam$",
                         full.names=TRUE)

footprints_2h.ctrl <- list()
footprints_2h.ifn <- list()
for (i in c2) {
  message(paste(">> Calculating Footprints for cluster", i))
  for (j in unique(clust.tf$TFs[clust.tf$cluster==i])) {
    mot <- matJaspar[grep(j, name(matJaspar))]

    if (length(mot) > 0) {
      mot <- mot[[1]]
      message(paste("----", j, ":", ID(mot)))

      mot <- Matrix(mot)
      mot <- prop.table(mot,2)

      ctrl <- suppressWarnings(factorFootprints(bam.ctrl_2h,
                                                pfm=mot,
                                                genome=Hsapiens,
                                                min.score="95%",
                                                regionSubset=cl.atac.2h.gr[cl.atac.2h.gr$cluster==i],
                                                upstream=100, downstream=100))

      ifn <- suppressWarnings(factorFootprints(bam.ifn_2h,
                                               pfm=mot,
                                               genome=Hsapiens,
                                               min.score="95%",
                                               regionSubset=cl.atac.2h.gr[cl.atac.2h.gr$cluster==i],
                                               upstream=100, downstream=100))
      footprints_2h.ctrl[[i]][[j]] <- ctrl
      footprints_2h.ifn[[i]][[j]] <- ifn
      save(ctrl, ifn, file=file.path(out_dir, paste0("footprint_", i, "_", j, "_2h.rda")))
    } else {
      message(paste("----", j, ": NOT FOUND"))
      footprints_2h.ctrl[[i]][[j]] <- NA
      footprints_2h.ifn[[i]][[j]] <- NA

      next
    }
  }
}

save(footprints_2h.ctrl, file=file.path(out_dir, "footprints_2h_ctrl.rda"))
save(footprints_2h.ifn, file=file.path(out_dir, "footprints_2h_ifn.rda"))

##--------------------------------------------
## 24 hours
##--------------------------------------------
## Define bam files
bam.ctrl_24h <- list.files("data/IFNa/ATAC/BAMs", pattern="ctrl-24h_[[:digit:]].bam$",
                           full.names=TRUE)
bam.ifn_24h <- list.files("data/IFNa/ATAC/BAMs", pattern="ifn-24h_[[:digit:]].bam$",
                          full.names=TRUE)

## Generate footprints
footprints_24h.ctrl <- list()
footprints_24h.ifn <- list()
for (i in c24) {
  message(paste(">> Calculating Footprints for cluster", i))
  for (j in unique(clust.tf$TFs[clust.tf$cluster==i])) {
    mot <- matJaspar[grep(j, name(matJaspar))]

    if (length(mot) > 0) {
      mot <- mot[[1]]
      message(paste("----", j, ":", ID(mot)))

      mot <- Matrix(mot)
      mot <- prop.table(mot,2)

      ctrl <- suppressWarnings(factorFootprints(bam.ctrl_24h,
                                                pfm=mot,
                                                genome=Hsapiens,
                                                min.score="95%",
                                                regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==i],
                                                upstream=100, downstream=100))

      ifn <- suppressWarnings(factorFootprints(bam.ifn_24h,
                                               pfm=mot,
                                               genome=Hsapiens,
                                               min.score="95%",
                                               regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==i],
                                               upstream=100, downstream=100))
      footprints_24h.ctrl[[i]][[j]] <- ctrl
      footprints_24h.ifn[[i]][[j]] <- ifn
      save(ctrl, ifn, file=file.path(out_dir, paste0("footprint_", i, "_", j, "_24h.rda")))
    } else {
      message(paste("----", j, ": NOT FOUND"))
      footprints_24h.ctrl[[i]][[j]] <- NA
      footprints_24h.ifn[[i]][[j]] <- NA

      next
    }
  }
}

save(footprints_24h.ctrl, file=file.path(out_dir, "footprints_24h_ctrl.rda"))
save(footprints_24h.ifn, file=file.path(out_dir, "footprints_24h_ifn.rda"))


##--------------------------------------------
## FOXA2
##--------------------------------------------
mot <- read.delim("fox2.motif", header=FALSE)
mot <- t(mot)
rownames(mot) <- c("A", "C", "G", "T")

mot[1,12] <- mot[1,12] + 0.001

cl <- clust.tf[grep("FOXA2", clust.tf$TFs),]

####### cluster 2
ctrl <- suppressWarnings(factorFootprints(bam.ctrl_2h,
                                          pfm=mot,
                                          genome=Hsapiens,
                                          min.score="95%",
                                          regionSubset=cl.atac.2h.gr[cl.atac.2h.gr$cluster==cl$cluster[1]],
                                          upstream=100, downstream=100))

ifn <- suppressWarnings(factorFootprints(bam.ifn_2h,
                                         pfm=mot,
                                         genome=Hsapiens,
                                         min.score="95%",
                                         regionSubset=cl.atac.2h.gr[cl.atac.2h.gr$cluster==cl$cluster[1]],
                                         upstream=100, downstream=100))
save(ctrl, ifn, file=file.path(out_dir, paste0("footprint_", cl$cluster[1], "_", "FOXA2", "_2h.rda")))

####### cluster 1B
ctrl <- suppressWarnings(factorFootprints(bam.ctrl_24h,
                                          pfm=mot,
                                          genome=Hsapiens,
                                          min.score="95%",
                                          regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==cl$cluster[2]],
                                          upstream=100, downstream=100))

ifn <- suppressWarnings(factorFootprints(bam.ifn_24h,
                                         pfm=mot,
                                         genome=Hsapiens,
                                         min.score="95%",
                                         regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==cl$cluster[2]],
                                         upstream=100, downstream=100))
save(ctrl, ifn, file=file.path(out_dir, paste0("footprint_", cl$cluster[2], "_", "FOXA2", "_24h.rda")))

####### cluster 2B
ctrl <- suppressWarnings(factorFootprints(bam.ctrl_24h,
                                          pfm=mot,
                                          genome=Hsapiens,
                                          min.score="95%",
                                          regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==cl$cluster[3]],
                                          upstream=100, downstream=100))

ifn <- suppressWarnings(factorFootprints(bam.ifn_24h,
                                         pfm=mot,
                                         genome=Hsapiens,
                                         min.score="95%",
                                         regionSubset=cl.atac.24h.gr[cl.atac.24h.gr$cluster==cl$cluster[3]],
                                         upstream=100, downstream=100))
save(ctrl, ifn, file=file.path(out_dir, paste0("footprint_", cl$cluster[3], "_", "FOXA2", "_24h.rda")))
