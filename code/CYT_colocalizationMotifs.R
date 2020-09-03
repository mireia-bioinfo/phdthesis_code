##################################################
## Colocalization of motifs in primed enhancers ##
##################################################

## Get location of motifs in primed IREs ---------------------------------
## Join all significant motifs into single file
path <- catSignMotifsHOMER(path_output="../RE_03_deNovoMotif/HOMER_IREs_endoc_fc1_padj0.05_primed_distal_mask/", 
                   num_sign=9)

## Get location for motifs

getPositionMotifs <- function(bedfile="../../data/bedfiles/IREs_endoc_fc1_padj0.05_Primed_distal.bed",
                              name="primed",
                              motifs=path) {
  findMotifs <- paste("~/tools/homer/bin/findMotifsGenome.pl", 
                      bedfile,
                      "hg19 test/ -p 6 -size given -mask",
                      "-find", motifs,
                      ">", paste0("data/MOTIFS_", name, "_findGenome.tsv"))
  system(findMotifs)
  
  tsv <- read.delim(paste0("data/MOTIFS_", name, "_findGenome.tsv"),
                    stringsAsFactors=FALSE)
  
  tsv$motif <- unlist(lapply(strsplit(tsv$Motif.Name, ":"), 
                             function(x) x[[2]]))
  
  tsv$Best.Match <- unlist(lapply(lapply(strsplit(tsv$motif, "(", fixed=TRUE),
                                         function(x) x[-length(x)]),
                                  function(x) paste(x, collapse="(")))
  tsv$Best.Match[tsv$Best.Match==""] <- tsv$motif[tsv$Best.Match==""]
  
  tsv$Motif.Name <- unlist(lapply(strsplit(tsv$motif, "/", fixed=TRUE),
                                  function(x) x[1]))
  tsv$Rank <- 0
  cont <- 1
  for (i in unique(tsv$motif)) {
    tsv$Rank[tsv$motif==i] <- cont
    cont <- cont + 1
  }
  
  tsv$Motif.Name <- paste0(tsv$Rank, "-", tsv$Motif.Name)
  
  table(tsv$Motif.Name)
  
  ## Obtain coordinates
  load("../../data/REs/REs_endoc_fc1_padj0.05_granges_subgroup.rda")
  re <- data.frame(re)[,c(1:3,6,23)]
  colnames(re)[4:5] <- c("PositionID", "type")
  mot <- unique(dplyr::left_join(tsv, re))
  mot <- mot[!is.na(mot$start),]
  
  ## Create GRanges
  mot.gr <- list()
  for (i in unique(mot$Motif.Name)) {
    mot.gr[[i]] <- regioneR::toGRanges(mot[mot$Motif.Name==i, c(10:12,1,13,4)])
  }
  
  mot.gr <- GRangesList(mot.gr)
  save(mot.gr, file=paste0("data/MOTIFS_", name, "_regions_granges.rda"))
}

getPositionMotifs(bedfile="../../data/bedfiles/IREs_endoc_fc1_padj0.05_Primed_distal.bed",
                  name="primed",
                  motifs=path)

getPositionMotifs(bedfile="../../data/bedfiles/SREs_endoc_fc1_padj0.05_distal.bed",
                  name="SRE",
                  motifs=path)

## Create dataframe, keeping redundant ---------------------------------
load("data/MOTIFS_primed_regions_granges.rda")

mot <- unlist(mot.gr)

ols <- findOverlaps(mot, mot)

coloc <- cbind(data.frame(mot)[queryHits(ols),c(6:8)],
               data.frame(mot)[subjectHits(ols),c(6:8)])
colnames(coloc) <- c(paste0("query_", colnames(coloc)[1:3]),
                     paste0("subj_", colnames(coloc)[4:6]))

# coloc <- coloc[!(coloc$query_Motif.Name==coloc$subj_Motif.Name &
#                  coloc$subj_PositionID==coloc$subj_PositionID),]

mat <- table(coloc$query_Motif.Name, coloc$subj_Motif.Name)
attributes(mat)$class <- "matrix" 
# mat[lower.tri(mat, diag = FALSE)] <- NA

df <- reshape2::melt(mat)
# df <- df[!is.na(df$value),]
# df <- df[!(df$Var1==df$Var2),]

colnames(df)[3] <- "primed"

mot.coloc <- df

###

load("data/MOTIFS_SRE_regions_granges.rda")

mot <- unlist(mot.gr)

ols <- findOverlaps(mot, mot)

coloc <- cbind(data.frame(mot)[queryHits(ols),c(6:8)],
               data.frame(mot)[subjectHits(ols),c(6:8)])
colnames(coloc) <- c(paste0("query_", colnames(coloc)[1:3]),
                     paste0("subj_", colnames(coloc)[4:6]))

# coloc <- coloc[!(coloc$query_Motif.Name==coloc$subj_Motif.Name &
#                    coloc$subj_PositionID==coloc$subj_PositionID),]

mat <- table(coloc$query_Motif.Name, coloc$subj_Motif.Name)
attributes(mat)$class <- "matrix" 
# mat[lower.tri(mat, diag = FALSE)] <- NA

df <- reshape2::melt(mat)
# df <- df[!is.na(df$value),]
# df <- df[!(df$Var1==df$Var2),]

colnames(df)[3] <- "SREs"

mot.coloc <- dplyr::left_join(mot.coloc, df)

mot.coloc <- reshape2::melt(mot.coloc,
                            id.vars=c(1:2),
                            value.vars=3:4,
                            variable.name="type",
                            value.name="overlaps")

load("../../data/REs/REs_endoc_fc1_padj0.05_granges_subgroup.rda")
mot.coloc$total <- NA
mot.coloc$total[grep("primed", mot.coloc$type)] <- length(re$atac.GeneID[re$atac.annotation=="Distal" &
                                                  grepl("Primed", re$subgroup2)])
mot.coloc$total[grep("SRE", mot.coloc$type)] <- length(re$atac.GeneID[re$atac.annotation=="Distal" &
                                               grepl("SRE", re$subgroup2)])

mot.coloc$percentage <- mot.coloc$overlaps/mot.coloc$total*100

rename <- data.frame(bestMatch=as.character(unique(mot.coloc$Var1)),
                     TFname=c("ISRE", "HNF1", "FOX", "RFX", "MAF", "PDX1", 
                              "NKX6.1", "STAT", "NEUROD1"),
                     type=c("Inflammatory", "Islet", "Other", "Other", "Islet", "Islet",
                            "Islet", "Inflammatory", "Islet"),
                     stringsAsFactors = FALSE)

mot.coloc$Var1 <- as.character(mot.coloc$Var1)
mot.coloc$Var2 <- as.character(mot.coloc$Var2)

for (i in 1:nrow(rename)) {
  mot.coloc$Var1[grep(rename$bestMatch[i], as.character(mot.coloc$Var1), fixed=TRUE)] <- rename$TFname[i]
  mot.coloc$Var2[grep(rename$bestMatch[i], as.character(mot.coloc$Var2), fixed=TRUE)] <- rename$TFname[i]
  
}

save(mot.coloc, file="data/MOT-COLOC_df.rda")

## Create dataframe, removing redundant ---------------------------------
load("data/MOTIFS_primed_regions_granges.rda")

mot <- unlist(mot.gr)

ols <- findOverlaps(mot, mot)

coloc <- cbind(data.frame(mot)[queryHits(ols),c(6:8)],
               data.frame(mot)[subjectHits(ols),c(6:8)])
colnames(coloc) <- c(paste0("query_", colnames(coloc)[1:3]),
                     paste0("subj_", colnames(coloc)[4:6]))

coloc <- coloc[!(coloc$query_Motif.Name==coloc$subj_Motif.Name &
                 coloc$subj_PositionID==coloc$subj_PositionID),]

mat <- table(coloc$query_Motif.Name, coloc$subj_Motif.Name)
attributes(mat)$class <- "matrix" 
mat[lower.tri(mat, diag = FALSE)] <- NA

df <- reshape2::melt(mat)
df <- df[!is.na(df$value),]
df <- df[!(df$Var1==df$Var2),]

colnames(df)[3] <- "primed"

mot.coloc <- df

###

load("data/MOTIFS_SRE_regions_granges.rda")

mot <- unlist(mot.gr)

ols <- findOverlaps(mot, mot)

coloc <- cbind(data.frame(mot)[queryHits(ols),c(6:8)],
               data.frame(mot)[subjectHits(ols),c(6:8)])
colnames(coloc) <- c(paste0("query_", colnames(coloc)[1:3]),
                     paste0("subj_", colnames(coloc)[4:6]))

coloc <- coloc[!(coloc$query_Motif.Name==coloc$subj_Motif.Name &
                   coloc$subj_PositionID==coloc$subj_PositionID),]

mat <- table(coloc$query_Motif.Name, coloc$subj_Motif.Name)
attributes(mat)$class <- "matrix" 
mat[lower.tri(mat, diag = FALSE)] <- NA

df <- reshape2::melt(mat)
df <- df[!is.na(df$value),]
df <- df[!(df$Var1==df$Var2),]

colnames(df)[3] <- "SREs"

mot.coloc <- dplyr::left_join(mot.coloc, df)

mot.coloc <- reshape2::melt(mot.coloc,
                            id.vars=c(1:2),
                            value.vars=3:4,
                            variable.name="type",
                            value.name="overlaps")

load("../../data/REs/REs_endoc_fc1_padj0.05_granges_subgroup.rda")
mot.coloc$total <- NA
mot.coloc$total[grep("primed", mot.coloc$type)] <- length(re$atac.GeneID[re$atac.annotation=="Distal" &
                                                                           grepl("Primed", re$subgroup2)])
mot.coloc$total[grep("SRE", mot.coloc$type)] <- length(re$atac.GeneID[re$atac.annotation=="Distal" &
                                                                        grepl("SRE", re$subgroup2)])

mot.coloc$percentage <- mot.coloc$overlaps/mot.coloc$total*100

rename <- data.frame(bestMatch=c(as.character(unique(mot.coloc$Var1)), "9-NeuroD1(bHLH)"),
                     TFname=c("ISRE", "HNF1", "FOX", "RFX", "MAF", "PDX1", 
                              "NKX6.1", "STAT", "NEUROD1"),
                     type=c("Inflammatory", "Islet", "Other", "Other", "Islet", "Islet",
                            "Islet", "Inflammatory", "Islet"),
                     stringsAsFactors = FALSE)

mot.coloc$Var1 <- as.character(mot.coloc$Var1)
mot.coloc$Var2 <- as.character(mot.coloc$Var2)

for (i in 1:nrow(rename)) {
  mot.coloc$Var1[grep(rename$bestMatch[i], as.character(mot.coloc$Var1), fixed=TRUE)] <- rename$TFname[i]
  mot.coloc$Var2[grep(rename$bestMatch[i], as.character(mot.coloc$Var2), fixed=TRUE)] <- rename$TFname[i]
  
}

save(mot.coloc, file="data/MOT-COLOC_df_rmRedundant.rda")
