dir.create("data/", F)
library(RnBeads)

###############################################
## Set up
##---------------------------------------------
## Directories
data.dir <- "idat/"
idat.dir <- data.dir
sample.annotation <- file.path("./", "SampleSheet.csv")
analysis.dir <- "data_analysis"
report.dir <- "reports"

## Options for the analysis
rnb.options(filtering.sex.chromosomes.removal=TRUE,
            identifiers.column="Sample_ID",
            export.to.csv=TRUE)
rnb.initialize.reports(report.dir)
logger.start(fname=NA) ## Restrict log messages to console
num.cores <- 10
parallel.setup(num.cores) ## Enable parallel programming


###############################################
## Load data
##---------------------------------------------
data.source <- c(idat.dir, sample.annotation)
result <- rnb.run.import(data.source=data.source,
                         data.type="infinium.idat.dir",
                         dir.reports=report.dir)

rnb.set <- result$rnb.set
rnb.set

save.rnb.set(rnb.set, path="data_analysis/rnb.set.raw", archive=F)

###############################################
## Quality Control
##---------------------------------------------
rnb.run.qc(rnb.set, report.dir)
save.rnb.set(rnb.set, path="data_analysis/rnb.set.qc", archive=F)


###############################################
## Filtering and normalization
##---------------------------------------------
rnb.set.unfiltered <- rnb.set
result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports=report.dir)
rnb.set <- result$rnb.set
save.rnb.set(rnb.set, path="data_analysis/rnb.set.norm", archive=F)


###############################################
## Exploratory Analysis
##---------------------------------------------
rnb.run.exploratory(rnb.set, report.dir)
save.rnb.set(rnb.set, path="data_analysis/rnb.set.expl", archive=F)


###############################################
## Differential Methylation
##---------------------------------------------
rnb.run.differential(rnb.set, report.dir)
save.rnb.set(rnb.set, path="data_analysis/rnb.set.dm", archive=F)

# cmp.cols <- "Sample_Group"
# rnb.execute.computeDiffMeth(rnb.set, cmp.cols)


## Export data
rnb.run.tnt(rnb.set, report.dir)
save.rnb.set(rnb.set, path="data_analysis/rnb.set.fin", archive=F)

