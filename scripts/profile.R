#!/usr/bin/env Rscript

library("proftools")
library("STITCH")

profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)
outputdir <- Sys.getenv("OUTPUTDIR")
nCores <- Sys.getenv("N_CORES")
if (Sys.getenv("USE") == "CRAMS") {
    bamlist <- ""
    cramlist = "cramlist.txt"
    reference = "mm10_2016_10_02.fa"
} else {
    bamlist <- "bamlist.txt"
    cramlist <- ""
    reference <- ""
}
     
##    bamlist = "bamlist.txt",
STITCH(
    chr = "chr19",
    bamlist =  bamlist,
    cramlist = cramlist,
    reference = reference,
    posfile = "pos.txt",
    genfile = "gen.txt",
    outputdir = outputdir,
    K = 4,
    nGen = 100,
    nCores = as.integer(nCores),
    tempdir = paste0(tempdir(), "/")
)
Rprof(NULL)
pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- Sys.getenv("OUTPUT_PLOT")
pdf(output_plot, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
title(title, outer=TRUE)
dev.off()
