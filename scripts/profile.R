#!/usr/bin/env Rscript

library("proftools")
library("STITCH")

profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)
outputdir <- Sys.getenv("OUTPUTDIR")
nCores <- as.integer(Sys.getenv("N_CORES"))
if (Sys.getenv("USE") == "CRAMS") {
    bamlist <- ""
    cramlist = "cramlist.txt"
    reference = "mm10_2016_10_02.fa"
} else {
    bamlist <- "bamlist.txt"
    cramlist <- ""
    reference <- ""
}

if (Sys.getenv("OPTION") != "") {
    option <- paste0(", ", Sys.getenv("OPTION"))
} else {
    option <- ""
}

if (Sys.getenv("K") != "") {
    K <- as.integer(Sys.getenv("K"))
} else {
    K <- 4
}

chr <- "chr19"
posfile <- "pos.txt"
genfile <- "gen.txt"
tempdir <- paste0(outputdir, "/")
##print("FIX ME")##tempdir()     'genfile = genfile,',    
nCores <- as.integer(nCores)
command <- paste0(
    'STITCH(',
    'chr = chr,',
    'bamlist = bamlist,',
    'posfile = posfile,',
    'outputdir = outputdir,',
    'K = ', K, ',',
    'nGen = 100,',
    'nCores = nCores,',
    'tempdir = tempdir',
    option
)
if (cramlist != "") {
    command <- paste0(command, ', cramlist = cramlist', ', reference = reference')
}
command <- paste0(command, ')')

print(command)
eval(parse(text = command))

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
