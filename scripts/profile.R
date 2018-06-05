#!/usr/bin/env Rscript

library("proftools")
library("STITCH")

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
        stitch_dir <- getwd()
    }
}

profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)

profile_start <- Sys.time()
##################
outputdir <- Sys.getenv("OUTPUTDIR")
if (outputdir == "") {
    outputdir <- file.path(getwd(), script_dir, "../", "test-results", "profile-one-off")
    dir.create(outputdir, showWarnings = FALSE)
}
nCores <- Sys.getenv("N_CORES")
if (nCores == "") {
    nCores <- 1
}
if (Sys.getenv("USE") == "CRAMS") {
    bamlist <- ""
    cramlist = "cramlist.txt"
    reference = "mm10_2016_10_02.fa"
} else {
    bamlist <- "bamlist.txt"
    cramlist <- ""
    reference <- ""
}

if (Sys.getenv("OPTION") != "" & Sys.getenv("OPTION") != "NA") {
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
nCores <- as.integer(nCores)
setwd("./test-data/mouse_data/")
command <- paste0(
    'STITCH(',
    'chr = chr,',
    'bamlist = bamlist,',
    'posfile = posfile,',
    'genfile = genfile,',    
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
##################
profile_end <- Sys.time()

Rprof(NULL)
pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- Sys.getenv("OUTPUT_PLOT")
if (output_plot == "") {
    setwd(stitch_dir)
    output_plot <- file.path("profile.pdf")
}
pdf(output_plot, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
title(title, outer=TRUE)
dev.off()

print(profile_end - profile_start)
