#!/usr/bin/env Rscript

## Script to summarize benchmarking runs

library("knitr")
## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
STITCH_home <- getwd()
setwd(STITCH_home)

setwd("benchmarks")
knitr::knit('summarize_benchmarking.Rmd', output = "summarize_benchmarking.md")
setwd("../")
system("rsync -av benchmarks florence:~/")
