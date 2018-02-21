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
print(getwd())

knitr::knit('benchmarks/summarize_benchmarking.Rmd', output = "benchmarks/summarize_benchmarking.md")
system("rsync -av benchmarks/summarize_benchmarking.md florence:~/")
