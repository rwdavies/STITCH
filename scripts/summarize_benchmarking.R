#!/usr/bin/env Rscript

## Script to summarize benchmarking runs

library("knitr")
STITCH_home <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/"
setwd(STITCH_home)

knitr::knit('benchmarks/summarize_benchmarking.Rmd', output = "benchmarks/summarize_benchmarking.md")
system("rsync -av benchmarks/summarize_benchmarking.md florence:~/")
