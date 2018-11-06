#!/usr/bin/env Rscript

required_packages <- c("proftools", "Rcpp", "RcppArmadillo", "optparse", "devtools", "testthat", "roxygen2", "data.table")
for(package in required_packages) {
    if (!suppressPackageStartupMessages(require(package, character.only = TRUE))) {
        out <- install.packages(package, repos="http://cran.rstudio.com/")
        out <- require(package)
        if (!out) {
            stop(paste0("Failed to install package:", package))
        }
    }
}
if (!suppressPackageStartupMessages(require("rrbgen")))
    install.packages("https://github.com/rwdavies/rrbgen/raw/master/releases/rrbgen_0.0.4.tar.gz", repos=NULL)

