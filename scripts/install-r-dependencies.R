#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("Rsamtools"))) {
    try(source("https://bioconductor.org/biocLite.R"))
    try(source("http://bioconductor.org/biocLite.R"))
    biocLite(ask = FALSE) # update packages
    biocLite("Rsamtools")
}
if (!suppressPackageStartupMessages(require("Rcpp")))
    install.packages("Rcpp", repos="http://cran.rstudio.com/")
if (!suppressPackageStartupMessages(require("RcppArmadillo")))
    install.packages("RcppArmadillo", repos="http://cran.rstudio.com/")
if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")
if (!suppressPackageStartupMessages(require("devtools")))
    install.packages("devtools", repos="http://cran.rstudio.com/")
if (!suppressPackageStartupMessages(require("testthat")))
    install.packages("testthat", repos="http://cran.rstudio.com/")
