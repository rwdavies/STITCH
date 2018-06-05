#!/usr/bin/env Rscript

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
if (!suppressPackageStartupMessages(require("roxygen2")))
    install.packages("roxygen2", repos="http://cran.rstudio.com/")
if (!suppressPackageStartupMessages(require("data.table")))
    install.packages("data.table", repos="http://cran.rstudio.com/")
#if (!suppressPackageStartupMessages(require("rrbgen")))
#    install.packages("https://github.com/rwdavies/rrbgen/raw/master/releases/rrbgen_0.0.0.tar.gz", repos=NULL)

