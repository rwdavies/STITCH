#!/usr/bin/env Rscript

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))

library("testthat")
library("parallel")
source("STITCH/R/cli.R")
source("STITCH/R/test-drivers.R")

## testthat doesn't do what I want outside of package form
## so don't bother wrappping, just fail

cli_function_build <- Sys.getenv("CLI_FUNCTION_BUILD")
if (cli_function_build != "") {
    print(paste0("Using ", cli_function_build))
    dir <- tempdir()
    system(paste0("cp ", cli_function_build, " ", dir, "/"))
    system(paste0("(cd ", dir, " && tar -zxvf ", dir, "/*tar.gz STITCH/R/functions.R)"))
    function_file <- file.path(dir, "STITCH/R/functions.R")
} else {
    function_file <- "STITCH/R/functions.R"
}

system(paste0("cp ", function_file, " ~/TEMP.R"))
## make CLI file
cli_output_file <- "STITCH.R"
make_STITCH_cli(
    function_file = function_file,
    cli_output_file = cli_output_file
)
system(paste0("chmod +x ", cli_output_file))



message("test that STITCH CLI produces help message")
## this is bad - optparse exits with code 1 on desired behaviour
## so can't check help message printed on error code alone
## do some super minimal parsing of the output
out <- suppressWarnings(system(
    paste0(cli_output_file, " --help "), intern = TRUE
))
expect_equal(grep("Options", out) > 0, TRUE)
expect_equal(attr(out, "status"), 1)




n_snps <- 5
chr <- 10
phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = 4,
    seed = 1,
    chr = chr,
    K = 2,
    phasemaster = phasemaster
)
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = chr
)



## this also test character, integer, double and NA variables
message("test that STITCH CLI can work")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        "--K=2",
        "--nGen=100",
        "--buffer=NA",
        "--nCores=1"
    ),
    stdout = stdout_file, stderr = stderr_file
)
stderr <- system(paste0("cat ", stderr_file), intern = TRUE)
stdout <- system(paste0("cat ", stdout_file), intern = TRUE)
expect_equal(0, out)


message("test that STITCH CLI stops when bad variable given")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        "--K=2",
        "--nGen=100AAAAAAAAAA"
    ),
    stdout = stdout_file, stderr = stderr_file
)
expect_equal(out > 0, TRUE)


cli_version <- Sys.getenv("CLI_VERSION")
error_check <- 2
if (cli_version != "") {
    ##message(paste0("Using CLI_VERSION=", cli_version))
    version <- as.numeric(strsplit(cli_version, ".", fixed = TRUE)[[1]])
    if ((version[2] >= 3 & version[3] >=7 ) | version[2] >= 4) {
        error_check <- 2
    } else {
        error_check <- 1
    }
}
##message(paste0("Using error_check = ", error_check))

message("test that STITCH CLI parses a logical variable correctly")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        "--outputInputInVCFFormat=TRUE",
        "--K=2",
        "--nGen=100"
    ),
    stdout = stdout_file, stderr = stderr_file
)
expect_equal(0, out)
## check this occured
stderr <- system(paste0("cat ", stderr_file), intern = TRUE)
stdout <- system(paste0("cat ", stdout_file), intern = TRUE)

if (error_check == 1) {
    expect_equal(length(grep("Build vcf from input", stdout)), 1)
    expect_equal(length(grep("teration", stdout)), 0)
} else if (error_check == 2) {
    expect_equal(length(grep("Build VCF from input", stderr)), 1)
    expect_equal(length(grep("teration", stderr)), 0)
} else {
    stop("bad CLI test")
}


message("test that STITCH CLI parses integer vector refillIterations correctly")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        "--K=2",
        "--nGen=100",
        "--refillIterations='c(4,30)'"
    ),
    stdout = stdout_file, stderr = stderr_file
)
expect_equal(0, out)
if (error_check == 1) {
    out_log <- system(paste0("cat ", stdout_file), intern = TRUE)
} else {
    out_log <- system(paste0("cat ", stderr_file), intern = TRUE)
}
a <- out_log[grep("refill infrequently used haplotypes", out_log)]
expect_equal(sum(sapply(sapply(c(4, 30), function(x) grep(x, a)), length) == 0), 0)



message("test that STITCH CLI uses default integer vector correctly")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        "--K=2",
        "--nGen=100"
    ),
    stdout = stdout_file, stderr = stderr_file
)
expect_equal(0, out)
## check this occured
if (error_check == 1) {
    out_log <- system(paste0("cat ", stdout_file), intern = TRUE)
} else {
    out_log <- system(paste0("cat ", stderr_file), intern = TRUE)
}
a <- out_log[grep("refill infrequently used haplotypes", out_log)]
## manually know default is 6, 10, 14, 18
expect_equal(sum(sapply(sapply(c(6, 10, 14, 18), function(x) grep(x, a)), length) == 0), 0)



message("test that STITCH CLI parses character vector correctly")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        paste0("--reference_haplotype_file=", refpack$reference_haplotype_file),
        paste0("--reference_legend_file=", refpack$reference_legend_file),
        paste0("--reference_sample_file=", refpack$reference_sample_file),
        "--reference_populations='c(\"CEU\",\"GBR\")'",
        "--K=2",
        "--nGen=100"
    ),
    stdout = stdout_file, stderr = stderr_file
)
expect_equal(0, out)
## check this occured
if (error_check == 1) {
    out_log <- system(paste0("cat ", stdout_file), intern = TRUE)
} else {
    out_log <- system(paste0("cat ", stderr_file), intern = TRUE)
}
check <- sapply(c("CEU", "GBR"), function(pop) {
    length(grep(pop, out_log))
})
expect_equal(sum(check == 0), 0)
