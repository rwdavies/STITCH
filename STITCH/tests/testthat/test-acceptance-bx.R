if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    ##    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

}


## one single test, useful for basic acceptance tests
test_that("can use bx tag", {

    n_snps <- 50
    reads_span_n_snps <- 6
    chr <- 1
    n_reads <- round(5 / (reads_span_n_snps / n_snps)) ## want about 5X / sample
    extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")
    phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemaster[2, ] <- c(1, 0)
    phasemaster[3, ] <- c(0, 0)
    phasemaster[7, ] <- c(1, 0)    
    data_package <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 77,
        chr = "chr5",
        K = 2,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster,
        use_bx_tag = TRUE
    )
    
    outputdir <- make_unique_tempdir()    
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
        S = 2,
        nGen = 100,
        nCores = 1
    )

    vcf <- read.table(
        file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})
