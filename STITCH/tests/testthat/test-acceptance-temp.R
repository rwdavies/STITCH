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



test_that("STITCH diploid works under default parameters with nCores = 40 and N = 25", {

    chr <- "chr3"
    n_snps <- 10
    phasemaster25 <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemaster25[3, ] <- c(1, 0)
    phasemaster25[7, ] <- c(0, 1)    
    phasemaster25[10, ] <- c(1, 1)    
    data_package25 <- make_acceptance_test_data_package(
        n_samples = 25,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = 2,
        phasemaster = phasemaster25
    )
    extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")
    
    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()
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
        STITCH(
            chr = data_package25$chr,
            bamlist = data_package25$bamlist,
            posfile = data_package25$posfile,
            genfile = data_package25$genfile,
            outputdir = outputdir,
            K = 4,
            nGen = 100,
            nCores = 40,
            output_format = output_format
        )

        data_package25$phase[, 1, ]        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package25$chr, extension[output_format])),
            data_package25,
            output_format,
            which_snps = NULL,
            tol = 0.2,
            min_info = 0.85
        )

    }

})
