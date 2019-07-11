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


n_snps <- 10
reads_span_n_snps <- 6
chr <- 1
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")

phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
phasemaster[2, ] <- c(1, 0)
phasemaster[3, ] <- c(0, 0)
phasemaster[4, ] <- c(1, 1)
phasemaster[7, ] <- c(1, 0)    
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 3,
    chr = "chr5",
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster
)


test_that("STITCH can generate interim plots", {

    for(method in get_available_methods()) {

        outputdir <- make_unique_tempdir()
        chr <- data_package$chr
        niterations <- 5
        STITCH(
            chr = chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            method = method,
            nCores = 1,
            niterations = niterations,
            plotHapSumDuringIterations = TRUE,
            shuffleHaplotypeIterations = NA,
	    refillIterations = NA
        )

        ## only care about plumbing and making plots
        d <- dir(file.path(outputdir, "plots"))
        expect_equal(length(grep(paste0("hapSum.", chr, ".iteration."), d)), niterations - 1)
        expect_equal(length(grep(paste0("alphaMat.", chr, ".iteration."), d)), 2 * (niterations - 1))

    }

})
