if (1 == 0) {

    ## source the environment
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
test_that("simple diploid method can phase", {

    skip("not currently exposed")
    
    n_snps <- 50
    L <- 31:80
    reads_span_n_snps <- 6
    chr <- 1
    n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
    extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")
    phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemaster[] <- sample(c(0, 1), n_snps * 2, replace = TRUE)
    data_package <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 77,
        chr = "chr5",
        K = 2,
        L = L,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster
    )
    regionStart <- 50
    regionEnd <- 70
    
    for(useTempdirWhileWriting in c(FALSE, TRUE)) {

        outputdir <- make_unique_tempdir()

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = 5,
            K = 2,
            S = 1,
            nGen = 100,
            nCores = 1,
            outputSNPBlockSize = 5,
            do_phasing = TRUE,
            useTempdirWhileWriting = useTempdirWhileWriting
        )

        regionName <- paste(data_package$chr, regionStart, regionEnd, sep = ".")
        which_snps <- match(regionStart:regionEnd, L)
        
        check_output_against_phase(
            file = file.path(outputdir, paste0("stitch.", regionName, ".vcf.gz")),
            data_package = data_package,
            which_snps = which_snps,
            tol = 0.2 ,
            allowed_pse = 0
        )
        
    }

})


test_that("something related to phasing", {

    skip("not currently exposed")
    ##
    ## simulate some data
    ##
    set.seed(10)
    
    n_snps <- 100
    reads_span_n_snps <- 3
    n_reads <- 1 / (reads_span_n_snps / n_snps) ## want about 1X / sample
    n_samples <- 50
    
    ## here we simulate 3 major ancestral haplotypes, but manually build a fourth one with a phase switch in it, to make things interesting
    K <- 4 
    phasemaster <- matrix(sample(c(0, 1), n_snps * K, replace = TRUE), ncol = K)
    phasemaster[, 1] <- c(phasemaster[1:40, 2], phasemaster[41:100, 3])
    ##
    haps_to_sample_all <- matrix(sample(2:4, 2 * n_samples, replace = TRUE), ncol = 2)
    haps_to_sample_all[1, ] <- c(1, 4) ## this sample has the recombined haplotype
    
    data_package <- make_acceptance_test_data_package(
        n_samples = n_samples,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 77,
        chr = "chr5",
        K = 2,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster,
        haps_to_sample_all = haps_to_sample_all
    )

    ## run STITCH normally
    outputdir <- make_unique_tempdir()    
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = K - 1, ## here we want to do using 3 not 4
        nGen = 100,
        nCores = 1,
        do_phasing = TRUE,
        phasing_method = 3
    )

    check_output_against_phase(
        file = file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        data_package = data_package,
        tol = 2, ## meh
        allowed_pse = 0,
        min_info = 0.85
    )
    
})    
