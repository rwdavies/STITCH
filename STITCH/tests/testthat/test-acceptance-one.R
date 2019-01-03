## one single test, useful for basic acceptance tests
test_that("simple diploid method can work", {

    n_snps <- 10
    reads_span_n_snps <- 6
    chr <- 1
    n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
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
        phasemaster = phasemaster
    )
    
    outputdir <- make_unique_tempdir()    
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
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
