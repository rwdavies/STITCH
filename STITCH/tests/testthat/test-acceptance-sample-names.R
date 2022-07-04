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



test_that("can change sample names and have them match genfile", {

    sampleNames_file <- tempfile()
    sampleNames <- data_package$sample_names
    sampleNames[2] <- "newname"
    write.table(
        sampleNames,
        file = sampleNames_file,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    
    genfile <- tempfile()
    system(paste0("sed 's/samp2/newname/' ", shQuote(data_package$genfile), " > ", genfile))
    
    outputdir <- make_unique_tempdir()    
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = genfile,
        outputdir = outputdir,
        sampleNames_file = sampleNames_file,
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


test_that("can change sample names using regenerateInput and have them match genfile", {

    sampleNames_file <- tempfile()
    sampleNames <- data_package$sample_names
    sampleNames[2] <- "newname2"
    write.table(
        sampleNames,
        file = sampleNames_file,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    
    outputdir <- make_unique_tempdir()    
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = outputdir,
        K = 2,
        S = 2,
        nGen = 100,
        nCores = 1,
        generateInputOnly = TRUE
    )

    genfile <- tempfile()
    system(paste0("sed 's/samp2/newname2/' ", shQuote(data_package$genfile), " > ", genfile))

    ## now need to modify manually
    load(file.path(outputdir, "RData", paste0("sampleNames.", data_package$chr, ".RData")))
    sampleNames[2] <- "newname2"
    save(sampleNames, file = file.path(outputdir, "RData", paste0("sampleNames.", data_package$chr, ".RData")))
    
    STITCH(
        chr = data_package$chr,
        posfile = data_package$posfile,
        genfile = genfile,
        outputdir = outputdir,
        K = 2,
        S = 2,
        nGen = 100,
        nCores = 1,
        regenerateInput = FALSE,
        regenerateInputWithDefaultValues = TRUE,
        originalRegionName = data_package$chr
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
