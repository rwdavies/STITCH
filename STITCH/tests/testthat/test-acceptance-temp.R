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
    seed = 3,
    chr = "chr5",
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster
)
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = chr,
    phasemaster = phasemaster
)



test_that("STITCH can initialize with reference data with three sizes of K vs number of haps", {

    ## these are diploid
    n_samples_per_pop <- 4
    LL <- c(2, 4, 6:8)
    refpackL <- make_reference_package(
        n_snps = length(LL),
        L = LL,
        n_samples_per_pop = n_samples_per_pop,
        reference_populations = c("CEU"),
        chr = 1,
        phasemaster = phasemaster[LL, ]
    )

    output_format <- "bgvcf"
    i_scenario <- 3

    for(i_scenario in 1:3) {

        outputdir <- make_unique_tempdir()
        ## scenarios are: too few, exactly right amount, too many
        K <- list(n_samples_per_pop * 4, n_samples_per_pop * 2, n_samples_per_pop)[[i_scenario]]
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackL$reference_haplotype_file,
            reference_legend_file = refpackL$reference_legend_file,
            K = K,
            S = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )
        ## AM HERE
        ## PAST REFERENCE FOR THE FIRST TIME (YAY!)
        ## NOW WORK ON PROPER DIPLOID ETC (hopefully straightforward given all prep work)

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }
})

