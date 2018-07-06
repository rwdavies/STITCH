if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/Google Drive/STITCH/"
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

n_samples <- 20
n_snps <- 10
K <- 4
reads_span_n_snps <- 6
chr <- 1
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")

## make 4 haplotyes
set.seed(22313)
phasemaster <- matrix(sample(c(0, 1), K * n_snps, replace = TRUE), ncol = K)

data_package_inbred <- make_acceptance_test_data_package(
    n_samples = n_samples,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 3,
    chr = "chr5",
    K = K,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster,
    samples_are_inbred = TRUE
)

test_that("STITCH diploid-inbred works under default parameters", {

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")

        outputdir <- make_unique_tempdir()                
        STITCH(
            chr = data_package_inbred$chr,
            bamlist = data_package_inbred$bamlist,
            posfile = data_package_inbred$posfile,
            genfile = data_package_inbred$genfile,
            outputdir = outputdir,
            method = "diploid-inbred",
            K = K,
            nGen = 100,
            nCores = 1,
            keepTempDir = TRUE,
            output_format = output_format
        )

        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package_inbred$chr, extension[output_format])),
            data_package_inbred,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )
        
    }

})


