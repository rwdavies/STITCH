if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

}


## make bigger
n_snps <- 1200
K <- 6
n_big_haps <- 100 ## 1000
chr <- 10

## make a spectrum or rare and common SNPs
phasemaster1 <- array(sample(c(0, 1), n_snps * K, replace = TRUE), c(n_snps, K))
phasemaster2 <- phasemaster1[, sample(1:K, n_big_haps, replace = TRUE)]

## make about 10% common, and do not touch
n_common <- round(n_snps / 10)
common <- sample(1:n_snps, n_common)
## make about 90% rare, and make about 1% freq
phasemaster2[-common, ] <- 0
## make remaining frequency about 1%
phasemaster2[-common, ] <- sample(c(0, 1), (n_snps - n_common) * n_big_haps, prob = c(0.99, 0.01), replace = TRUE)


reads_span_n_snps <- 50
## heck it is fine
n_reads <- 100
## deliberately make a few "SNPs" fail
L <- 1:n_snps
x <- replicate(n_snps, sample(c("A", "C", "G", "T"), 2))
refs <- x[1, ]
alts <- x[2, ]


data_package_crams <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 3,
    chr = chr,
    K = 2,
    L = L,
    refs = refs,
    alts = alts,
    phasemaster = phasemaster2,
    reads_span_n_snps = reads_span_n_snps,
    use_crams = TRUE
)

test_that("STITCH diploid works under default parameters using CRAM files", {

    
    outputdir <- make_unique_tempdir()    

    set.seed(843)
    
    ## alread moved
    ref <- data_package_crams$ref

    ## check an error thrown when trying to access WITHOUT ref
    expect_error(
        system(paste0("samtools view -h ", shQuote(data_package_crams$cram_files[1])),
               fixed = TRUE)
    )
    expect_equal(    system(paste0("samtools view -T ", shQuote(ref), " -h ", shQuote(data_package_crams$cram_files[1]), " > /dev/null")), 0L)

    print(paste0("The ref being passed in is:", ref))
    ## run
    STITCH(
        chr = data_package_crams$chr,
        cramlist = data_package_crams$cramlist,
        reference = ref,
        posfile = data_package_crams$posfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1
    )
    
    check_output_against_phase(
        file.path(outputdir, paste0("stitch.", data_package_crams$chr, ".vcf.gz")),
        data_package_crams,
        output_format,
        which_snps = NULL,
        tol = 0.2
    )

})
