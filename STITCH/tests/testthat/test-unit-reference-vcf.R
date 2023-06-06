if ( 1 == 0 ) {

    library("testthat")
    library("STITCH")
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)

}





test_that("speed check", {

    skip("not needed often")

    af_cutoff <- 0.01 ## here, 1%. note maf < 1% corresponding to af > 99% is very rare so could ignore for simplicity?
    ref_error <- 0.001
    vcffile <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/ukbb.20.20000000.20100000.vcf.gz"
    start <- 20000000

    ## check computational complexity
    ## check here
    widths <- round(10 ** seq(2, 4, length.out = 5))
    out <- mclapply(widths, mc.cores = 10, function(width) {
        print(paste0("width = ", width, ", ", date()))
        region <- paste0("20:", start, "-", start + width)
        c(
            system.time({
            outRcpp <- Rcpp_get_hap_info_from_vcf(
                vcffile = vcffile,
                af_cutoff = af_cutoff,
                region = region
            )
            })[["elapsed"]],
               sum(outRcpp$snp_is_common),
            sum(!outRcpp$snp_is_common),
            sapply(outRcpp, object.size)
            )
    })
    mat <- t(sapply(out, I))

    mat
    cbind(mat, mat[, 1] / mat[, 2])
    cbind(mat, mat[, 1] / mat[, 3])
    ## seems OK at first glance
    ## 10 minutes for 100,000 at scale

    ## OK, seems to work
    ## no problems yet with push back stuff, would eventually kick in, maybe?
    ## might need some more testing at some point

    width <- 10000
    region <- paste0("20:", start, "-", start + width)
    outRcpp <- Rcpp_get_hap_info_from_vcf(
        vcffile = vcffile,
        af_cutoff = af_cutoff,
        region = region
    )
    
    
})



set.seed(919)

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


reads_span_n_snps <- 3
## want about 4X here
n_reads <- round(4 * n_snps / reads_span_n_snps)
data_package <- make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 3,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = K,
    phasemaster = phasemaster2
)
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = n_big_haps,
    reference_populations = c("GBR"),
    chr = chr,
    phasemaster = phasemaster2
)







## mostly for QUILT
test_that("can load and split into rare and common", {

    ## to test
    ## normal
    ## samples / populations
    ## regions
    ## samples, populations, regions

    af_cutoff <- 0.01 ## here, 1%. note maf < 1% corresponding to af > 99% is very rare so could ignore for simplicity?
    ref_error <- 0.001
    region <- paste0(chr, ":", 100, "-", 1100)
    vcffile <- refpack$reference_vcf_file
    
    outR <- R_get_hap_info_from_vcf(
        vcffile = vcffile,
        af_cutoff = af_cutoff,
        region = region
    )

    outRcpp <- Rcpp_get_hap_info_from_vcf(
        vcffile = vcffile,
        af_cutoff = af_cutoff,
        region = region
    )

    expect_equal(outR[["rare_per_hap_info"]], outRcpp[["rare_per_hap_info"]])
    expect_equal(outR[["pos"]], outRcpp[["pos"]])
    expect_equal(outR[["rhb_t"]], outRcpp[["rhb_t"]])
    expect_equal(outR[["snp_is_common"]], outRcpp[["snp_is_common"]])


})


## speed check?

## 

## test_that("scratch", {


## Rcpp::cppFunction('
## Rcpp::DataFrame test() {
## IntegerVector L(2);
## CharacterVector ref(2);
## CharacterVector alt(2);
## L(0)=1;
## L(1)=2;
## ref(0)="A";
## ref(1)="A";
## alt(0)="T";
## alt(1)="G";
##   Rcpp::DataFrame posX = Rcpp::DataFrame::create(Rcpp::Named("V1") = clone(L), Rcpp::Named("V2") = clone(L));
## return(posX);
## }')
##     test()
    
     
## })
