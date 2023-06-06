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

    af_cutoff <- 0.01 ## here, 1%. note maf < 1% corresponding to af > 99% is very rare so could ignore for simplicity?
    ref_error <- 0.001
    vcffile <- "/data/smew1/rdavies/ukbb_gel_2023_01_26/ukbb.20.20000000.20100000.vcf.gz"
    start <- 20000000


    ## check computational complexity
    ## check here
    for(iii in 1:3) {
        if (iii == 1) {af_cutoff <- 0.01} ## make ALL of them common
        if (iii == 2) {af_cutoff <- 0} ## make ALL of them common
        if (iii == 3) {af_cutoff <- 1} ## make ALL of them rare
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
            w = width,
               n_c = sum(outRcpp$snp_is_common),
            n_r = sum(!outRcpp$snp_is_common),
            n_s = outRcpp[["n_skipped"]],
            sapply(outRcpp, object.size)
            )
    })
    mat <- t(sapply(out, I))

        mat <- cbind(mat, mat[, 1] / mat[, 2])
        mat <- cbind(mat, mat[, 1] / mat[, 3])
        print(mat)

    }
    ## seems OK at first glance
    ## 10 minutes for 100,000 at scale

    ## OK, seems to work
    ## no problems yet with push back stuff, would eventually kick in, maybe?
    ## might need some more testing at some point
    
    
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
## deliberately make a few "SNPs" fail
L <- 1:n_snps
refs <- rep("A", n_snps)
alts <- rep("G", n_snps)
## make not a SNP
refs[500] <- "AG"
alts[510] <- "AG"
## make not bi-allelic
alts[520] <- "A,G"
## make same position
L[530] <- L[529]
L[540] <- L[539]
L[541] <- L[540]


refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = n_big_haps,
    reference_populations = c("GBR"),
    chr = chr,
    phasemaster = phasemaster2,
    L = L,
    refs = refs,
    alts = alts
)



## TODO
##
## ?be able to use only SNPs in desired SNP list
## ?be able to skip unwanted SNPs on the fly



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
    expect_equal(outR[["n_skipped"]], outRcpp[["n_skipped"]])
    expect_equivalent(outR[["ref_alleleCount"]], outRcpp[["ref_alleleCount"]]) ## not sure, doesnt' seem to matter


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
    
## Rcpp::cppFunction('
## Rcpp::NumericMatrix test() {
## int nsnps = 3;
## int nhaps = 10;
## Rcpp::NumericVector af(3);
## af(0) = 0.3;
## af(1) = 0.2;
## af(2) = 0.1;
##   Rcpp::NumericMatrix ref_alleleCount(nsnps, 3);
## double nhapsd = double(nhaps);
##   for(int is = 0; is < nsnps; is++) {
##     ref_alleleCount(is, 0) = af(is);
##   ref_alleleCount(is, 1) = nhapsd;
##   }
##   ref_alleleCount.column(2) = ref_alleleCount.column(0)  / ref_alleleCount.column(1);
## return(ref_alleleCount);
## }')
##     test()

## })


