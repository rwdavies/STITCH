test_that("can calculate likelihoods of flips in C++", {

    n_snps <- 5
    T <- n_snps
    K <- 2
    eHapsCurrent_t <- rbind(
        rep(0.99, n_snps),
        rep(0.01, n_snps)
    )
    alphaMatCurrent_t <- rbind(
        rep(0.5, n_snps),
        rep(0.5, n_snps)
    )
    priorCurrent <- rep(1 / K, K)
    sigmaCurrent <- rep(0.99, n_snps)
    transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)

    sampleReads <- list(
        list(
            2, 1,
            matrix(rep(-20, 3), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        ),
        list(
            2, 2,
            matrix(rep(20, 3), ncol = 1),
            matrix(c(1, 2, 3), ncol = 1)
        ),
        list(
            0, 3,
            matrix(-20, ncol = 1),
            matrix(3, ncol = 1)
        ),
        list(
            0, 3,
            matrix(20, ncol = 1),
            matrix(3, ncol = 1)
        )
    )
    reads_at_SNPs <- get_reads_at_SNP(sampleReads, T)

    swap_mat <- rbind(
        rep(1, 4),
        c(0, 1, 0, 1),
        rep(0, 4)
    )

    if (length(sampleReads) != ncol(swap_mat))
        stop("bad assumptions")

    eMatHap_t <- rcpp_make_eMatHap_t(
        sampleReads = sampleReads,
        nReads = length(sampleReads),
        eHaps_t = eHapsCurrent_t,
        maxDifferenceBetweenReads = 1000,
        Jmax = 10
    )

    ## 1 - how to sample
    ## 2 - where reads are
    out <- rcpp_calculate_many_likelihoods(
        swap_mat = swap_mat,
        reads_at_SNPs = reads_at_SNPs,
        eMatHap_t = eMatHap_t,
        sampleReads = sampleReads,
        nReads = length(sampleReads),
        eHaps_t = eHapsCurrent_t,
        maxDifferenceBetweenReads = 1000,
        Jmax = 10,
        pi = priorCurrent,
        transMatRate_t = transMatRate_t,
        alphaMat_t = alphaMatCurrent_t
    )

    ## require 2nd option is best
    cFinal <- rowSums(out)
    expect_equal(which.min(cFinal), 2)


})
