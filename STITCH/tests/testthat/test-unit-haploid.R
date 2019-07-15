test_that("can calculate likelihoods of flips in C++", {

    n_snps <- 5
    T <- n_snps
    K <- 2
    S <- 3
    eHapsCurrent_tc <- array(runif(K * n_snps * S), c(K, n_snps, S))
    alphaMatCurrent_tc <- array(runif(K * (n_snps - 1) * S), c(K, n_snps - 1, S))    
    for(s in 1:S) {
        eHapsCurrent_tc[, , s] <- rbind(
            rep(0.99, n_snps),
            rep(0.01, n_snps)
        )
        alphaMatCurrent_tc[, , s] <- rbind(
            rep(0.5, n_snps - 1),
            rep(0.5, n_snps - 1)
        )
    }
    priorCurrent <- array(1 / K, c(K, S))
    sigmaCurrent_m <- array(0.99, c(n_snps, S))
    transMatRate_t <- get_transMatRate_m(method = "diploid-inbred", sigmaCurrent_m = sigmaCurrent_m)

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
    nReads <- length(sampleReads)

    swap_mat <- rbind(
        rep(1, 4),
        c(0, 1, 0, 1),
        rep(0, 4)
    )

    if (length(sampleReads) != ncol(swap_mat)) {
        stop("bad assumptions")
    }

    ## right - did I rename this
    eMatRead_t <- array(1, c(K, nReads))
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = eHapsCurrent_tc,
        maxDifferenceBetweenReads = 1000,
        Jmax = 10,
	eMatHapOri_t = array(0, c(1, 1)), ## argh?
	pRgivenH1 = array(NA, c(1, 1)),
	pRgivenH2 = array(NA, c(1, 1)),
	run_pseudo_haploid = FALSE,
        s = 1,
        prev = 0,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text"
    )

    ## 1 - how to sample
    ## 2 - where reads are
    ## out <- rcpp_calculate_many_likelihoods(
    ##     swap_mat = swap_mat,
    ##     reads_at_SNPs = reads_at_SNPs,
    ##     eMatHap_t = eMatHap_t,
    ##     sampleReads = sampleReads,
    ##     nReads = length(sampleReads),
    ##     eHaps_t = eHapsCurrent_t,
    ##     maxDifferenceBetweenReads = 1000,
    ##     Jmax = 10,
    ##     pi = priorCurrent,
    ##     transMatRate_t = transMatRate_t,
    ##     alphaMat_t = alphaMatCurrent_t
    ## )

    ## ## require 2nd option is best
    ## cFinal <- rowSums(out)
    ## expect_equal(which.min(cFinal), 2)


})
