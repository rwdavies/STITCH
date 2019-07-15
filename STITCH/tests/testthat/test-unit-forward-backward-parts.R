## test compability between development R and production C++ code
## ensure has features expected of HMM e.g. certain sums

test_that("R and c++ run forward haploid algorithms are the same", {

    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 10,
        gridWindowSize = 3
    )
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    eMatGrid_t <- test_package$list_of_eMatGrid_t[[1]]
    priorCurrent_m <- test_package$priorCurrent_m

    ##
    ## c++ - pass by reference
    ##
    alphaHat_t <- array(0, c(K, nGrids))
    c <- array(0, nGrids)
    Rcpp_run_forward_haploid(
        alphaHat_t = alphaHat_t,
        c = c,
        eMatGrid_t = eMatGrid_t,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        priorCurrent_m = priorCurrent_m,
        s = 0, ## 0-based
        alphaStart = array(0, 1),
        run_fb_subset = FALSE
    )

    ##
    ## R (returns alphaHat_t, c)
    ##
    out_R <- R_run_forward_haploid(
        alphaHat_t = array(0, c(K, nGrids)),
        c = array(0, nGrids),
        eMatGrid_t = eMatGrid_t,
        alphaMat_t = alphaMatCurrent_tc[, , 1],
        transMatRate_t_H = transMatRate_tc_H[, , 1],
        T = nGrids,
        K = K,
        pi = priorCurrent_m[, 1],
        alphaStart = 0,
        run_fb_subset = FALSE
    )

    expect_equal(alphaHat_t, out_R$alphaHat_t)
    expect_equal(c, out_R$c)

    ##
    ## second part, test backward algorithms
    ##

    ##
    ## R (returns alphaHat_t, c)
    ##
    betaHat_t <- array(0, c(K, nGrids))
    Rcpp_run_backward_haploid(
        betaHat_t = betaHat_t,
        c = c,
        eMatGrid_t = eMatGrid_t,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        s = 0
    )

    out_R2 <- R_run_backward_haploid(
        betaHat_t = array(0, c(K, nGrids)),
        c = out_R$c,
        eMatGrid_t = eMatGrid_t,
        alphaMat_t = alphaMatCurrent_tc[, , 1],
        transMatRate_t_H = transMatRate_tc_H[, , 1]
    )

    expect_equal(betaHat_t, out_R2$betaHat_t)

})
