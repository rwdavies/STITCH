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
    transMatRate_t_H <- test_package$transMatRate_t_H
    alphaMatCurrent_t <- test_package$alphaMatCurrent_t
    eMatHapSNP_t <- test_package$eMatHapSNP_t
    priorCurrent <- test_package$priorCurrent

    ## 
    ## c++ - pass by reference
    ##
    alphaHat_t <- array(0, c(K, nGrids))
    c <- array(0, nGrids)
    Rcpp_run_forward_haploid(
        alphaHat_t = alphaHat_t,
        c = c,
        eMatHapSNP_t = eMatHapSNP_t,
        alphaMat_t = alphaMatCurrent_t,
        transMatRate_t_H = transMatRate_t_H,
        T = nGrids,
        K = K,
        pi = priorCurrent,
        alphaStart = array(0, 100),
        run_fb_subset = FALSE
    )

    ## 
    ## R (returns alphaHat_t, c)
    ## 
    out_R <- R_run_forward_haploid(
        alphaHat_t = array(0, c(K, nGrids)),
        c = array(0, nGrids),
        eMatHapSNP_t = eMatHapSNP_t,
        alphaMat_t = alphaMatCurrent_t,
        transMatRate_t_H = transMatRate_t_H,
        T = nGrids,
        K = K,
        pi = priorCurrent,
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
        eMatHapSNP_t = eMatHapSNP_t,
        alphaMat_t = alphaMatCurrent_t,    
        transMatRate_t_H = transMatRate_t_H
    )
    
    out_R2 <- R_run_backward_haploid(
        betaHat_t = array(0, c(K, nGrids)),
        c = out_R$c,
        eMatHapSNP_t = eMatHapSNP_t,
        alphaMat_t = alphaMatCurrent_t,    
        transMatRate_t_H = transMatRate_t_H
    )
    
    expect_equal(betaHat_t, out_R2$betaHat_t)

})
