test_that("can do single reference iteration", {

    ## hmm, need to dummy up a bit
    ## should be able to do with N=1
    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 10,
        gridWindowSize = 3,
        S = 2
    )
    S <- test_package$S
    sampleReads <- test_package$sampleReads
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    eMatHapSNP_t <- test_package$list_of_eMatHapSNP_t[[1]]
    sigmaCurrent_m <- test_package$sigmaCurrent_m
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    grid <- test_package$grid
    grid_distances <- test_package$grid_distances

    ## dummy up a bit
    tempdir <- tempfile("tempdir")
    dir.create(tempdir)
    N_haps <- 2
    nCores <- 1
    reference_bundling_info <- NULL
    nGen <- 10000

    regionName <- "jimmy"
    for(iBam in 1:2) {
        save(sampleReads, file = file_referenceSampleReads(tempdir, iBam, regionName))
    }

    out <- single_reference_iteration(eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, sigmaCurrent_m = sigmaCurrent_m, priorCurrent_m = priorCurrent_m, N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, L = L, grid = grid, grid_distances = grid_distances, nGen = nGen)


})



test_that("speed test", {

    skip("WER")
    test_package <- make_fb_test_package(
        K = 100,
        nReads = 100,
        nSNPs = 100000,
        gridWindowSize = NA
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
    Rcpp_run_forward_haploid_normal(
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

})

test_that("speed test", {

    skip("WER")
    ## so seems the cube approach is best
    ## hopefully not that complicated
    ## eHapsCurrent_tc <- suffix "tc" for transposed c?
    ## alphaMatCurrent_tc
    ## hapSumCurrent_tc
    ## sigmaCurrent_t ## i imagine... for consistency?
    K <- 100
    nSNPs <- 5000
    S <- 10
    cube_eHaps_t <- array(runif(K * nSNPs * S), c(K, nSNPs, S))
    list_of_eHaps_t <- lapply(1:S, function(i) {
        return(cube_eHaps_t[, , i])
    })
    gamma_t <- array(runif(K * nSNPs), c(K, nSNPs))
    eHaps_input <- array(runif(K * nSNPs), c(K, nSNPs))

    f <- function(option) {
        test_eHaps_options(
            cube_eHaps_t,
            list_of_eHaps_t,
            gamma_t,
            eHaps_input,
            option,
            nSNPs,
            K,
            S
        )
    }

    expect_equal(f("best_possible"), f("direct_slice"))
    library("microbenchmark")
    m <- microbenchmark(
        f("direct"), f("direct_slice"), f("list"), f("pre_get_slice"), f("best_possible"),
        times = 10
    )
    print(m)


})
