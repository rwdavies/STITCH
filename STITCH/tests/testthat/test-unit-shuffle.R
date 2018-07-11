test_that("can determine switches", {

    nGrids <- 100    
    break_results <- cbind(
        left_break = c(5, 20, 50),
        right_break = c(10, 30, 60)
    )
    nbreaks <- nrow(break_results)
    
    K <- 4
    fromMat <- array(0, c(nbreaks, K, K))
    ## now make simple
    for(k in 1:K) {
        fromMat[, k, k] <- 1
    }
    switchOrder <- determine_switch_order(fromMat, nbreaks, K)

    expect_equal(
        switchOrder,
        matrix(rep(1:4, 3), nrow = 3, byrow = TRUE)
    )

})


test_that("can smooth rate", {

    set.seed(10)            
    shuffle_bin_radius <- 5000
    L <- sort(sample(100000, 100))
    gridWindowSize <- NA
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids

    ## now simulate a rate
    sigmaSum_unnormalized <- runif(nGrids - 1)
    sigma_rate <- -log(sigmaSum_unnormalized) / grid_distances ## by definition

    ## just check the same
    results1 <- make_smoothed_rate(
        sigmaSum_unnormalized,
        sigma_rate,
        L_grid,
        grid_distances,
        nGrids,
        shuffle_bin_radius
    )

    results2 <- rcpp_make_smoothed_rate(
        sigmaSum_unnormalized,
        sigma_rate,
        L_grid,
        grid_distances,
        nGrids,
        shuffle_bin_radius
    )

    expect_equal(as.numeric(results1), as.numeric(results2))
    
})


test_that("can define breaks on tiny region", {

    set.seed(10)        
    tempdir <- tempdir()
    regionName <- "blargh"
    L <- sort(sample(100000, 100))
    shuffle_bin_radius <- 2000
    ## make sure it tests going to the end
    L <- unique(sort(c(1, shuffle_bin_radius, L)))
    L <- unique(sort(c(L, 100000 - shuffle_bin_radius + 1, 100000)))
    gridWindowSize <- NA
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    nGen <- 100
    minRate <- 0.1
    maxRate <- 100

    for(i in 1:2) {
        
        if (i == 1) {
            sigmaSum_unnormalized <- runif(nGrids - 1)
        } else {
            sigmaSum_unnormalized <- rep(0.8, nGrids - 1)
        }

        define_and_save_breaks_to_consider(
            tempdir,
            regionName,
            sigmaSum_unnormalized,
            L_grid,
            grid_distances,
            nGrids,
            nGen,
            minRate,
            maxRate,
            shuffle_bin_radius = 2000,
            plot_shuffle_haplotype_attempts = FALSE
        )

        load(file_break_results(tempdir, regionName))
        ## have an output
        expect_equal(length(break_results) > 0, TRUE)
    }
        
})


test_that("can choose points to break if highest on the left", {

    smoothed_rate <- runif(100)
    smoothed_rate[1:10] <- NA
    smoothed_rate[11] <- 100
    nGrids <- 100
    
    out <- choose_points_to_break(
        smoothed_rate,
        nGrids
    )
    
    expect_equal(length(out) > 0, TRUE)    
    
})


test_that("can choose points to break if highest on the right", {

    smoothed_rate <- runif(100)
    smoothed_rate[90:100] <- NA
    smoothed_rate[89] <- 100
    nGrids <- 100
    
    out <- choose_points_to_break(
        smoothed_rate,
        nGrids
    )
    expect_equal(length(out) > 0, TRUE)

})
