test_that("can determine switches", {

    nGrids <- 100    
    break_results <- cbind(
        left_grid_break_0_based = c(5, 20, 50),
        left_grid_focal_0_based = c(10, 24, 58),
        right_grid_focal_0_based = c(11, 25, 59),
        right_grid_break_0_based = c(15, 30, 60)
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



test_that("can run getBetterSwitchesSimple with / without grid", {

    K <- 4
    nSNPs <- 100
    set.seed(2342)
    L <- sort(sample(10000, 100))

    for(gridWindowSize in c(NA, 500)) {

        out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids
        snps_in_grid_1_based <- out$snps_in_grid_1_based
        
        eHapsFuture_t <- matrix(runif(K * nSNPs), nrow = K, ncol = nSNPs)
        alphaMatFuture_t <- matrix(runif(K * (nGrids - 1)), nrow = K, ncol = nGrids - 1)
        
        ## I think these are grids but are they 1 or 0 based
        break_results <- cbind(
            left_grid_break_0_based = c(0, nGrids - 4),
            left_grid_focal_0_based = c(1, nGrids - 3),
            right_grid_focal_0_based = c(2, nGrids - 2),
            right_grid_break_0_based = c(3, nGrids - 1)
        )
        if (nSNPs == nGrids) {
            break_results <- rbind(break_results, c(40, 43, 44, 50))
        }
        break_results <- break_results[order(break_results[, 2]), ]        
        nbreaks <- nrow(break_results)

        ## from -> to
        fromMat <- array(0, c(nbreaks, K, K))
        ## make wacky!
        for(k in 1:K) {
            fromMat[, k, c(K, 1:(K - 1))[k]] <- 1
        }

        out <- getBetterSwitchesSimple(
            fromMat = fromMat,
            nbreaks = nbreaks,
            break_results = break_results,
            K = K,
            eHapsFuture_t = eHapsFuture_t,
            alphaMatFuture_t = alphaMatFuture_t,
            grid = grid,
            snps_in_grid_1_based = snps_in_grid_1_based
        )
        test_that(ncol(out$eHapsCurrent_t), nSNPs)
        test_that(ncol(out$alphaMatCurrent_t), nGrids - 1)

    }

})


test_that("can smooth rate", {

    set.seed(10)            
    shuffle_bin_radius <- 5000
    L <- sort(sample(100000, 100))

    for(gridWindowSize in c(NA, 10000)) {
        
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

    }
    
})


test_that("can define breaks on tiny region", {

    set.seed(10)        
    tempdir <- tempdir()
    regionName <- "blargh"
    L <- sort(sample(1000000, 500))
    shuffle_bin_radius <- 2000
    ## make sure it tests going to the end
    L <- unique(sort(c(1, shuffle_bin_radius, L)))
    L <- unique(sort(c(L, 100000 - shuffle_bin_radius + 1, 100000)))
    nGen <- 100
    minRate <- 0.1
    maxRate <- 100

    for(gridWindowSize in c(NA, 10000)) {

        out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids

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
                plot_shuffle_haplotype_attempts = TRUE
            )

            load(file_break_results(tempdir, regionName))
            ## have an output
            expect_equal(length(break_results) > 0, TRUE)
            
        }

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

