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

})


test_that("can define breaks on tiny region", {

    skip("wer")
    
    set.seed(10)        
    tempdir <- tempdir()
    regionName <- "blargh"
    L <- sort(sample(100000, 100))
    gridWindowSize <- NA
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    nGen <- 100
    minRate <- 0.1
    maxRate <- 100
    shuffle_bin_radius <- 2000

    ## just generate something 
    sigmaSum_unnormalized <- runif(nGrids)    

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
    
})
