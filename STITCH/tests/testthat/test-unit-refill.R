test_that("can refill simple", {

    nSNPs <- 20
    K <- 10
    N <- 1000
    L <- 1:nSNPs
    out <- assign_positions_to_grid(L, gridWindowSize = NA)
    grid <- out$grid
    nGrids <- out$nGrids
    L_grid <- out$L_grid
    grid_distances <- out$grid_distances
    distance_between_check <- 3
    
    hapSum <- array(1, c(nGrids, K))
    gammaSum <- array(0, c(nSNPs, K))
    for(k in 1:K) 
        gammaSum[, k] <- (k - 0.5)/ K

    ## make one bad hapSum
    hapSum[, 1] <- 0
    hapSum <- hapSum / rowSums(hapSum) * N    
    
    gammaSum_new <- refillSimple(
        hapSum_t = t(hapSum),
        nGrids = nGrids,
        K = K,
        gammaSum_t = t(gammaSum),
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,
        grid = grid
    )$gammaSum

    expect_equal(sum(gammaSum_new[, 1] == gammaSum[, 1]), 0)

})


test_that("can refill simple", {

    K <- 10
    N <- 1000
    L <- c(5, 8, 9, 13, 14, 18, 19, 21, 24)
    nSNPs <- length(L)    
    out <- assign_positions_to_grid(L, gridWindowSize = NA)
    grid <- out$grid
    nGrids <- out$nGrids
    L_grid <- out$L_grid
    grid_distances <- out$grid_distances
    distance_between_check <- 10
    
    hapSum <- array(1, c(nGrids, K))
    gammaSum <- array(0, c(nSNPs, K))
    for(k in 1:K) 
        gammaSum[, k] <- (k - 0.5)/ K

    hapSum[1:3, 1] <- 0
    hapSum[8:9, 1] <- 0    
    hapSum <- hapSum / rowSums(hapSum) * N    
    
    gammaSum_new_t <- refillSimple(
        hapSum_t = t(hapSum),
        nGrids = nGrids,
        K = K,
        gammaSum_t = t(gammaSum),
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,        
        grid = grid
    )$gammaSum_t
    gammaSum_new <- t(gammaSum_new_t)

    expect_equal(sum(gammaSum_new[1:3, 1] == gammaSum[1:3, 1]), 0)
    expect_equal(sum(gammaSum_new[4:7, 1] != gammaSum[4:7, 1]), 0)        
    expect_equal(sum(gammaSum_new[8:9, 1] == gammaSum[8:9, 1]), 0)    

})


test_that("can refill simple with grid", {

    nSNPs <- 20
    K <- 10
    N <- 1000
    L <- 1:nSNPs
    out <- assign_positions_to_grid(L, gridWindowSize = 5)
    grid <- out$grid
    nGrids <- out$nGrids
    L_grid <- out$L_grid
    grid_distances <- out$grid_distances
    distance_between_check <- 3
    
    hapSum <- array(1, c(nGrids, K))
    gammaSum <- array(0, c(nSNPs, K))
    for(k in 1:K) 
        gammaSum[, k] <- (k - 0.5)/ K

    ## make one bad hapSum
    hapSum[, 1] <- 0
    hapSum <- hapSum / rowSums(hapSum) * N    
    
    gammaSum_new <- refillSimple(
        hapSum_t = t(hapSum),
        nGrids = nGrids,
        K = K,
        gammaSum_t = t(gammaSum),
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,                
        grid = grid
    )$gammaSum

    expect_equal(sum(gammaSum_new[, 1] == gammaSum[, 1]), 0)

})



test_that("can refill simple with grid for only part of hapSum", {

    nSNPs <- 20
    K <- 10
    N <- 1000
    L <- 1:nSNPs
    gridWindowSize <- 5
    out <- assign_positions_to_grid(L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    nGrids <- out$nGrids
    L_grid <- out$L_grid
    grid_distances <- out$grid_distances
    distance_between_check <- 5
    
    hapSum <- array(1, c(nGrids, K))
    gammaSum <- array(0, c(nSNPs, K))
    for(k in 1:K) 
        gammaSum[, k] <- (k - 0.5)/ K

    ## make two bad blocks
    hapSum[1, 1] <- 0
    hapSum[4, 1] <- 0    
    hapSum <- hapSum / rowSums(hapSum) * N    
    
    gammaSum_new_t <- refillSimple(
        hapSum_t = t(hapSum),
        nGrids = nGrids,
        K = K,
        gammaSum_t = t(gammaSum),
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,                        
        grid = grid
    )$gammaSum_t
    gammaSum_new <- t(gammaSum_new_t)

    expect_equal(sum(gammaSum_new[grid == 0, 1] == gammaSum[grid == 0, 1]), 0)
    expect_equal(sum(gammaSum_new[grid == 1, 1] != gammaSum[grid == 1, 1]), 0)
    expect_equal(sum(gammaSum_new[grid == 2, 1] != gammaSum[grid == 2, 1]), 0)        
    expect_equal(sum(gammaSum_new[grid == 3, 1] == gammaSum[grid == 3, 1]), 0)    

})
