test_that("can refill simple", {

    for(gridWindowSize in c(NA, 3)) {
        test_package <- make_fb_test_package(
            K = 4,
            nReads = 8,
            nSNPs = 20,
            gridWindowSize = gridWindowSize,
            S = 2
        )
        S <- test_package$S
        sampleReads <- test_package$sampleReads
        nSNPs <- test_package$nSNPs
        nGrids <- test_package$nGrids
        K <- test_package$K
        transMatRate_tc_H <- test_package$transMatRate_tc_H
        transMatRate_tc_D <- test_package$transMatRate_tc_D    
        alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
        eMatGrid_t <- test_package$list_of_eMatGrid_t[[1]]
        priorCurrent_m <- test_package$priorCurrent_m
        eHapsCurrent_tc <- test_package$eHapsCurrent_tc
        grid <- test_package$grid
        N <- test_package$N
        L_grid <- test_package$L_grid
        
        distance_between_check <- 3 ## how often to check
        
        hapSum_tc <- array(1, c(K, nGrids, S))
        gammaSum_tc <- array(0, c(K, nSNPs, S))
        for(s in 1:S) {
            for(k in 1:K) {
                gammaSum_tc[k,, s] <- (k - 0.5)/ K
            }
            hapSum_tc[1, , s] <- 0
            hapSum_tc[, , s] <- hapSum_tc[, , s] / colSums(hapSum_tc[, , s]) * N            
        }
        
        ## make one bad hapSum
        gammaSum_tc_new <- refillSimple(
            hapSum_tc = hapSum_tc,
            gammaSum_tc = gammaSum_tc,
            N = N,
            distance_between_check = distance_between_check,
            L_grid = L_grid,
            grid = grid
        )$gammaSum_tc
        
        expect_equal(sum(gammaSum_tc_new[1, , 1] == gammaSum_tc[1, , 1]), 0)

    }

})


test_that("can refill simple again", {

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
    S <- 2
    
    hapSum_tc <- array(1, c(K, nGrids, S))
    gammaSum_tc <- array(0, c(K, nSNPs, S))
    for(s in 1:S) {
        for(k in 1:K) {
            gammaSum_tc[k, , s] <- (k - 0.5)/ K
        }
        hapSum_tc[1, 1:3, s] <- 0
        hapSum_tc[1, 8:9, s] <- 0    
        hapSum_tc[, , s] <- hapSum_tc[, , s] / rowSums(hapSum_tc[, , s]) * N
    }
    
    gammaSum_tc_new <- refillSimple(
        hapSum_tc = hapSum_tc,
        gammaSum_tc = gammaSum_tc,
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,        
        grid = grid
    )$gammaSum_tc

    expect_equal(sum(gammaSum_tc_new[1, 1:3, 1] == gammaSum_tc[1, 1:3, 1]), 0)
    expect_equal(sum(gammaSum_tc_new[1, 4:7, 1] != gammaSum_tc[1, 4:7, 1]), 0)        
    expect_equal(sum(gammaSum_tc_new[1, 8:9, 1] == gammaSum_tc[1, 8:9, 1]), 0)    

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
    S <- 3
    
    hapSum_tc <- array(1, c(K, nGrids, S))
    gammaSum_tc <- array(0, c(K, nSNPs, S))
    for(s in 1:S) {
        for(k in 1:K) {
            gammaSum_tc[k, , s] <- (k - 0.5)/ K
            ## make two bad blocks
            hapSum_tc[1, 1, 1] <- 0
            hapSum_tc[1, 4, 1] <- 0    
            hapSum_tc[, , s] <- hapSum_tc[, , s]/ rowSums(hapSum_tc[, , s]) * N    
        }
    }

    gammaSum_tc_new <- refillSimple(
        hapSum_tc = hapSum_tc,
        gammaSum_tc = gammaSum_tc,
        N = N,
        distance_between_check = distance_between_check,
        L_grid = L_grid,                        
        grid = grid
    )$gammaSum_tc

    expect_equal(sum(gammaSum_tc_new[1, grid == 0, 1] == gammaSum_tc[1, grid == 0, 1]), 0)
    expect_equal(sum(gammaSum_tc_new[1, grid == 1, 1] != gammaSum_tc[1, grid == 1, 1]), 0)
    expect_equal(sum(gammaSum_tc_new[1, grid == 2, 1] != gammaSum_tc[1, grid == 2, 1]), 0)        
    expect_equal(sum(gammaSum_tc_new[1, grid == 3, 1] == gammaSum_tc[1, grid == 3, 1]), 0)    

})
