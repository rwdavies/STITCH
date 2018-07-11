test_that("can determine how assign SNPs to blocks for output on toy example", {
    
    grid <- 0:9 ## 0-based
    start_and_end_minus_buffer <- c(1, 10) ## 1-based
    outputSNPBlockSize <- 1000
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )

    expect_equivalent(blocks_for_output[, "snp_start_1_based"], 1)
    expect_equivalent(blocks_for_output[, "snp_end_1_based"], 10)    
    expect_equivalent(blocks_for_output[, "grid_start_0_based"], 0)
    expect_equivalent(blocks_for_output[, "grid_end_0_based"], 9)    
    
})

test_that("can determine how assign SNPs to blocks for output on tiny example", {

    ## yes this is supported
    grid <- c(0, 1)
    start_and_end_minus_buffer <- c(2, 2)
    outputSNPBlockSize <- 1000
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )

    expect_equivalent(blocks_for_output[, "snp_start_1_based"], 2)
    expect_equivalent(blocks_for_output[, "snp_end_1_based"], 2)    
    expect_equivalent(blocks_for_output[, "grid_start_0_based"], 1)
    expect_equivalent(blocks_for_output[, "grid_end_0_based"], 1)    
    
})



test_that("can determine how assign SNPs to blocks for output on another toy example", {
    
    grid <- rep(0, 10)
    start_and_end_minus_buffer <- c(1, 10) ## 1-based
    outputSNPBlockSize <- 1000
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )

    expect_equivalent(blocks_for_output[, "snp_start_1_based"], 1)
    expect_equivalent(blocks_for_output[, "snp_end_1_based"], 10)    
    expect_equivalent(blocks_for_output[, "grid_start_0_based"], 0)
    expect_equivalent(blocks_for_output[, "grid_end_0_based"], 0)    
    
})

test_that("can determine how assign SNPs to blocks for tiny outputSNPBlockSize", {
    
    L <- 1:50
    regionStart <- 10
    regionEnd <- 30
    buffer <- 2
    gridWindowSize <- 4
    outputSNPBlockSize <- 1
    out <- assign_positions_to_grid(L, gridWindowSize)
    grid <- out$grid
    start_and_end_minus_buffer <- shrink_region(regionStart, regionEnd, buffer, L, validate = FALSE)$start_and_end_minus_buffer
    
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )
    
    ## expect each block has >1 grid
    ## only exception is for single region zones
    expect_equal(
        sum(!(blocks_for_output[, "grid_end_0_based"] - blocks_for_output[, "grid_start_0_based"]) > 0),
        0
    )
    
})


test_that("can determine how assign SNPs to blocks for output on yet another toy example", {

    L <- 1:50
    regionStart <- 10
    regionEnd <- 30
    buffer <- 2
    gridWindowSize <- 4
    outputSNPBlockSize <- 8
    out <- assign_positions_to_grid(L, gridWindowSize)
    grid <- out$grid
    start_and_end_minus_buffer <- shrink_region(regionStart, regionEnd, buffer, L, validate = FALSE)$start_and_end_minus_buffer
    
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )
    expect_equivalent(blocks_for_output[1, "snp_start_1_based"], head(start_and_end_minus_buffer, 1))
    expect_equivalent(blocks_for_output[nrow(blocks_for_output), "snp_end_1_based"], tail(start_and_end_minus_buffer, 1))

})


test_that("can determine how assign SNPs to blocks for output", {

    for(gridWindowSize in c(24, NA)) {
        L <- 1:1000
        gridWindowSize <- 24
        outputSNPBlockSize <- 57
        start_and_end_minus_buffer <- c(50, 950) ## 1-based which SNPs in region to output
    
        out <- assign_positions_to_grid(L, gridWindowSize)
        grid <- out$grid
        
        blocks_for_output <- determine_snp_and_grid_blocks_for_output(
            grid = grid,
            start_and_end_minus_buffer = start_and_end_minus_buffer,
            outputSNPBlockSize = outputSNPBlockSize
        )
        
        temp <- array(NA, 1000)
        for(i in 1:nrow(blocks_for_output)) {
            w <- blocks_for_output[i, 1:2]
            temp[w[1]:w[2]] <- i
        }
        
        ## check that no grid spans blocks
        for(i in 0:(out$nGrids - 1)) {
            w <- unique(temp[which(grid == i)])
            expect(length(w), 1)
        }
        
        ## check things
        expect_equal(sum(is.na(temp[1:(start_and_end_minus_buffer[1] - 1)]) == FALSE), 0)
        expect_equal(sum(is.na(temp[(start_and_end_minus_buffer[2] + 1):length(temp)]) == FALSE), 0)

    }

})

test_that("can run forwardBackward, then re-run using list of forward and backward probabilities", {

    n_snps <- 20
    L <- 1:n_snps
    gridWindowSize <- 3
    K <- 2
    regionStart <- 4
    regionEnd <- 12
    buffer <- 2
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    nGrids <- out$nGrids ## 1-based
    
    sampleRead <- function(x) {
        y <- x + c(0, 1, 2, 3)
        y <- y[y < (n_snps -1)]
        return(list(
            length(y) -1, round(median(grid[y + 1])),
            matrix(sample(c(10, -10), length(y), replace = TRUE), ncol = 1),
            matrix(y, ncol = 1)
        ))
    }
    
    sampleReads <- lapply(1:17, sampleRead)
    outputSNPBlockSize <- 4
    
    eHapsCurrent_t <- array(runif(n_snps * K), c(K, n_snps))
    sigmaCurrent <- rep(0.999, nGrids - 1)
    alphaMatCurrent_t <- array(1 / K, c(K, nGrids - 1))
    priorCurrent <- runif(K) / K

    ## run with, without grids
    for (method in c("diploid", "pseudoHaploid", "diploid-inbred")) {

        if (method == "pseudoHaploid") {
            pRgivenH1 <- runif(length(sampleReads))
            pRgivenH2 <- runif(length(sampleReads))
        }
        transMatRate_t <- get_transMatRate(
            method = method,
            sigmaCurrent
        )

        out1 <- run_forward_backwards(
            sampleReads = sampleReads,
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,            
            method = method,
            K = K,
            priorCurrent = priorCurrent,
            alphaMatCurrent_t = alphaMatCurrent_t,
            eHapsCurrent_t = eHapsCurrent_t,
            transMatRate_t = transMatRate_t,
            whatToReturnOriginal = 0
        )

        grid_starts <- c(1, 2)
        grid_ends <- c(3, 5)
        alphaBetaBlock <- lapply(out1$fbsoL, function(fbso) {
            return(
                list(
                    alphaHatBlocks_t = fbso$alphaHat_t[, grid_starts + 1, drop = FALSE],
                    betaHatBlocks_t = fbso$betaHat_t[, grid_ends + 1, drop = FALSE]
                )
            )
        })
        
        ## do between 0-based grids 2 and 4 inclusive
        for(i_snp_block_for_alpha_beta in 1:2) {
            grid_start <- grid_starts[i_snp_block_for_alpha_beta]
            grid_end <- grid_ends[i_snp_block_for_alpha_beta]            
            x <- sapply(sampleReads, function(x) x[[2]])
            whichReads <- which((grid_start <= x) & (x <= grid_end))

            out2 <- run_forward_backwards(
                sampleReads = sampleReads[whichReads],
                pRgivenH1 = pRgivenH1[whichReads],
                pRgivenH2 = pRgivenH2[whichReads],            
                method = method,
                K = K,
                priorCurrent = priorCurrent,
                alphaMatCurrent_t = alphaMatCurrent_t[, 1 + grid_start:(grid_end - 1), drop = FALSE],
                eHapsCurrent_t = eHapsCurrent_t,
                transMatRate_t = transMatRate_t[, 1 + grid_start:(grid_end - 1), drop = FALSE],
                run_fb_subset = TRUE,
                alphaBetaBlock = alphaBetaBlock,
                run_fb_grid_offset = grid_start,
                suppressOutput = 1,
                i_snp_block_for_alpha_beta = i_snp_block_for_alpha_beta
            )
            
            for(iNor in 1:length(out1[[1]])) {
                expect_equal(
                    out1$fbsoL[[1]]$gamma_t[, grid_start:grid_end + 1],
                    out2$fbsoL[[1]]$gamma_t
                )
            }
        }
        
    }

})
