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



test_that("can properly assign reads to those blocks in which they are needed", {


    L <- 1:300
    regionStart <- 50
    regionEnd <- 224
    buffer <- 14
    gridWindowSize <- NA
    outputSNPBlockSize <- 27

    f <- function(grid) {
        lapply(grid, function(i) {
            return(list(0, i, matrix(10, ncol = 1), matrix(i, ncol = 1)))
        })
    }
    
    for(gridWindowSize in c(NA, 15)) {

        out <- assign_positions_to_grid(
            L = L,
            gridWindowSize = gridWindowSize
        )
        grid <- out$grid
        nGrids <- out$nGrids ## 1-based

        inRegion <- L >= (regionStart) & L <= (regionEnd)
        w <- which(inRegion)
        start_and_end_minus_buffer <- c(head(w, 1), tail(w, 1)) ## first SNP
        blocks_for_output <- determine_snp_and_grid_blocks_for_output(
            grid = grid,
            start_and_end_minus_buffer = start_and_end_minus_buffer,
            outputSNPBlockSize = outputSNPBlockSize
        )
        
        blocks_in_vector_form <- get_blocks_in_vector_form(blocks_for_output, nGrids)

        ## make fake sampleReads with reads everywhere!
        sampleReads <- f(grid)
        out <- determine_starts_ends_whats_for_sampleReads(sampleReads, blocks_in_vector_form)
        central_snp_in_read <- sapply(sampleReads, function(x) x[[2]]) + 1 ## 1-based

        ## check that all reads in each output block are represented
        for(i_output_block in 1:nrow(blocks_for_output)) {
            s <- out$starts[i_output_block]
            e <- out$ends[i_output_block]
            w <- out$whats[i_output_block]
            s2 <- blocks_for_output[i_output_block, "grid_start_0_based"] + 1
            e2 <- blocks_for_output[i_output_block, "grid_end_0_based"] + 1
            expect_equal(
                which((s2 <= central_snp_in_read) & (central_snp_in_read <= e2)),
                s:e
            )
        }

    }

})


