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

    n_snps <- 100
    L <- 1:n_snps
    K <- 2
    regionStart <- 4
    regionEnd <- 48
    buffer <- 2
    inRegion <- L >= (regionStart) & L <= (regionEnd)
    w <- which(inRegion)
    start_and_end_minus_buffer <- c(head(w, 1), tail(w, 1)) ## first SNP
    outputSNPBlockSize <- 4    
    
    sampleRead <- function(x, grid) {
        y <- x + c(0, 1, 2, 3)
        y <- y[y < (n_snps -1)]
        return(list(
            length(y) -1, round(median(grid[y + 1])),
            matrix(sample(c(10, -10), length(y), replace = TRUE), ncol = 1),
            matrix(y, ncol = 1)
        ))
    }

    for(gridWindowSize in c(NA, 3)) {

        out <- assign_positions_to_grid(
            L = L,
            gridWindowSize = gridWindowSize
        )
        grid <- out$grid
        nGrids <- out$nGrids ## 1-based

        blocks_for_output <- determine_snp_and_grid_blocks_for_output(
            grid = grid,
            start_and_end_minus_buffer = start_and_end_minus_buffer,
            outputSNPBlockSize = outputSNPBlockSize
        )
        
        ## depends on grid
        sampleReads <- lapply(1:(n_snps - 3), sampleRead, grid = grid)
        eHapsCurrent_t <- array(runif(n_snps * K), c(K, n_snps))
        sigmaCurrent <- rep(0.999, nGrids - 1)
        alphaMatCurrent_t <- array(1 / K, c(K, nGrids - 1))
        priorCurrent <- runif(K) / K
        

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
                whatToReturnOriginal = 0,
                blocks_for_output = blocks_for_output,
                generate_fb_snp_offsets = TRUE
            )$fbsoL
            gp_t_all <- calculate_gp_t_from_fbsoL(
                eHapsCurrent_t = eHapsCurrent_t,
                grid = grid,
                method = method,
                fbsoL = out1
            )
            alphaBetaBlock <- lapply(out1, function(x) x$alphaBetaBlocks)
            
            blocks_in_vector_form <- array(NA, nGrids) ## want 0, 1, 2, etc, depending on which output block
            for(i in 1:nrow(blocks_for_output)) {
                s2 <- blocks_for_output[i, "grid_start_0_based"]
                e2 <- blocks_for_output[i, "grid_end_0_based"]
                blocks_in_vector_form[1 + s2:e2] <- i
            }
            rse <- determine_starts_ends_whats_for_sampleReads(
                sampleReads = sampleReads,
                blocks_in_vector_form = blocks_in_vector_form
            )
            
            for(i_output_block in 1:nrow(blocks_for_output)) {
                
                first_snp_in_region <- blocks_for_output[i_output_block, "snp_start_1_based"]
                last_snp_in_region <- blocks_for_output[i_output_block, "snp_end_1_based"]
                snps_in_output_block <- first_snp_in_region:last_snp_in_region
                first_grid_in_region <- blocks_for_output[i_output_block, "grid_start_0_based"]
                last_grid_in_region <- blocks_for_output[i_output_block, "grid_end_0_based"]
                if (first_grid_in_region < last_grid_in_region) {
                    grids_to_use <- first_grid_in_region:(last_grid_in_region - 1)
                    alphaMatCurrentLocal_t <- alphaMatCurrent_t[, 1 + grids_to_use, drop = FALSE]
                    transMatRateLocal_t <- transMatRate_t[, 1 + grids_to_use, drop = FALSE]
                }
                whichReads <- rse$starts[i_output_block]:rse$ends[i_output_block]

                out2 <- run_forward_backwards(
                    sampleReads = sampleReads[whichReads],
                    pRgivenH1 = pRgivenH1[whichReads],
                    pRgivenH2 = pRgivenH2[whichReads],            
                    method = method,
                    K = K,
                    priorCurrent = priorCurrent,
                    alphaMatCurrent_t = alphaMatCurrentLocal_t,
                    eHapsCurrent_t = eHapsCurrent_t,
                    transMatRate_t = transMatRateLocal_t,
                    run_fb_subset = TRUE,
                    alphaBetaBlock = alphaBetaBlock,
                    run_fb_grid_offset = first_grid_in_region,
                    suppressOutput = 1,
                    i_snp_block_for_alpha_beta = i_output_block
                )$fbsoL
                
                gp_t_local <- calculate_gp_t_from_fbsoL(
                    eHapsCurrent_t = eHapsCurrent_t,
                    grid = grid,
                    method = method,
                    fbsoL = out2,
                    snp_start_1_based = first_snp_in_region,
                    snp_end_1_based = last_snp_in_region,
                    grid_offset_0_based = first_grid_in_region
                )
            
                for(iNor in 1:length(out1[[1]])) {
                    expect_equal(
                        out1[[1]]$gamma_t[, 1 + first_grid_in_region:last_grid_in_region, drop = FALSE],
                        out2[[1]]$gamma_t
                    )
                    expect_equivalent(
                        gp_t_all[, first_snp_in_region:last_snp_in_region], 
                        gp_t_local
                    )
                }
            }
        }
        
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
