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
        
        transMatRate_t_H <- get_transMatRate(method = "diploid-inbred", sigmaCurrent)
        transMatRate_t_D <- get_transMatRate(method = "diploid", sigmaCurrent)
        
        for(output_haplotype_dosages in c(FALSE, TRUE)) {
            for (method in c("diploid", "pseudoHaploid", "diploid-inbred")) {            
                
                if (method == "pseudoHaploid") {
                    pRgivenH1 <- runif(length(sampleReads))
                    pRgivenH2 <- runif(length(sampleReads))
                }
                transMatRate_t <- get_transMatRate(method = method, sigmaCurrent)
                
                out1 <- run_forward_backwards(
                    sampleReads = sampleReads,
                    pRgivenH1 = pRgivenH1,
                    pRgivenH2 = pRgivenH2,            
                    method = method,
                    K = K,
                    priorCurrent = priorCurrent,
                    alphaMatCurrent_t = alphaMatCurrent_t,
                    eHapsCurrent_t = eHapsCurrent_t,
                    transMatRate_t_H = transMatRate_t_H,
                    transMatRate_t_D = transMatRate_t_D,
                    blocks_for_output = blocks_for_output,
                    generate_fb_snp_offsets = TRUE,
                    return_genProbs = TRUE, ## might not be for all methods
                    grid = grid,
                    output_haplotype_dosages = output_haplotype_dosages
                )$fbsoL
                if (method == "pseudoHaploid") {
                    gammaEK_t_all <- out1[[1]][["gammaEK_t"]] + out1[[2]][["gammaEK_t"]]
                } else {
                    gammaEK_t_all <- out1[[1]][["gammaEK_t"]]                            
                }
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
                        transMatRate_t_H = transMatRateLocal_t, ## ugh it will decide what to do
                        transMatRate_t_D = transMatRateLocal_t, ## ugh                    
                        run_fb_subset = TRUE,
                        alphaBetaBlock = alphaBetaBlock,
                        run_fb_grid_offset = first_grid_in_region,
                        suppressOutput = 1,
                        i_snp_block_for_alpha_beta = i_output_block,
                        return_genProbs = TRUE, ## might not be for all methods
                        grid = grid,
                        snp_start_1_based = first_snp_in_region,
                        snp_end_1_based = last_snp_in_region,
                        output_haplotype_dosages = TRUE
                    )$fbsoL
                    if (output_haplotype_dosages) {
                        if (method == "pseudoHaploid") {
                            gammaEK_t_local <- out2[[1]][["gammaEK_t"]] + out2[[2]][["gammaEK_t"]]
                        } else {
                            gammaEK_t_local <- out2[[1]][["gammaEK_t"]]
                        }
                    }
                    
                    gp_t_local <- calculate_gp_t_from_fbsoL(
                        eHapsCurrent_t = eHapsCurrent_t,
                        grid = grid,
                        method = method,
                        fbsoL = out2,
                        snp_start_1_based = first_snp_in_region,
                        snp_end_1_based = last_snp_in_region,
                        grid_offset_0_based = first_grid_in_region
                    )
                    ## annoying? do same thing for state probabilities?
                    ## these 
                    
                    for(iNor in 1:length(out1[[1]])) {
                        expect_equal(
                            out1[[1]]$gamma_t[, 1 + first_grid_in_region:last_grid_in_region, drop = FALSE],
                            out2[[1]]$gamma_t
                        )
                        expect_equivalent(
                            gp_t_all[, first_snp_in_region:last_snp_in_region], 
                            gp_t_local
                        )
                        if (output_haplotype_dosages) {
                            expect_equivalent(
                                gammaEK_t_all[, first_snp_in_region:last_snp_in_region],
                                gammaEK_t_local
                            )
                        }
                    }
                }
            }
        }
        
    }

})


