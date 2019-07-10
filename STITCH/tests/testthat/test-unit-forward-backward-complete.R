test_that("can run forwardBackward, then re-run using list of forward and backward probabilities", {

    nSNPs <- 100
    L <- 1:nSNPs
    K <- 2
    S <- 3
    regionStart <- 4
    regionEnd <- 48
    buffer <- 2
    inRegion <- L >= (regionStart) & L <= (regionEnd)
    w <- which(inRegion)
    start_and_end_minus_buffer <- c(head(w, 1), tail(w, 1)) ## first SNP
    outputSNPBlockSize <- 4    
    
    sampleRead <- function(x, grid) {
        y <- x + c(0, 1, 2, 3)
        y <- y[y < (nSNPs -1)]
        return(list(
            length(y) -1, round(median(grid[y + 1])),
            matrix(sample(c(10, -10), length(y), replace = TRUE), ncol = 1),
            matrix(y, ncol = 1)
        ))
    }
    gridWindowSize <- NA    

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
        sampleReads <- lapply(1:(nSNPs - 3), sampleRead, grid = grid)

        eHapsCurrent_tc <- array(NA, c(K, nSNPs, S))
        alphaMatCurrent_tc <- array(NA, c(K, nGrids - 1, S))
        eHapsCurrent_t <-
            0.1 * array(runif(K * nSNPs), c(K, nSNPs)) +
            0.9 * array(sample(c(0, 1), K * nSNPs, replace = TRUE), c(K, nSNPs))
        for(s in 1:S) {
            eHapsCurrent_tc[, , s] <- 0.90 * eHapsCurrent_t + 0.10 * runif(nSNPs * K)
            m <- array(runif(K * (nGrids - 1)), c(K, (nGrids - 1)))
            alphaMatCurrent_tc[, , s] <- t(t(m) / colSums(m))
        }
        sigmaCurrent_m <- array(0.9 + 0.1 * runif((nGrids - 1) * S), c(nGrids - 1, S))
        priorCurrent_m <- array(1 / K, c(K, S))
        
        transMatRate_tc_H <- get_transMatRate_m(method = "diploid-inbred", sigmaCurrent_m)
        transMatRate_tc_D <- get_transMatRate_m(method = "diploid", sigmaCurrent_m)

        output_haplotype_dosages <- FALSE
        
        for(output_haplotype_dosages in c(FALSE, TRUE)) {

            method <- "diploid"
            
            for (method in c("diploid", "pseudoHaploid", "diploid-inbred")) {

                if (method == "pseudoHaploid") {
                    out <- get_default_hapProbs(
                        pseudoHaploidModel = 9,
                        sampleReads = sampleReads,
                        S = S
                    )
                    pRgivenH1_m <- out$pRgivenH1_m
                    pRgivenH2_m <- out$pRgivenH2_m
                }

                fbsoL1 <- run_forward_backwards(
                    sampleReads = sampleReads,
                    pRgivenH1_m = pRgivenH1_m,
                    pRgivenH2_m = pRgivenH2_m,
                    method = method,
                    priorCurrent_m = priorCurrent_m,
                    alphaMatCurrent_tc = alphaMatCurrent_tc,
                    eHapsCurrent_tc = eHapsCurrent_tc,
                    transMatRate_tc_H = transMatRate_tc_H,
                    transMatRate_tc_D = transMatRate_tc_D,
                    blocks_for_output = blocks_for_output,
                    generate_fb_snp_offsets = TRUE,
                    return_genProbs = TRUE, ## might not be for all methods
                    return_hapDosage = TRUE,
                    return_gamma = TRUE,
                    return_extra = TRUE,                    
                    grid = grid,
                    output_haplotype_dosages = output_haplotype_dosages
                )
                
                if (method == "pseudoHaploid") {
                    gammaEK_t_all <- fbsoL1[[1]][["gammaEK_t"]] + fbsoL1[[2]][["gammaEK_t"]]
                } else {
                    gammaEK_t_all <- fbsoL1[[1]][["gammaEK_t"]]                            
                }
                gp_t_all <- calculate_gp_t_from_fbsoL(
                    fbsoL = fbsoL1,
                    method = method
                )
                list_of_alphaBetaBlocks <- lapply(fbsoL1, function(x) x$list_of_alphaBetaBlocks)

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

                i_output_block <- 1
                
                for(i_output_block in 1:nrow(blocks_for_output)) {
                    
                    first_snp_in_region <- blocks_for_output[i_output_block, "snp_start_1_based"]
                    last_snp_in_region <- blocks_for_output[i_output_block, "snp_end_1_based"]
                    snps_in_output_block <- first_snp_in_region:last_snp_in_region
                    first_grid_in_region <- blocks_for_output[i_output_block, "grid_start_0_based"]
                    last_grid_in_region <- blocks_for_output[i_output_block, "grid_end_0_based"]
                    if (first_grid_in_region < last_grid_in_region) {
                        grids_to_use <- first_grid_in_region:(last_grid_in_region - 1)
                        alphaMatCurrentLocal_tc <- alphaMatCurrent_tc[, 1 + grids_to_use, , drop = FALSE]
                        transMatRateLocal_tc_H <- transMatRate_tc_H[, 1 + grids_to_use, , drop = FALSE]
                        transMatRateLocal_tc_D <- transMatRate_tc_D[, 1 + grids_to_use, , drop = FALSE]
                    }
                    whichReads <- rse$starts[i_output_block]:rse$ends[i_output_block]

                    list_of_alphaBetaBlocks <- lapply(fbsoL1, function(x) x[["list_of_alphaBetaBlocks"]])
                    fbsoL2 <- run_forward_backwards(
                        sampleReads = sampleReads[whichReads],
                        pRgivenH1_m = pRgivenH1_m[whichReads, , drop = FALSE],
                        pRgivenH2_m = pRgivenH2_m[whichReads, , drop = FALSE],
                        method = method,
                        priorCurrent_m = priorCurrent_m,
                        alphaMatCurrent_tc = alphaMatCurrentLocal_tc,
                        eHapsCurrent_tc = eHapsCurrent_tc,
                        transMatRate_tc_H = transMatRateLocal_tc_H,
                        transMatRate_tc_D = transMatRateLocal_tc_D,
                        run_fb_subset = TRUE,
                        list_of_alphaBetaBlocks = list_of_alphaBetaBlocks,
                        run_fb_grid_offset = first_grid_in_region,
                        suppressOutput = 1,
                        i_snp_block_for_alpha_beta = i_output_block,
                        return_genProbs = TRUE, ## might not be for all methods
                        return_hapDosage = TRUE,                        
                        return_gamma = TRUE,
                        return_extra = TRUE,
                        grid = grid,
                        snp_start_1_based = first_snp_in_region,
                        snp_end_1_based = last_snp_in_region,
                        output_haplotype_dosages = TRUE
                    )
                    
                    if (output_haplotype_dosages) {
                        if (method == "pseudoHaploid") {
                            gammaEK_t_local <- fbsoL2[[1]][["gammaEK_t"]] + fbsoL2[[2]][["gammaEK_t"]]
                        } else {
                            gammaEK_t_local <- fbsoL2[[1]][["gammaEK_t"]]
                        }
                    }

                    ## OK - I have not done this at all!
                    ## save to alphaStartM and same for betaEndM
                    ## then pass through!
                    ## save(fbsoL1, fbsoL2, file = "~/temp.RData")

                    ## ##
                    ## print(list_of_alphaBetaBlocks[[1]][[3]]$alphaHatBlocks_t[, 1])
                    ## print("What it definitely is")
                    ## print(fbsoL1[[1]]$alphaHat_t[, 1 + first_grid_in_region:last_grid_in_region, drop = FALSE])
                    ## print("What seems to have been used")
                    ## print(fbsoL2[[1]]$alphaHat_t) ## so this is not as expected (for s=4) - why
                    
                    ##stop("FFFFFFFFF")
                    ##load("~/temp.RData")

                    gp_t_local <- calculate_gp_t_from_fbsoL( 
                        fbsoL = fbsoL2,
                        method = method
                    )
                    
                    for(iNor in 1:length(fbsoL1)) {
                        
                        expect_equal(
                            fbsoL1[[iNor]]$list_of_gamma_t[[1]][, 1 + first_grid_in_region:last_grid_in_region, drop = FALSE],
                            fbsoL2[[iNor]]$list_of_gamma_t[[1]]
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


