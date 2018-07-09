test_that("can determine how assign SNPs to blocks for output on toy example", {
    
    grid <- 0:9 ## 0-based
    start_and_end_minus_buffer <- c(1, 10) ## 1-based
    outputSNPBlockSize <- 1000
    snp_blocks_for_output <- determine_snp_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )

    expected_output <- array(c(1, 10), c(1, 2))
    colnames(expected_output) <- c("snp_start", "snp_end")
    expect_equal(snp_blocks_for_output, expected_output)
    
})

test_that("can determine how assign SNPs to blocks for output", {
    
    L <- 1:1000
    gridWindowSize <- 24
    outputSNPBlockSize <- 57
    start_and_end_minus_buffer <- c(50, 950) ## 1-based which SNPs in region to output
    
    out <- assign_positions_to_grid(L, gridWindowSize)
    grid <- out$grid
    
    snp_blocks_for_output <- determine_snp_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )

    temp <- array(NA, 1000)
    for(i in 1:nrow(snp_blocks_for_output)) {
        w <- snp_blocks_for_output[i, 1:2]
        temp[w[1]:w[2]] <- i
    }
    
    ## check that no grid spans blocks
    for(i in 0:(out$nGrids - 1)) {
        w <- unique(temp[which(grid == i)])
        expect(length(w), 1)
    }

    ## check that
    expect_equal(sum(is.na(temp[1:(start_and_end_minus_buffer[1] - 1)]) == FALSE), 0)
    expect_equal(sum(is.na(temp[(start_and_end_minus_buffer[2] + 1):length(temp)]) == FALSE), 0)    

})

test_that("can run forwardBackward, then re-run using list of forward and backward probabilities", {

    n_snps <- 12
    L <- 1:n_snps
    K <- 2
    sampleRead <- function(x) {
        return(list(
            3, x + 1,
            matrix(c(10, 10, -10, -10), ncol = 1),
            matrix(x + c(0, 1, 2, 3), ncol = 1)
        ))
    }
    sampleReads <- list(sampleRead(0), sampleRead(4), sampleRead(8))
    eHapsCurrent_t <- array(runif(n_snps * K), c(K, n_snps))
    sigmaCurrent <- rep(0.999, n_snps - 1)
    alphaMatCurrent_t <- array(1 / K / K, c(K, n_snps - 1))
    priorCurrent <- runif(K) / K

    outputSNPBlockSize <- 6

    ## run with, without grids
    for (method in c("diploid", "pseudoHaploid", "diploid-inbred")) {

        if (method == "pseudoHaploid") {
            pRgivenH1 <- runif(3)
            pRgivenH2 <- runif(3)
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

        whichReads <- 2
        snp_start <- 5
        snp_end <- 8
        alphaStart <- out1$fbsoL[[1]]$alphaHat_t[, snp_start]
        betaEnd <- out1$fbsoL[[1]]$betaHat_t[, snp_end]

        out2 <- run_forward_backwards(
            sampleReads = sampleReads[whichReads],
            pRgivenH1 = pRgivenH1[whichReads],
            pRgivenH2 = pRgivenH2[whichReads],            
            method = method,
            K = K,
            priorCurrent = priorCurrent,
            alphaMatCurrent_t = alphaMatCurrent_t[, snp_start:(snp_end - 1)],
            eHapsCurrent_t = eHapsCurrent_t,
            transMatRate_t = transMatRate_t[, snp_start:(snp_end - 1)],
            run_fb_subset = TRUE,
            alphaStart = alphaStart,
            betaEnd = betaEnd,
            run_fb_snp_offset = snp_start - 1, ## this is 0-based. snp_start is 1-based
            suppressOutput = 1
        )

        expect_equal(
            out1$fbsoL[[1]]$gamma_t[, snp_start:snp_end],
            out2$fbsoL[[1]]$gamma_t
        )
        
    }


})
