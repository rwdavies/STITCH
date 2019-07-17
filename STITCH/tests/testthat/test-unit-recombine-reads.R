test_that("can split a read", {

    L <- c(1, 10, 100, 1000)
    grid <- c(0, 1, 2, 3)
    sampleRead <- list(
        3, 0,
        matrix(c(10, 10, -10, -10), ncol = 1),
        matrix(c(0, 1, 2, 3), ncol = 1)
    )
    
    set.seed(1)
    out <- get_sampleRead_from_SNP_i_to_SNP_j(
        sampleRead,
        1,
        2,
        L,
        grid
    )
    expect_equal(
        out,
        list(
            1, 0,
            matrix(c(10, 10), ncol = 1),
            matrix(c(0, 1), ncol = 1)
        )
    )

    set.seed(1)
    out <- get_sampleRead_from_SNP_i_to_SNP_j(
        sampleRead,
        3,
        4,
        L,
        grid
    )
    expect_equal(
        out,
        list(
            1, 2,
            matrix(c(-10, -10), ncol = 1),
            matrix(c(2, 3), ncol = 1)
        )
    )
    
})



test_that("a good read does not want to get split", {

    ## currently only working on length 4 and above
    K <- 2
    nSNPs <- 4
    eHapsCurrent_t <- rbind(rep(0.001, nSNPs), rep(0.999, nSNPs))
    ## - is ref, + is alt
    a_read <- list(
        3, 0,
        matrix(c(-10, -10, -10, -10), ncol = 1),
        matrix(c(0, 1, 2, 3), ncol = 1)            
    )
    sampleReads <- list(a_read)
    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    expect_equal(
        w,
        integer(0)
    )

})


test_that("a bad read wants to get split", {

    K <- 2
    nSNPs <- 4
    eHapsCurrent_t <- rbind(rep(0.001, nSNPs), rep(0.999, nSNPs))
    a_read <- list(
        3, 0,
        matrix(c(-10, 10, -10, 10), ncol = 1),
        matrix(c(0, 1, 2, 3), ncol = 1)            
    )
    sampleReads <- list(a_read)
    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    expect_equal(
        w,
        1
    )

})

test_that("correctly split a bad read", {

    iSample <- 1
    regionName <- "a_region"
    K <- 2
    S <- 1
    nSNPs <- 4
    L <- 1:nSNPs
    gridWindowSize <- NA
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids

    
    eHapsCurrent_t <- rbind(rep(0.001, nSNPs), rep(0.999, nSNPs))
    gammaK_t <- rbind(rep(0.001, nSNPs), rep(0.999, nSNPs))
    gammaK_t[1, 3:4] <- 0.999
    gammaK_t[2, 3:4] <- 0.001
    readS <- list(
        0, 0, matrix(-10, ncol = 1), matrix(0, ncol = 1)
    )
    readE <- list(
        0, 3, matrix(-10, ncol = 1), matrix(3, ncol = 1)
    )
    sampleRead <- list(
        3, 0,
        matrix(c(10, 10, -10, -10), ncol = 1),
        matrix(c(0, 1, 2, 3), ncol = 1)            
    )
    sampleReads <- list(readS, sampleRead, readE)

    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    expect_equal(w, 2)

    set.seed(1)

    method <- "pseudoHaploid"
    out <- get_default_hapProbs(pseudoHaploidModel = 9, sampleReads = sampleReads, S = S)
    pRgivenH1_m <- out$pRgivenH1_m
    pRgivenH2_m <- out$pRgivenH2_m
    pRgivenH1_m[] <- 0.9
    pRgivenH2_m[] <- 0.4    
    srp <- out$srp

    tempdir <- tempdir()
    
    out <- findRecombinedReadsPerSample(
        gammaK_t = gammaK_t,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K,
        L = L,
        iSample = iSample,
        sampleReads = sampleReads,
        tempdir = tempdir,
        regionName = regionName,
        grid = grid,
        method = method,
        pRgivenH1_m = pRgivenH1_m,
        pRgivenH2_m = pRgivenH2_m,
        srp = srp
    )
    
    expect_equal(out$readsSplit, 1)
    expect_equal(out$readsTotal, 4)

    load(file_sampleReads(tempdir, iSample, regionName))
    load(file_sampleProbs(tempdir, iSample, regionName))
    
    ## trust get_sampleRead_from_SNP_i_to_SNP_j function
    ## as it is unit tested above
    set.seed(1)
    expect_equal(
        sampleReads[[2]],
        get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, 1, 2, L, grid
        )
    )

    ## do not bother setting seed - is second sample call
    new_read_2 <- get_sampleRead_from_SNP_i_to_SNP_j(
        sampleRead, 3, 4, L, grid
    )

    if (length(RNGkind()) > 2 && RNGkind()[3] == "Rejection") {
        new_read_2[[2]] <- 3
        ## >= 3.6.0
        print(sampleReads)
        print(new_read_2)
        expect_equal(
            sampleReads[[4]], ## this test involves a sample call - is still fine with >=R3.6.0new_
            new_read_2
        )
    } else {
        expect_equal(
            sampleReads[[3]], ## this test involves a sample call - is still fine with >=R3.6.0new_
            new_read_2
        )
    }

    expect_equal(srp, sapply(sampleReads, function(x) x[[2]]))
    expect_equal(nrow(pRgivenH1_m), length(sampleReads))
    expect_equal(nrow(pRgivenH2_m), length(sampleReads))    
    
})


test_that("correctly split a bad read with grid mode", {

    iSample <- 1
    regionName <- "a_region"
    K <- 2
    nSNPs <- 8
    gridWindowSize <- 2
    L <- 1:nSNPs    
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids

    
    eHapsCurrent_t <- rbind(
        rep(0.001, nSNPs),
        rep(0.999, nSNPs)
    )
    gammaK_t <- rbind(
        rep(0.001, nGrids),
        rep(0.999, nGrids)
    )
    gammaK_t[1, 3:4] <- 0.999
    gammaK_t[2, 3:4] <- 0.001
    readS <- list(
        0, 0, matrix(-10, ncol = 1), matrix(0, ncol = 1)
    )
    readE <- list(
        0, 7, matrix(-10, ncol = 1), matrix(7, ncol = 1)
    )
    sampleRead <- list(
        3, 3,
        matrix(c(10, 10, -10, -10), ncol = 1),
        matrix(c(1, 2, 3, 4), ncol = 1)            
    )
    sampleReads <- list(readS, sampleRead, readE)
    sampleReads <- snap_sampleReads_to_grid(sampleReads, grid)

    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    expect_equal(w, 2)

    set.seed(1)

    tempdir <- tempdir()
    out <- findRecombinedReadsPerSample(
        gammaK_t = gammaK_t,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K,
        L = L,
        iSample = iSample,
        sampleReads = sampleReads,
        tempdir = tempdir,
        regionName = regionName,
        grid = grid
    )
    
    expect_equal(out$readsSplit, 1)
    expect_equal(out$readsTotal, 4)

    load(file_sampleReads(tempdir, iSample, regionName))
    ## trust get_sampleRead_from_SNP_i_to_SNP_j function
    ## as it is unit tested above
    set.seed(1)
    expect_equal(
        sampleReads[[2]],
        get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, 1, 2, L, grid
        )
    )

    if (length(RNGkind()) > 2 && RNGkind()[3] == "Rejection") {
        set.seed(4) ## not the same seed as above - is second call
    } else {
        set.seed(1)
    }
    new_read_2 <- get_sampleRead_from_SNP_i_to_SNP_j(
        sampleRead, 3, 4, L, grid
    )    
    expect_equal(
        sampleReads[[3]],
        new_read_2
    )

})
