n_snps <- 5
chr <- 1

phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = 4,
    seed = 1,
    chr = chr,
    K = 2,
    phasemaster = phasemaster
)



test_that("can validate gridWindowSize", {
    expect_null(validate_gridWindowSize(as.numeric(1000)))
    expect_null(validate_gridWindowSize(as.integer(1000)))    
    expect_null(validate_gridWindowSize(10))    
    expect_null(validate_gridWindowSize(3))
    expect_null(validate_gridWindowSize(NA))    
    expect_error(validate_gridWindowSize(1))
    expect_error(validate_gridWindowSize(10.5))    
    expect_error(validate_gridWindowSize("1000"))
})

test_that("can assign physical positions to grid", {

    out <- assign_positions_to_grid(
        L = c(1, 10, 11, 20),
        gridWindowSize = 10
    )
    expect_equal(out$grid, c(0, 0, 1, 1))
    expect_equal(out$grid_distances, 10)
    expect_equal(out$L_grid, c(5, 15))
    expect_equal(out$nGrids, 2)

})

test_that("distances on the grid are OK", {

    out <- assign_positions_to_grid(
        L = 1e6 + c(1, 15456, 95123),
        gridWindowSize = 1e4
    )
    expect_equal(out$grid, c(0, 1, 2))
    expect_equal(out$grid_distances, c(1e4, 8 * 1e4))
    expect_equal(out$L_grid, 1e6 + c(5000, 15000, 95000))
    expect_equal(out$nGrids, 3)

})



test_that("removal of buffer from grid makes sense", {

    L <- 1:20
    n_snps <- 20
    regionStart <- 5
    regionEnd <- 16
    buffer <- 5
    K <- 4

    for(gridWindowSize in c(NA, 1, 3)) {
        
        out <- assign_positions_to_grid(
            L = L,
            gridWindowSize = gridWindowSize
        )
    
        alphaMatCurrent <- array(0, c(out$nGrids, K))
        
        out <- remove_buffer_from_variables(
            L = L,
            regionStart = regionStart,
            regionEnd = regionEnd,
            grid = out$grid,
            grid_distances = out$grid_distances,
            alphaMatCurrent = alphaMatCurrent,
            L_grid = out$L_grid,
            nGrids = out$nGrids,
            gridWindowSize = gridWindowSize,
            verbose = FALSE
        )
        
        expect_equal(length(out$grid), length(regionStart:regionEnd))
        ## 5-6, 7-9, 10-12, 13-15, 16
        if (is.na(gridWindowSize) == FALSE) {
            if (gridWindowSize == 3) {
                expect_equal(out$nGrids, 5)
                expect_equal(out$grid_distances, rep(gridWindowSize, 4)) ## no spacing
                expect_equal(out$L_grid, 4.5 + 3 * 0:4)
                expect_equal(nrow(out$alphaMatCurrent), 4)
            }
        }
        
    }
        

})




test_that("can use grid", {

    n_snps <- 10 ## set to 10000 to check times better
    K <- 20

    phasemaster <- matrix(
        c(rep(0, n_snps), rep(1, n_snps)),
        ncol = K
    )
    data_package <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_snps * 2,
        seed = 2,
        chr = 10,
        K = K,
        phasemaster = phasemaster,
        reads_span_n_snps = 3,
        n_cores = 1
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = data_package$L,
        pos = data_package$pos,
        nSNPs = data_package$nSNPs,
        bam_files = data_package$bam_files,
        N = 1,
        sampleNames = data_package$sample_names,
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = data_package$chr,
        chrStart = 1,
        chrEnd = max(data_package$pos[, 2]) + 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    L <- data_package$pos[, 2]
    eHaps <- array(runif(n_snps * K), c(n_snps, K))

    ## for now, based on physical position?
    gridWindowSize <- 3
    out <- assign_positions_to_grid(L, gridWindowSize)
    grid <- out$grid
    nGrids <- out$nGrids
    sigma <- runif(nGrids - 1)
    alphaMat <- array(runif((nGrids - 1) * K), c(nGrids - 1, K))
    x <- sigma
    transMatRate <- cbind(x ** 2, x * (1 - x), (1 - x) ** 2)
    pi <- runif(K) / K

    ## also, update?
    sampleReads <- snap_sampleReads_to_grid(sampleReads, grid)

    transMatRate_t <- get_transMatRate(
        method = "diploid",
        sigmaCurrent = sigma
    )

    out <- run_forward_backwards(
        sampleReads = sampleReads,
        priorCurrent = pi,
        transMatRate_t = transMatRate_t,
        alphaMatCurrent_t = t(alphaMat),
        eHapsCurrent_t = t(eHaps),
        method = "diploid"
    )$fbsoL[[1]]
    
    expect_equal(ncol(out$gamma_t), nGrids)
    expect_equal(min(out$gamma_t) >= 0, TRUE)
    expect_equal(max(out$gamma_t) <= 1, TRUE)    

    pRgivenH1L <- runif(length(sampleReads))
    pRgivenH2L <- runif(length(sampleReads))

    transMatRate_t <- get_transMatRate(
        method = "pseudoHaploid",
        sigmaCurrent = sigma
    )

    for(method in c("pseudoHaploid", "diploid-inbred")) {
        out <- run_forward_backwards(
            sampleReads = sampleReads,
            priorCurrent = pi,
            transMatRate_t = t(transMatRate),
            alphaMatCurrent_t = t(alphaMat),
            eHapsCurrent_t = t(eHaps),
            pRgivenH1 = pRgivenH1L,
            pRgivenH2 = pRgivenH2L,
            method = method
        )$fbsoL[[1]]
        
        expect_equal(ncol(out$gamma_t), nGrids)
        expect_equal(min(out$gamma_t) >= 0, TRUE)
        expect_equal(max(out$gamma_t) <= 1, TRUE)
    }
        
})


test_that("can calculate fbd dosage", {

    nSNPs <- 10
    K <- 3
    KK <- K * K
    eHapsCurrent_t <- array(runif(K * nSNPs), c(K, nSNPs))
    gamma_t <- array(runif(KK * nSNPs), c(KK, nSNPs))
    grid <- 0:(nSNPs - 1)
    nGrids <- nSNPs ## since no grid

    out <- rcpp_calculate_fbd_dosage(
        eHapsCurrent_t = eHapsCurrent_t,
        gamma_t = gamma_t,
        grid = grid,
        snp_start_1_based = 1,
        snp_end_1_based = length(grid)
    )

    expect_equal(length(out$dosage), nSNPs)
    expect_equal(ncol(out$genProbs_t), nSNPs)
    expect_equal(nrow(out$genProbs_t), 3)
    expect_equal(sum(out$dosage == 0), 0)
    expect_equal(sum(out$genProbs_t== 0), 0)


})


test_that("can calculate fbd dosage using grid", {

    nSNPs <- 10
    nGrids <- 3
    K <- 3
    KK <- K * K
    eHapsCurrent_t <- array(runif(K * nSNPs), c(K, nSNPs))
    gamma_t <- array(runif(KK * nGrids), c(KK, nGrids))
    grid <- c(0, 0, 0, 1, 1, 1, 2, 2, 2, 2)

    out <- rcpp_calculate_fbd_dosage(
        eHapsCurrent_t = eHapsCurrent_t,
        gamma_t = gamma_t,
        grid = grid,
        snp_start_1_based = 1,
        snp_end_1_based = length(grid)
    )
    expect_equal(length(out$dosage), nSNPs)
    expect_equal(nrow(out$genProbs_t), 3)
    expect_equal(ncol(out$genProbs_t), nSNPs)
    expect_equal(sum(out$dosage == 0), 0)
    expect_equal(sum(out$genProbs_t == 0), 0)

})

test_that("can downsample for gridding appropriately", {

    sampleReads <- list(
        list(0, 0, matrix(0, ncol = 1), matrix(0, ncol = 1)),
        list(0, 1, matrix(0, ncol = 1), matrix(0, ncol = 1)),
        list(0, 1, matrix(0, ncol = 1), matrix(0, ncol = 1)),
        list(0, 2, matrix(0, ncol = 1), matrix(0, ncol = 1))        
    )
    sampleNames <- "jimmy"
    iBam <- 1
    regionName <- "someRegion"
    tempdir <- tempdir()

    for(downsampleToCov in c(1, 2)) {

        out <- downsample_snapped_sampleReads(
            sampleReads,
            iBam,
            downsampleToCov,
            sampleNames,
            verbose = FALSE
        )
        
        expect_equal(
            length(out),
            c(3, 4)[downsampleToCov]
        )
    }

    ## check this works through high level function
    N <- 1
    nCores <- 1
    downsampleToCov <- 1
    inputBundleBlockSize <- NA
    L <- c(1, 2, 3)
    gridWindowSize <- NA
    out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
    grid <- out$grid
    grid_distances <- out$grid_distances
    bundling_info <- get_bundling_position_information(
        N = N,
        nCores = nCores,
        blockSize = inputBundleBlockSize
    )

    ## expect that when it goes it, it has length 3
    expect_equal(
        length(sampleReads),
        4
    )
    save(sampleReads, file = file_sampleReads(tempdir, 1, regionName))

    ## now perform high level function
    snap_reads_to_grid(
        N = N,
        nCores = nCores,
        regionName = regionName,
        tempdir = tempdir,
        bundling_info = bundling_info,
        grid = grid,
        downsampleToCov = downsampleToCov,
        verbose = FALSE
    )
    load(file = file_sampleReads(tempdir, 1, regionName))
    expect_equal(
        length(sampleReads),
        3
    )

})
