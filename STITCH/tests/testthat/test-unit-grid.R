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
    S <- 3

    for(gridWindowSize in c(NA, 1, 3)) {
        
        out <- assign_positions_to_grid(
            L = L,
            gridWindowSize = gridWindowSize
        )
    
        alphaMatCurrent_tc <- array(0, c(K, out$nGrids, S))
        
        out <- remove_buffer_from_variables(
            L = L,
            regionStart = regionStart,
            regionEnd = regionEnd,
            grid = out$grid,
            grid_distances = out$grid_distances,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
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
                expect_equal(ncol(out$alphaMatCurrent_tc), 4)
            }
        }
        
    }
        

})




## test_that("can use grid, and calculate gp_t", {
## this is deprecated, and used in test-unit-cpp


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
            length(out$sampleReads),
            c(3, 4)[downsampleToCov]
        )
        expect_equal(
            out$remove_stats,
            c(c(1, 0)[downsampleToCov], 4)
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




test_that("can understand genetic map format", {
    
    ## from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
    ## from the top of one of the genetic maps
    genetic_map <- rbind(
        c(150118, 1.13462264157027, 0),
        c(154675, 1.12962782559127, 0.00517047537763574),
        c(154753, 1.13654510133156, 0.00525858634803186),
        c(168567, 1.58657526542862, 0.0209588203778261)
    )
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## add back in last column, check OK
    expect_equal(
        fill_in_genetic_map_cm_column(genetic_map)[, 3],
        genetic_map[, 3]
    )

})
