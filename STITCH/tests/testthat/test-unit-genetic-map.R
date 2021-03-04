test_that("can understand genetic map format", {

    ## from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
    ## from the top of one of the genetic maps
    genetic_map <- rbind(
        c(150118, 1.13462264157027, 0),
        c(154675, 1.12962782559127, 0.00517047537763574),
        c(154753, 1.13654510133156, 0.00525858634803186),
        c(168567, 1.58657526542862, 0.0209588203778261)
    )
    genetic_map[4, 2] <- 0 ## pretend entire file
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## add back in last column, check OK
    expect_equal(
        fill_in_genetic_map_cm_column(genetic_map)[, 3],
        genetic_map[, 3]
    )
    expect_equal(
        fill_in_genetic_map_rate_column(genetic_map)[, 2],
        genetic_map[, 2]
    )

    ## check alt form - argh, needs stuff not here, just remove
    ## genetic_map[, 2] <- 1
    ## genetic_map[, 3] <- NA
    ## expect_equal(
    ##     fill_in_genetic_map_cm_column(genetic_map)[, 3],
    ##     c(0, cumsum(diff(genetic_map[, 1]) * expRate)) / 1e6
    ## )
    

})


test_that("can understand genetic map format second file", {

    ## Position(bp)    Rate(cM/Mb)     Map(cM) Filtered
    genetic_map <- rbind(
        c(63231,       3.86179280588,   0.0),
        c(63244,       3.87400693386,   5.02025586727e-05),
        c(63799,       3.88023863181,   0.00220027640697),
        c(68749,       7.72268441237,   0.0214074585932),
        c(69094,       7.81094442978,   0.0240717836608),
        c(71093,       0.443658545019,  0.0396858627648),
        c(74347,       0.443390102646,  0.0411295266732)
    )
    genetic_map[nrow(genetic_map), 2] <- 0
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## add back in last column, check OK
    expect_equal(
        fill_in_genetic_map_cm_column(genetic_map)[, 3],
        genetic_map[, 3]
    )
    ## weird - very slightly off
    expect_equal(
        max(fill_in_genetic_map_rate_column(genetic_map)[, 2] -         genetic_map[, 2]) < 1e-4,
        TRUE
    )

})



## file looks like
## position COMBINED_rate.cM.Mb. Genetic_Map.cM.
## 82590               3.8618       0.0000000
## 82603               3.8740       0.0000502
## 83158               3.8802       0.0022003

## or alternatively like
##position COMBINED_rate(cM/Mb) Genetic_Map(cM)
## 150118 1.13462264157027 0
## 154675 1.12962782559127 0.00517047537763574
test_that("can load and validate reference genetic map", {

    refpack <- make_reference_package()
    reference_genetic_map_file <- refpack$reference_genetic_map_file
    genetic_map <- read.table(reference_genetic_map_file, header = TRUE)
    expect_null(
        validate_genetic_map(genetic_map)
    )

})

test_that("can error invalid genetic reference map", {

    L <- 1:10
    n_snps <- 10
    genetic_map <- make_genetic_map_file(L, n_snps, expRate = 0.5)
    genetic_map[5, "Genetic_Map.cM."] <- 2 * genetic_map[5, "Genetic_Map.cM."]
    expect_error(validate_genetic_map(genetic_map, verbose = FALSE))

})


test_that("can simply match genetic map to desired positions", {

    ## make a (fake) truth per-bp genetic map
    ## that fully matches what we've sampled
    genetic_map_L <- c(101, 4001, 5001, 10001)
    mult <- 1000 ## for debugging - make differences larger
    expRate <- mult * 0.8 ## cM / Mb
    truth_rate <- rep(expRate, 20000)
    truth_rate[101:4000] <- 0.7 * mult
    truth_rate[4001:5000] <- 0.6 * mult
    truth_rate[5001:10000] <- 0.5 * mult
    truth_map <- c(0, cumsum(truth_rate / 1e6))

    ## now, using the truth genetic map, sample some points
    genetic_map <- array(0, c(length(genetic_map_L), 3))
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    genetic_map[, "position"] <- genetic_map_L
    genetic_map[, "Genetic_Map.cM."] <- truth_map[genetic_map_L]
    genetic_map <- fill_in_genetic_map_rate_column(genetic_map)

    ## now, have the points we want
    L <- c(5, 15, 4500, 5000, 5600, 9000, 11000)
    cM <- match_genetic_map_to_L(genetic_map, L, expRate = expRate)
    ## now get truth
    truth_genetic_map <- array(0, c(length(L), 3))
    colnames(truth_genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    truth_genetic_map[, "position"] <- L
    truth_genetic_map[, "Genetic_Map.cM."] <- truth_map[L]
    truth_genetic_map <- fill_in_genetic_map_rate_column(truth_genetic_map)    

    difference <- truth_genetic_map[, "Genetic_Map.cM."] - cM
    expect_lt(max(difference), 1e-4)

})


test_that("can initialize sigmaCurrent the same way with constant rate", {

    nGen <- 100
    expRate <- 0.79
    L <- sort(sample(10000, 100))
    truth_rate <- rep(expRate, 10000)
    truth_map <- c(0, cumsum(truth_rate / 1e6))
    dl <- diff(L)
    S <- 1
    nGrids <- length(L)
    cM <- truth_map[L]
    cM_grid <- cM
    
    original <- initialize_sigmaCurrent_m(cM_grid = NULL, nGen = nGen, nGrids = nGrids, S = S, dl = dl, expRate = expRate)
    new <- initialize_sigmaCurrent_m(cM_grid = cM_grid, nGen = nGen, nGrids = nGrids, S = S, dl = dl, expRate = expRate)    
    expect_equal(original, new)
    

})


test_that("can assign grid positions when genetic map is given", {

    ## so make a genetic map with a few peaks
    truth_rate <- array(0, 10000)
    ## make a sharp peak around 3000
    truth_rate[3001 + -500:500] <- 10 * dnorm(x = seq(-1, 1, length.out = 1001), sd = 0.1)
    ## make a diffuse peak around 6000
    truth_rate[6001 + -1000:1000] <- 10 * dnorm(x = seq(-1, 1, length.out = 2001), sd = 0.25)
    truth_rate <- truth_rate + runif(length(truth_rate))
    truth_map <- c(0, cumsum(truth_rate / 1e6))
    ##
    ## want grid to focus around those regions
    L <- seq(101, 9901, 25)
    L_grid <- L
    ## L <- c(223, 334, 503, 723, 757, 823, 1061, 1145, 1247, 1285, 1331, 1430, 1624, 1631, 1752, 1808, 1817, 1867, 1976, 1990, 2329, 2347, 2469, 2478, 2496, 2514, 2691, 2779, 2842, 3052, 3218, 3326, 3386, 3430, 3696, 3727, 3742, 3786, 4090, 4266, 4318, 4483, 4652, 4756, 5067, 5074, 5094, 5108, 5238, 5278, 5287, 5305, 5483, 5507, 5614, 5674, 5704, 5850, 5912, 5985, 5998, 6120, 6395, 6468, 6515, 6682, 6711, 6772, 6988, 6990, 7044, 7047, 7064, 7087, 7175, 7216, 7301, 7545, 7886, 8031, 8093, 8173, 8192, 8208, 8480, 8520, 8539, 8648, 8662, 8797, 9071, 9139, 9140, 9141, 9353, 9410, 9612, 9641, 9732, 9991)
    cM <- truth_map[L]
    ##
    ##
    gridWindowSize <- 1000
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize,
        cM = cM
    )

    ## check the peak areas have proportionately more grids
    expect_gt(
        length(unique(out$grid[2501 <= L & L <= 3501])) / 1000,
        out$nGrids / 10000
    )
    expect_gt(
        length(unique(out$grid[5001 <= L & L <= 7000])) / 1000,
        out$nGrids / 10000
    )

    smoothed_cM <- make_smoothed_cM_rate(
        cM = cM,
        L_grid = L_grid,
        shuffle_bin_radius = 250
    )
    
    ## optional - plot to verify manually
    if (1 == 0) {

        par(mfrow = c(2, 1))
        ##
        ## rate
        ##
        plot(truth_rate, type = "l")
        ## also - plot smoothed
        segments(
            x0 = L[-length(L)],
            x1 = L[-1] - 1,
            y0 = diff(smoothed_cM) / diff(L) * 1e6,
            y1 = diff(smoothed_cM) / diff(L) * 1e6,
            col = "orange",
            lwd = 3
        )
        
        ##
        ## cumsum
        ##
        y <- cumsum(truth_rate / 1e6)
        ## plot rates at SNPs
        plot(x = L, y = cM, col = c("orange", "green")[(out$grid %% 2) + 1], type = "o")
        ylim <- range(y)
        ## grids
        snps_in_grid_1_based <- out$snps_in_grid_1_based
        for(iRow in 1:nrow(snps_in_grid_1_based)) {
            s <- snps_in_grid_1_based[iRow, 1]
            e <- snps_in_grid_1_based[iRow, 2]
            y0 <- ylim[1] + diff(ylim) * (0.25 + 0.5 * iRow %% 2)
            rect(
                xleft = L[s], xright = L[e],
                ybottom = ylim[1] + 0.25 * diff(ylim),
                ytop = ylim[1] + 0.75 * diff(ylim)
            )
        }
        points(x = out$L_grid, y = out$cM_grid, col = "blue", type = "o")

    }
    
})


