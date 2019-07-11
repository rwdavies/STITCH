test_that("can determine switches", {

    nGrids <- 100    
    break_results <- cbind(
        left_grid_break_0_based = c(5, 20, 50),
        left_grid_focal_0_based = c(10, 24, 58),
        right_grid_focal_0_based = c(11, 25, 59),
        right_grid_break_0_based = c(15, 30, 60)
    )
    nbreaks <- nrow(break_results)
    
    K <- 4
    fromMat <- array(0, c(nbreaks, K, K))
    ## now make simple
    for(k in 1:K) {
        fromMat[, k, k] <- 1
    }
    switchOrder <- determine_switch_order(fromMat, nbreaks, K)

    expect_equal(
        switchOrder,
        matrix(rep(1:4, 3), nrow = 3, byrow = TRUE)
    )

})



test_that("can run getBetterSwitchesSimple with / without grid", {

    ## dont worry about s here
    K <- 4
    ##    S <- 3
    nSNPs <- 100
    set.seed(2342)
    L <- sort(sample(10000, 100))
    gridWindowSize <- NA
    
    for(gridWindowSize in c(NA, 500)) {

        out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids
        snps_in_grid_1_based <- out$snps_in_grid_1_based
        
        eHapsFuture_t <- array(runif(K * nSNPs), c(K, nSNPs))
        alphaMatFuture_t <- array(runif(K * (nGrids - 1)), c(K, nGrids - 1))
        
        ## I think these are grids but are they 1 or 0 based
        break_results <- cbind(
            left_grid_break_0_based = c(0, nGrids - 4),
            left_grid_focal_0_based = c(1, nGrids - 3),
            right_grid_focal_0_based = c(2, nGrids - 2),
            right_grid_break_0_based = c(3, nGrids - 1)
        )
        if (nSNPs == nGrids) {
            break_results <- rbind(break_results, c(40, 43, 44, 50))
        }
        break_results <- break_results[order(break_results[, 2]), ]        
        nbreak <- nrow(break_results)

        ## from -> to
        fromMat <- array(0, c(nbreak, K, K))
        ## make wacky!
        for(k in 1:K) {
            fromMat[, k, c(K, 1:(K - 1))[k]] <- 1
        }

        out <- getBetterSwitchesSimple(
            fromMat = fromMat,
            nbreak = nbreak,
            break_results = break_results,
            eHapsFuture_t = eHapsFuture_t,
            alphaMatFuture_t = alphaMatFuture_t,
            grid = grid,
            snps_in_grid_1_based = snps_in_grid_1_based
        )
        test_that(ncol(out$eHapsCurrent_t), nSNPs)
        test_that(ncol(out$alphaMatCurrent_t), nGrids - 1)

    }

})


test_that("can smooth rate", {

    set.seed(10)            
    shuffle_bin_radius <- 5000
    L <- sort(sample(100000, 100))

    for(gridWindowSize in c(NA, 10000)) {
        
        out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids
        
        ## now simulate a rate
        sigmaSum_unnormalized <- runif(nGrids - 1)
        sigma_rate <- -log(sigmaSum_unnormalized) / grid_distances ## by definition
        
        ## just check the same
        results1 <- make_smoothed_rate(
            sigmaSum_unnormalized,
            sigma_rate,
            L_grid,
            grid_distances,
            nGrids,
            shuffle_bin_radius
        )
        
        results2 <- rcpp_make_smoothed_rate(
            sigmaSum_unnormalized,
            sigma_rate,
            L_grid,
            grid_distances,
            nGrids,
            shuffle_bin_radius
        )
        
        expect_equal(as.numeric(results1), as.numeric(results2))

    }
    
})


test_that("can define breaks on tiny region", {

    set.seed(11)
    tempdir <- tempdir()
    regionName <- "blargh"
    L <- sort(sample(1000000, 500))
    shuffle_bin_radius <- 2000
    ## make sure it tests going to the end
    L <- unique(sort(c(1, shuffle_bin_radius, L)))
    L <- unique(sort(c(L, 1000000 - shuffle_bin_radius + 1, 1000000)))
    ## carve out a hole - guarantees, given fixed sigmaSum below, will pass test
    L <- L[L <= (500000 - 20000) | (500000 + 20000) <= L]
    nGen <- 100
    minRate <- 0.1
    maxRate <- 100
    S <- 3

    for(gridWindowSize in c(NA, 10000)) {

        out <- assign_positions_to_grid(L = L, gridWindowSize = gridWindowSize)
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids

        for(i in 1:2) {
        
            if (i == 1) {
                sigmaSum_m_unnormalized <- array(runif(S * (nGrids - 1)), c(nGrids - 1, S))
            } else {
                sigmaSum_m_unnormalized <- array(0.8, c(nGrids - 1, S))
            }
            
            define_and_save_breaks_to_consider(
                tempdir = tempdir,
                regionName = regionName,
                sigmaSum_m_unnormalized = sigmaSum_m_unnormalized,
                L_grid = L_grid,
                grid_distances = grid_distances,
                nGrids = nGrids,
                nGen = nGen,
                minRate = minRate,
                maxRate = maxRate,
                shuffle_bin_radius = 2000,
                plot_shuffle_haplotype_attempts = TRUE
            )

            load(file_break_results(tempdir, regionName))
            ## have an output
            expect_equal(length(break_results) > 0, TRUE)
            
        }

    }
    
})


test_that("can choose points to break if highest on the left", {

    nGrids <- 100
    L_grid <- 1:100
    smoothed_rate <- runif(nGrids)
    smoothed_rate[1:10] <- NA
    smoothed_rate[11] <- 100
    shuffle_bin_radius <- 10
    grid_distances <- diff(L_grid)

    out <- choose_points_to_break(
        smoothed_rate,
        nGrids,
        L_grid,
        shuffle_bin_radius,
        grid_distances
    )

    expect_equal(length(out) > 0, TRUE)    
    
})


test_that("can choose points to break if highest on the right", {

    nGrids <- 100
    L_grid <- 1:100
    grid_distances <- diff(L_grid)    
    smoothed_rate <- runif(nGrids)
    smoothed_rate[90:100] <- NA
    smoothed_rate[89] <- 100
    nGrids <- 100
    shuffle_bin_radius <- 10    
    
    out <- choose_points_to_break(
        smoothed_rate,
        nGrids,
        L_grid,
        shuffle_bin_radius,
        grid_distances
    )
    expect_equal(length(out) > 0, TRUE)

})

test_that("can nuke SNPs, accounting for differences in grid distances", {

    shuffle_bin_radius <- 10
   
    for(x in c(1, 4)) {

        L_grid <- seq(1, 500, x)
        grid_distances <- diff(L_grid)
        nGrids <- length(L_grid)
        
        ## check normal performance
        for(i_where in 1:3) {
            if (i_where == 1) {
                snp_best <- 100
            } else if (i_where == 2) {
                snp_best <- 2
            } else if (i_where == 3) {
                snp_best <- nGrids - 2
            }
            snp_left <- snp_best - 1
            snp_right <- snp_best + 1
            out <- get_snps_to_nuke(grid_distances, nGrids, shuffle_bin_radius, snp_best, snp_left, snp_right)
            e <- range(which(abs(mean(L_grid[snp_best + 0:1]) - L_grid) < (2 * shuffle_bin_radius)))
            e[2] <- e[2] - 1
            expect_equal(out[["nuke_left"]], e[1])
            expect_equal(out[["nuke_right"]], e[2])
        }

    }
        
})



test_that("can smooth rate with lots of small peaks", {

    ## from some human data
    smoothed_rate <- c(NA, 0.015442, 0.008321, 0.017876, 0.041003, 0.059737, 0.064171, 0.048742, 0.026267, 0.011444, 0.011549, 0.013898, 0.009484, 0.00977, 0.010118, 0.026074, 0.044993, 0.033122, 0.025366, 0.028651, 0.023531, 0.017467, 0.015018, 0.02795, 0.038438, 0.028658, 0.024437, 0.027357, 0.022311, 0.01097, 0.005014, 0.007553, 0.007554, 0.003408, 0.001296, 0.002222, 0.005499, 0.00863, 0.006101, 0.001935, 0.000806, 0.000806, 0.000785, 0.000804, 0.006186, 0.01982, 0.023825, 0.011826, 0.003326, 0.004388, 0.01227, 0.02641, 0.036658, 0.039786, 0.033915, 0.01731, 0.010189, 0.011444, 0.006264, 0.003138, 0.007011, 0.009902, 0.006055, 0.001673, 0.000822, 0.001307, 0.002539, 0.003097, 0.002452, 0.003193, 0.018705, 0.032085, 0.016704, 0.006522, 0.010918, 0.009187, 0.01048, 0.02064, 0.035764, 0.037283, 0.02221, 0.010872, 0.005437, 0.004629, 0.010243, 0.013006, 0.012254, 0.015208, 0.015867, 0.012694, 0.009785, 0.006441, 0.005936, 0.009198, 0.011125, 0.014293, 0.016029, 0.010694, 0.007276, 0.006402, 0.015751, 0.02593, 0.015458, 0.006125, 0.005897, 0.003122, 0.002552, 0.013242, 0.025577, 0.020003, 0.00872, 0.005058, 0.003864, 0.002985, 0.002579, 0.003486, 0.017944, 0.031849, 0.018007, 0.003711, 0.003711, 0.003425, 0.002489, 0.001346, 0.000684, 0.000898, 0.007784, 0.017451, 0.015095, 0.007295, 0.003537, 0.001633, 0.001052, 0.001066, 0.001163, 0.002309, 0.020467, 0.040819, 0.032187, 0.016, 0.006832, 0.002023, 0.00146, 0.002134, 0.003296, 0.006852, 0.010624, 0.006144, 0.002228, 0.004989, 0.01136, 0.015713, 0.011125, 0.005139, 0.003953, 0.003086, 0.001431, 0.006951, 0.018598, 0.027343, 0.025602, 0.011728, 0.00367, 0.006208, 0.007979, 0.011988, 0.021989, 0.02442, 0.014042, 0.011465, 0.021787, 0.024689, 0.017985, 0.023467, 0.032068, 0.011743, 0.026588, 0.022205, 0.013469, 0.015223, 0.013626, 0.005207, 0.002165, 0.003318, 0.004675, 0.015247, 0.026167, 0.019929, 0.018465, 0.034921, 0.03637, 0.015298, 0.002891, 0.009799, 0.027004, 0.03738, 0.037052, 0.048004, 0.052004, 0.031631, 0.021818, 0.018155, 0.008239, 0.010048, 0.014152, 0.008168, 0.010741, 0.018988, 0.018029, 0.017533, 0.018979, 0.018765, 0.01487, 0.01021, 0.007762, 0.005358, 0.002864, 0.002395, 0.003228, 0.00316, 0.010024, 0.018845, 0.015556, 0.010033, 0.007995, 0.003896, 0.001308, 0.001653, 0.00432, 0.007042, 0.005912, 0.003481, 0.007436, NA)
    L_grid <- c(902500, 907500, 912500, 917500, 922500, 927500, 932500, 937500, 942500, 947500, 952500, 957500, 962500, 967500, 972500, 977500, 982500, 987500, 992500, 997500, 1002500, 1007500, 1012500, 1017500, 1022500, 1027500, 1032500, 1037500, 1042500, 1047500, 1052500, 1057500, 1062500, 1067500, 1072500, 1077500, 1082500, 1087500, 1092500, 1097500, 1102500, 1107500, 1112500, 1117500, 1122500, 1127500, 1132500, 1137500, 1142500, 1147500, 1152500, 1157500, 1162500, 1167500, 1172500, 1177500, 1182500, 1187500, 1192500, 1197500, 1202500, 1207500, 1212500, 1217500, 1222500, 1227500, 1232500, 1237500, 1242500, 1247500, 1252500, 1257500, 1262500, 1267500, 1272500, 1277500, 1282500, 1287500, 1292500, 1297500, 1302500, 1307500, 1312500, 1317500, 1322500, 1327500, 1332500, 1337500, 1342500, 1347500, 1352500, 1357500, 1362500, 1367500, 1372500, 1377500, 1382500, 1387500, 1392500, 1397500, 1402500, 1407500, 1412500, 1417500, 1422500, 1427500, 1432500, 1437500, 1442500, 1447500, 1452500, 1457500, 1462500, 1467500, 1472500, 1477500, 1482500, 1487500, 1492500, 1497500, 1502500, 1507500, 1512500, 1517500, 1522500, 1527500, 1532500, 1537500, 1542500, 1547500, 1552500, 1557500, 1562500, 1587500, 1592500, 1597500, 1602500, 1607500, 1612500, 1617500, 1622500, 1627500, 1632500, 1637500, 1642500, 1647500, 1652500, 1657500, 1662500, 1667500, 1672500, 1677500, 1682500, 1687500, 1692500, 1697500, 1702500, 1707500, 1712500, 1717500, 1722500, 1727500, 1732500, 1737500, 1742500, 1747500, 1752500, 1757500, 1762500, 1767500, 1772500, 1777500, 1782500, 1787500, 1792500, 1797500, 1807500, 1812500, 1817500, 1822500, 1827500, 1832500, 1837500, 1842500, 1847500, 1852500, 1857500, 1862500, 1867500, 1872500, 1877500, 1882500, 1887500, 1892500, 1897500, 1902500, 1907500, 1912500, 1917500, 1922500, 1927500, 1932500, 1937500, 1942500, 1947500, 1952500, 1957500, 1962500, 1967500, 1972500, 1977500, 1982500, 1987500, 1992500, 1997500, 2002500, 2007500, 2012500, 2017500, 2022500, 2027500, 2032500, 2037500, 2042500, 2047500, 2052500, 2057500, 2062500, 2067500, 2072500, 2077500, 2082500, 2087500, 2092500, 2097500)

    nGrids <- 235
    grid_distances <- diff(L_grid)
    shuffle_bin_radius <- 5000
    out <- choose_points_to_break(
        smoothed_rate = smoothed_rate,
        nGrids = nGrids,
        L_grid = L_grid,
        shuffle_bin_radius = shuffle_bin_radius,
        grid_distances = grid_distances
    )
    ## there are lots of little breaks here. try to capture them all / many
    results <- out$results

    expect_equal(nrow(results), 22) ## not really smart, but can manually interrogate
    
})


test_that("can smooth rate with lots of big peaks", {

    ##w <- L_grid >= 2000000 & L_grid <= 3000000
    ##smoothed_rate <- smoothed_rate[w]
    ##L_grid <- L_grid[w]
    ## save(smoothed_rate, L_grid, file = "STITCH/tests/testthat/data-shuffle.RData")
    
    load("data-shuffle.RData") ## load up data
    grid_distances <- diff(L_grid)
    nGrids <- length(L_grid)
    shuffle_bin_radius <- 5000
    out <- choose_points_to_break(
        smoothed_rate = smoothed_rate,
        nGrids = nGrids,
        L_grid = L_grid,
        shuffle_bin_radius = shuffle_bin_radius,
        grid_distances = grid_distances
    )
    ## there are lots of little breaks here. try to capture them all / many
    results <- out$results

    expect_equal(nrow(results), 17) ## can manually interrogate it if breaks, determine if appropriate

})
