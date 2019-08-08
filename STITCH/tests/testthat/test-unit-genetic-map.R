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
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## add back in last column, check OK
    expect_equal(
        fill_in_genetic_map_cm_column(genetic_map)[, 3],
        genetic_map[, 3]
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


test_that("can simply fill in sigmaCurrent from genetic map", {

    genetic_map_L <- c(100, 4000, 5000, 10000)
    n_snps <- length(genetic_map_L)
    genetic_map <- make_genetic_map_file(genetic_map_L, n_snps, expRate = 0.5)
    
    L <- c(5, 15, 4500, 5000, 5600, 9000, 11000)
    convert_genetic_map_to_sigmaCurrent_m(genetic_map, L, expRate = 0.5, minRate = 0.1, maxRate = 100)
        
})
