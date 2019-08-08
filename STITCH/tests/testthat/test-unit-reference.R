test_that("can validate reference iterations", {

    expect_error(validate_reference_iterations("wer"))
    expect_error(validate_reference_iterations("two"))
    expect_error(validate_reference_iterations(NULL))
    expect_error(validate_reference_iterations(NA))
    expect_error(validate_reference_iterations(0.1))
    expect_error(validate_reference_iterations(11.2))    
    
    expect_null(validate_reference_iterations(0))
    expect_null(validate_reference_iterations(1))  
    expect_null(validate_reference_iterations(10))      
})


test_that("reference position SNPs can be valid for niterations = 1", {
    expect_equal(
        validate_pos_and_legend_snps_for_niterations_equals_1(
            legend_snps = c("100-A-C", "101-A-G"),
            pos_snps = c("100-A-C", "101-A-G"),
            niterations = 1
        ),
        NULL
    )
})

test_that("reference position SNPs can be invalid for niterations = 1", {
    expect_error(
        validate_pos_and_legend_snps_for_niterations_equals_1(
            legend_snps = c("100-A-C", "101-A-C"), pos_snps = "100-A-G", niterations = 1
        ),
        "You have selected to use reference haplotypes with niterations=1, which requires exact matching of reference legend SNPs and posfile SNPs. However, reference legend SNP with pos-ref-alt 100-A-C was not found in posfile"
    )
})

test_that("reference position SNPs are ignored for niterations > 1", {
    expect_equal(
        validate_pos_and_legend_snps_for_niterations_equals_1(
            legend_snps = c("100-A-C", "101-A-C"), pos_snps = "100-A-G", niterations = 2
        ),
        NULL
    )
})



test_that("reference sample file requires POP column", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")    ,
        reference_sample_header = c("ID", "POPXX", "GROUP", "SEX")
    )

    sink("/dev/null")
    expect_error(
        get_haplotypes_from_reference(
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_sample_file = refpack$reference_sample_file,
            reference_populations = refpack$reference_populations,
            pos = refpack$pos
        ),
        paste0("Cannot find column POP in reference_sample_file:", refpack$reference_sample_file)
    )
    sink()

})



test_that("reference data can be loaded for an autosome", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")
    )
    sink("/dev/null")
    reference_haps <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = refpack$pos,
        tempdir = tempdir(),
        regionName = "test"
    )
    sink()

    expect_equal(refpack$reference_haplotypes, reference_haps)

})

test_that("haps matrix without NA is passed through", {
    haps <- array(c(0, 1, 0, 1, 0, 1), c(3, 2))
    expect_equal(
        remove_NA_columns_from_haps(haps),
        haps
    )
})


test_that("haps matrix with bad number of NAs is caught", {
    haps <- array(c(0, 1, 0, 1, 0, 1), c(3, 2))
    haps[2, 2] <- NA
    expect_error(
        remove_NA_columns_from_haps(haps),
        "The missing haplotype entry '-' for a reference male sample on chromosome X does not make up the entirety of a sample haplotype for at least one sample for the region of interest"
    )
})


test_that("haps matrix all NA is caught", {
    haps <- array(c(0, 1, 0, 1, 0, 1), c(3, 2))
    haps[, 1] <- NA
    haps[, 2] <- NA
    expect_error(
        remove_NA_columns_from_haps(haps),
        "There are no viable haplotypes from the reference samples for the region of interest"
    )
})

test_that("haps matrix with one column NA is OK and that column gets removed", {
    haps <- array(c(0, 1, 0, 1, 0, 1), c(3, 2))
    haps[, 2] <- NA
    expect_equal(
        remove_NA_columns_from_haps(haps),
        haps[, 1, drop = FALSE]
    )
})



test_that("reference data can be loaded for chromosome X", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        chr = "X"
    )
    reference_haps <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = refpack$pos,
        tempdir = tempdir(),
        regionName = "test"
    )

    sex <- refpack$reference_samples[, "SEX"]
    keep <- array(TRUE, length(sex) * 2)
    keep[2 * which(sex == "male")] <- FALSE
    expect <- matrix(as.integer(refpack$reference_haplotypes[, keep]), ncol = sum(keep))
    expect_equal(
        expect,
        reference_haps
    )

})


test_that("reference data can be loaded for an autosome specifying regionStart, regionEnd and buffer", {

    regionStart <- 5
    regionEnd <- 7
    buffer <- 1
    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")
    )
    sink("/dev/null")
    reference_haps <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = refpack$pos,
        tempdir = tempdir(),
        regionName = "test",
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )
    sink()

    expect <- refpack$reference_haplotypes
    not_NA <- refpack$pos[, "POS"] >= (regionStart - buffer) & refpack$pos[, "POS"] <= (regionEnd + buffer)
    expect[not_NA == FALSE, ] <- NA
    expect_equal(
        expect,
        reference_haps
    )

})


test_that("reference data can be loaded for chromosome X specifying regionStart, regionEnd and buffer", {

    regionStart <- 5
    regionEnd <- 7
    buffer <- 1
    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        chr = "X"
    )
    pos <- refpack$pos
    sink("/dev/null")

    reference_haps <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = pos,
        tempdir = tempdir(),
        regionName = "test",
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    sink()

    expect <- refpack$reference_haplotypes
    not_NA <- pos[, "POS"] >= (regionStart - buffer) & pos[, "POS"] <= (regionEnd + buffer)
    expect[not_NA == FALSE, ] <- NA
    sex <- refpack$reference_samples[, "SEX"]
    keep <- array(TRUE, length(sex) * 2)
    keep[2 * which(sex == "male")] <- FALSE
    expect <- matrix(as.integer(expect[, keep]), ncol = sum(keep))

    expect_equal(
        expect,
        reference_haps
    )

})

test_that("reference data can be loaded for a population", {

    regionStart <- 5
    regionEnd <- 7
    buffer <- 1
    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")
    )
    ref_pops_to_load <- c("GBR", "CHB")
    pos <- refpack$pos
    sink("/dev/null")
    reference_haps <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = ref_pops_to_load,
        pos = pos,
        tempdir = tempdir(),
        regionName = "test",
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )
    sink()

    expect <- refpack$reference_haplotypes
    not_NA <- pos[, "POS"] >= (regionStart - buffer) & pos[, "POS"] <= (regionEnd + buffer)
    expect[not_NA == FALSE, ] <- NA
    keep_samples <- is.na(match(refpack$reference_samples[, "POP"], ref_pops_to_load)) == FALSE
    expect <- expect[, 2 * rep(which(keep_samples), each = 2) + -1:0]

    expect_equal(
        expect,
        reference_haps
    )

})


test_that("can make fake reference samples in C++", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")    ,
        reference_sample_header = c("ID", "POPXX", "GROUP", "SEX")
    )
    reference_haps <- refpack$reference_haplotypes
    non_NA_cols <- which(is.na(reference_haps[ , 1]) == FALSE)
    reference_phred <- 20

    iSample <- 1
    sampleReads1 <- make_sampleReads_from_hap(
        non_NA_cols,
        reference_phred,
        reference_hap = reference_haps[, iSample]
    )
    sampleReads2 <- rcpp_make_sampleReads_from_hap(
        non_NA_cols,
        reference_phred,
        reference_hap = reference_haps[, iSample]
    )
    expect_equal(sampleReads1, sampleReads2)

})


test_that("can do single reference iteration", {

    ## hmm, need to dummy up a bit
    ## should be able to do with N=1
    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 10,
        gridWindowSize = 3,
        S = 2
    )
    S <- test_package$S
    sampleReads <- test_package$sampleReads
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    eMatHapSNP_t <- test_package$list_of_eMatHapSNP_t[[1]]
    sigmaCurrent_m <- test_package$sigmaCurrent_m
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    grid <- test_package$grid
    grid_distances <- test_package$grid_distances

    ## dummy up a bit
    tempdir <- tempfile("tempdir")
    dir.create(tempdir)
    N_haps <- 2
    nCores <- 1
    reference_bundling_info <- NULL
    nGen <- 10000

    regionName <- "jimmy"
    for(iBam in 1:2) {
        save(sampleReads, file = file_referenceSampleReads(tempdir, iBam, regionName))
    }

    out <- single_reference_iteration(eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, sigmaCurrent_m = sigmaCurrent_m, priorCurrent_m = priorCurrent_m, N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, L = L, grid = grid, grid_distances = grid_distances, nGen = nGen)


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


## something like 
## position COMBINED_rate.cM.Mb. Genetic_Map.cM.
## 82590               3.8618       0.0000000
## 82603               3.8740       0.0000502
## 83158               3.8802       0.0022003

## alternatively like
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
