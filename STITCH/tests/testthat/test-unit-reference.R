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

