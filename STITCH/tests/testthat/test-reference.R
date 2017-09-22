test_that("reference legend file has expected columns of id, position, a0 and a1", {

    good_reference_legend_header <- c("id", "position", "a0", "a1")

    ## still do not require ID to exist
    ## require uniqueness on position instead
    sink("/dev/null")    
    for (i_bad_header in 2:4) {
        
        reference_legend_header <- good_reference_legend_header
        reference_legend_header[i_bad_header] <- "something_else"
        refpack <- make_reference_package(
            n_snps = 10,
            n_samples_per_pop = 4,
            reference_populations = c("CEU", "GBR", "CHB"),
            reference_legend_header = reference_legend_header
        )

        expect_error(
            get_haplotypes_from_reference(
                reference_haplotype_file = refpack$reference_haplotype_file,
                reference_legend_file = refpack$reference_legend_file,
                reference_populations = refpack$reference_populations,
                pos = refpack$pos
            ),
            paste0(
                "Cannot find '",
                good_reference_legend_header[i_bad_header],
                "' column ",
                "in reference_legend_file"
            )
        )
    }
    sink()
    
})


test_that("reference legend position is required to be an integer", {

    n_snps <- 10
    L <- c(1:5, "A", 7:10)
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        refs = rep("A", n_snps),
        alts = rep("G", n_snps),
        L = L
    )

    expect_error(
        get_haplotypes_from_reference(
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_populations = refpack$reference_populations,
            pos = refpack$pos
        ),
        paste0("reference_legend_file position needs to be integer valued positions with values between 1 and ", .Machine$integer.max, " but in row 6 A was observed")
    )
    

})

test_that("reference legend position is required to be sorted", {

    n_snps <- 10
    L <- c(10, 1:9)
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        refs = rep("A", n_snps),
        alts = rep("G", n_snps),
        L = L
    )
    
    expect_error(
        get_haplotypes_from_reference(
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_populations = refpack$reference_populations,
            pos = refpack$pos
        ),
        "reference_legend_file position needs to be sorted in strictly increasing order but row number 1 has position 10 and row number 2 has position 1. Please further note that STITCH uses IMPUTE2 style reference legend and haplotype files that have to be split by chromosome, and that no checking of correct chromosome is done programatically by STITCH as chromosome is not a column in the reference legend file, so please match your files with care. Sorry for the inconvenience"
    )
    
})


test_that("reference legend position repeated more than once is OK if different alleles", {

    n_snps <- 10
    L <- c(1:5, 5:9)
    refs <- rep("A", n_snps)
    alts <- rep("G", n_snps)
    alts[6] <- "T"
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        refs = refs,
        alts = alts,
        L = L
    )
    
    expect_equal(
        get_haplotypes_from_reference(
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_populations = refpack$reference_populations,
            pos = refpack$pos
        ),
        refpack$reference_haplotypes
    )
    
})


test_that("reference legend file with duplicate entries throws an error", {
    reference_legend <- data.frame(
        position = c(10, 10, 20),
        a0 = c("A", "A", "T"),
        a1 = c("A", "A", "G")
    )
    expect_error(
        validate_reference_legend(reference_legend),
        "There are 1 duplicate row ids. One such example is 10 A A"
    )
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

