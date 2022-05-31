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




test_that("sensible error message thrown when asked to load reference out of range", {

    refpack <- make_reference_package(
        n_snps = 10,
        L = 11:20,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")    ,
        reference_sample_header = c("ID", "POPXX", "GROUP", "SEX")
    )


    reference_legend_file <- refpack$reference_legend_file
    regionStart <- 30
    regionEnd <- 40
    buffer <- 3
    
    legend_header <- as.character(unlist(read.table(reference_legend_file, nrow = 1, sep = " ")))
    validate_legend_header(legend_header)

    ## no error for >=1 variant loaded
    out <- load_reference_legend(
        legend_header = legend_header,
        reference_legend_file = reference_legend_file,
        regionStart = 14,
        regionEnd = 15,
        buffer = 6
    )
    expect_true(length(out) > 0)


    ## error for out of range
    expect_error(
        load_reference_legend(
            legend_header = legend_header,
            reference_legend_file = reference_legend_file,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer
        )
    )


})







test_that("reference sample file requires POP column", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")    ,
        reference_sample_header = c("ID", "POPXX", "GROUP", "SEX")
    )

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

})



test_that("reference data can be loaded for an autosome", {

    refpack <- make_reference_package(
        n_snps = 10,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB")
    )
    out <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = refpack$pos,
        tempdir = tempdir(),
        regionName = "test"
    )
    
    reference_haps <- out[["reference_haps"]]
    expect_equal(refpack[["reference_haplotypes"]], reference_haps)
    expect_equal(out[["rhb1"]], out[["rhb2"]])

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
    out <- get_haplotypes_from_reference(
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations,
        pos = refpack$pos,
        tempdir = tempdir(),
        regionName = "test"
    )
    reference_haps <- out[["reference_haps"]]

    sex <- refpack$reference_samples[, "SEX"]
    keep <- array(TRUE, length(sex) * 2)
    keep[2 * which(sex == "male")] <- FALSE
    expect <- matrix(as.integer(refpack$reference_haplotypes[, keep]), ncol = sum(keep))
    expect_equal(expect, reference_haps)
    expect_equal(out[["rhb1"]], out[["rhb2"]])
    expect_equal(out[["rhb1"]], out[["rhb3"]])

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
    out <- get_haplotypes_from_reference(
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
    reference_haps <- out[["reference_haps"]]

    ##
    not_NA <- refpack$pos[, "POS"] >= (regionStart - buffer) & refpack$pos[, "POS"] <= (regionEnd + buffer)
    expect <- refpack$reference_haplotypes[not_NA, ]
    expect_equal(expect, reference_haps)
    expect_equal(out[["rhb1"]], out[["rhb2"]])
    expect_equal(out[["rhb1"]], out[["rhb3"]])
    
    
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

    out <- get_haplotypes_from_reference(
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
    reference_haps <- out[["reference_haps"]]

    not_NA <- pos[, "POS"] >= (regionStart - buffer) & pos[, "POS"] <= (regionEnd + buffer)    
    expect <- refpack$reference_haplotypes[not_NA, ]
    sex <- refpack$reference_samples[, "SEX"]
    keep <- array(TRUE, length(sex) * 2)
    keep[2 * which(sex == "male")] <- FALSE
    expect <- matrix(as.integer(expect[, keep]), ncol = sum(keep))

    expect_equal(expect, reference_haps)
    expect_equal(out[["rhb1"]], out[["rhb2"]])
    expect_equal(out[["rhb1"]], out[["rhb3"]])

    
})

test_that("reference data can be loaded for a population", {

    ## specifically ensure either 31, 32 or 33 SNPS are to be loaded
    regionStart <- 8
    buffer <- 2

    x <- ((regionStart) - (2 * buffer + 1))
    test_lengths <- c(c(31, 32, 33), 32 + c(31, 32, 33))
    for(regionEnd in (test_lengths + x)) {

        refpack <- make_reference_package(
            n_snps = regionEnd + buffer + 5,
            n_samples_per_pop = 4,
            reference_populations = c("CEU", "GBR", "CHB")
        )
        ref_pops_to_load <- c("GBR", "CHB")
        pos <- refpack$pos
        out <- get_haplotypes_from_reference(
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
        reference_haps <- out[["reference_haps"]]
        
        not_NA <- pos[, "POS"] >= (regionStart - buffer) & pos[, "POS"] <= (regionEnd + buffer)    
        expect <- refpack$reference_haplotypes[not_NA, ]
        keep_samples <- is.na(match(refpack$reference_samples[, "POP"], ref_pops_to_load)) == FALSE
        expect <- expect[, 2 * rep(which(keep_samples), each = 2) + -1:0]
        
        expect_equal(expect, reference_haps)
        expect_equal(out[["rhb1"]], out[["rhb2"]])
        expect_equal(out[["rhb1"]], out[["rhb3"]])
    }
    
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


test_that("can convert a reference hap into eMatGrid_t directly", {

    for(has_NA in c(FALSE, TRUE)) {
        for(gridWindowSize in c(NA, 3)) {
        
            test_package <- make_fb_test_package(
                K = 4,
                nReads = 8,
                nSNPs = 10,
                gridWindowSize = gridWindowSize,
                S = 2
            )
            nSNPs <- test_package$nSNPs
            nGrids <- test_package$nGrids
            K <- test_package$K
            S <- test_package$S
            
            reference_haps <- t(round(test_package$eHapsCurrent_tc[, , 1]))
            reference_haps_full <- reference_haps
            if (has_NA) {
                reference_haps <- reference_haps[-c(4, 7), ]
                rh_in_L <- (1:nSNPs)[-c(4, 7)]
            } else {
                rh_in_L <- 1:nrow(reference_haps)
            }

            grid <- test_package$grid
            eHapsCurrent_tc <- test_package$eHapsCurrent_tc
            
            eMatGrid_t1 <- array(1, c(K, nGrids))
            eMatGrid_t2 <- array(1, c(K, nGrids))        
            s <- 0
            iSample <- 0
            reference_phred <- 20
            maxEmissionMatrixDifference <- 1e10

            ehh_h1_S <- array(0, c(K, nSNPs, S))
            ehh_h1_D <- array(0, c(K, nSNPs, S))
            ehh_h0_S <- array(0, c(K, nSNPs, S))
            ehh_h0_D <- array(0, c(K, nSNPs, S))
            
            ref_make_ehh(
                eHapsCurrent_tc = eHapsCurrent_tc,
                ehh_h1_S = ehh_h1_S,
                ehh_h1_D = ehh_h1_D,
                ehh_h0_S = ehh_h0_S,
                ehh_h0_D = ehh_h0_D,
                reference_phred = reference_phred
            )

            reference_hap <- reference_haps[, iSample + 1]
            rcpp_ref_make_eMatGrid_t(
                eMatGrid_t = eMatGrid_t1,
                reference_hap = reference_hap,
                rh_in_L = rh_in_L,
                eHapsCurrent_tc = eHapsCurrent_tc,
                grid = grid,
                reference_phred = reference_phred,
                s = s,
                maxEmissionMatrixDifference = maxEmissionMatrixDifference,
                ehh_h1_S = ehh_h1_S,
                ehh_h0_S = ehh_h0_S,
                rescale = TRUE,
                bound = TRUE
            )

            ## compare against old way
            sampleReads <- rcpp_make_sampleReads_from_hap(
                rh_in_L = rh_in_L,
                reference_phred = reference_phred,
                reference_hap = reference_haps_full[, iSample + 1] ## old method, requires full
            )
            ## 
            sampleReads <- snap_sampleReads_to_grid(
                sampleReads = sampleReads,
                grid = grid
            )
            ## 
            eMatRead_t <- array(1, c(K, length(sampleReads)))
            rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads = 1e10, Jmax = 100, eMatHapOri_t = array(0, c(1, 1)), pRgivenH1 = array(0, 1), pRgivenH2 = array(0, 1), prev = 0, suppressOutput = 1, prev_section ="", next_section = "", run_pseudo_haploid = FALSE);
            ## 
            rcpp_make_eMatGrid_t(eMatGrid_t2, eMatRead_t, 1, sampleReads, 1, nGrids, prev = 0, suppressOutput = 1, prev_section = "",next_section = "", run_fb_grid_offset = 0, TRUE, TRUE, maxEmissionMatrixDifference = maxEmissionMatrixDifference, rescale = TRUE)

            expect_equal(eMatGrid_t1, eMatGrid_t2)
            
        }

    }
    
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

    ## close enough!
    reference_haps <- t(round(test_package$eHapsCurrent_tc[, , 1]))
    rhi <- reference_haps
    rhi_t <- t(rhi)
    rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
    rhb <- t(rhb_t)    
    rh_in_L <- which(is.na(reference_haps[ , 1]) == FALSE)

    
    ## dummy up a bit
    tempdir <- tempfile("tempdir")
    dir.create(tempdir)
    N_haps <- 2
    nCores <- 1
    nGen <- 10000

    out <- single_reference_iteration(eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, sigmaCurrent_m = sigmaCurrent_m, priorCurrent_m = priorCurrent_m, N_haps = N_haps, nCores = nCores, tempdir = tempdir, regionName = regionName, L = L, grid = grid, grid_distances = grid_distances, nGen = nGen, rhb = rhb, rh_in_L = rh_in_L)

    ## placeholder
    expect_equal(class(out), "list")

})


test_that("can make reference ematgrid with large grid", {

    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 1000,
        gridWindowSize = 1000,
        S = 2
    )
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    S <- test_package$S
    grid <- test_package$grid
    reference_haps <- t(round(test_package$eHapsCurrent_tc[, , 1]))
    reference_haps_full <- reference_haps
    rh_in_L <- 1:nrow(reference_haps)
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
            
    eMatGrid_t1 <- array(1, c(K, nGrids))
    eMatGrid_t2 <- array(1, c(K, nGrids))        
    s <- 0
    iSample <- 0
    reference_phred <- 20
    maxEmissionMatrixDifference <- 1e10
    
    ehh_h1_S <- array(0, c(K, nSNPs, S))
    ehh_h1_D <- array(0, c(K, nSNPs, S))
    ehh_h0_S <- array(0, c(K, nSNPs, S))
    ehh_h0_D <- array(0, c(K, nSNPs, S))
    
    ref_make_ehh(
        eHapsCurrent_tc = eHapsCurrent_tc,
        ehh_h1_S = ehh_h1_S,
        ehh_h1_D = ehh_h1_D,
        ehh_h0_S = ehh_h0_S,
        ehh_h0_D = ehh_h0_D,
        reference_phred = reference_phred
    )

    ## take bad hap
    reference_hap <- sample(c(0, 1), nSNPs, replace = TRUE)

    eMatGrid_t <- rcpp_ref_make_eMatGrid_t(
        eMatGrid_t = eMatGrid_t1,
        reference_hap = reference_hap,
        rh_in_L = rh_in_L,
        eHapsCurrent_tc = eHapsCurrent_tc,
        grid = grid,
        reference_phred = reference_phred,
        s = s,
        maxEmissionMatrixDifference = maxEmissionMatrixDifference,
        ehh_h1_S = ehh_h1_S,
        ehh_h0_S = ehh_h0_S,
        rescale = TRUE,
        bound = TRUE
    )
    expect_equal(sum(is.na(eMatGrid_t1)), 0)
    
    
})
