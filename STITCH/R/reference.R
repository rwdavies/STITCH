get_and_initialize_from_reference <- function(
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    hapSumCurrent_tc,
    sigmaCurrent_m,
    priorCurrent_m,
    reference_haplotype_file,
    reference_legend_file,
    reference_sample_file,
    reference_populations,
    reference_phred,
    reference_iterations,
    nSNPs,
    K,
    S,
    L,
    pos,
    inputBundleBlockSize,
    nCores,
    regionName,
    alleleCount,
    expRate,
    nGen,
    tempdir,
    outputdir,
    pseudoHaploidModel,
    emissionThreshold,
    alphaMatThreshold,
    minRate,
    maxRate,
    regionStart,
    regionEnd,
    buffer,
    niterations,
    grid,
    grid_distances,
    nGrids,
    reference_shuffleHaplotypeIterations,
    L_grid,
    plot_shuffle_haplotype_attempts,
    shuffle_bin_radius,
    snps_in_grid_1_based,
    plotHapSumDuringIterations
) {

    print_message("Begin initializing paramters using reference haplotypes")
    ## get reference haplotypes matched to posfile
    ## NA's where there are no match
    ## only for populations of interest
    out <- get_haplotypes_from_reference(
        reference_haplotype_file = reference_haplotype_file,
        reference_legend_file = reference_legend_file,
        reference_sample_file = reference_sample_file,
        reference_populations = reference_populations,
        pos = pos,
        tempdir = tempdir,
        regionName = regionName,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        chr = chr,
        niterations = niterations,
        extraction_method = "hap_v3" ## to do - make both, then only the one I want
    )
    ## change here - make rh compact

    ## rhb is the old "reference_haps" in binary form, where "SNPs" are 32-compressed
    ##   consecutive SNPs, stored in integer form
    ## reference_haps is a row = SNP, column = sample/K, where 0 = ref, 1 = alt
    ## defined at reference haplotype SNPs, a subset of pos
    ## rh_in_L is 1-based mapping of where reference_haplotypes are in L
    ## i.e. rh_in_L being 1, 3, 5 means that reference_haplotypes has 3 rows,
    ##     mapping to SNPs 1, 3 and 5 in pos (with 2 and 4 not being in the reference)
    ##
    rhb <- out[["rhb3"]]
    rh_in_L <- out[["rh_in_L"]]
    ref_alleleCount <- out[["ref_alleleCount3"]] ## defined at all SNPs

    reference_panel_SNPs <- array(FALSE, nSNPs)
    reference_panel_SNPs[rh_in_L] <- TRUE
    ## note - non_NA_cols is deprecated, being the same thing
    ## non_NA_cols <- which(is.na(reference_haps[ , 1]) == FALSE)

    if (is.null(alleleCount) == FALSE) {
        compare_reference_haps_against_alleleCount(
            alleleCount = alleleCount,
            ref_alleleCount = ref_alleleCount,
            outputdir = outputdir,
            regionName = regionName
        )
    }

    N_haps <- ncol(rhb)
    nRefSNPs <- length(rh_in_L)
    e <- emissionThreshold
    rhb_t <- t(rhb)    
    
    if (K > N_haps) {

        print_message("K is set to be more than the number of reference haplotypes. The rest of the K ancestral haplotypes will be filled with noise to start")

        ## hmm, not sure how efficient
        rm(rhb); if (N_haps > 1000) { gc(reset = TRUE); gc(reset = TRUE);}        
        eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1] <- inflate_fhb_t(rhb_t = rhb_t, haps_to_get = 0:(N_haps - 1), nSNPs = nRefSNPs)
        eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1][eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1] == 0] <- e
        eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1][eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1] == 1] <- (1 - e)

        ## better probably
        if (S > 1) {
            for(s in 2:S) {
                eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, s] <- eHapsCurrent_tc[1:N_haps, reference_panel_SNPs, 1]
            }
        }
        

    } else if (K == N_haps) {

        print_message("There are exactly as many reference haplotypes as K. Using these haplotypes directly as the initial estimate of the ancestral haplotypes")
        
        rm(rhb); if (N_haps > 1000) { gc(reset = TRUE); gc(reset = TRUE);}
        eHapsCurrent_tc[, reference_panel_SNPs, 1] <- inflate_fhb_t(rhb_t = rhb_t, haps_to_get = 0:(N_haps - 1), nSNPs = nRefSNPs)
        eHapsCurrent_tc[, reference_panel_SNPs, 1][eHapsCurrent_tc[, reference_panel_SNPs, 1] == 0] <- e
        eHapsCurrent_tc[, reference_panel_SNPs, 1][eHapsCurrent_tc[, reference_panel_SNPs, 1] == 1] <- (1 - e)
        if (S > 1) {
            for(s in 2:S) {
                eHapsCurrent_tc[, reference_panel_SNPs, s] <- eHapsCurrent_tc[, reference_panel_SNPs, 1]
            }
        }

    } else {

        ## chose some haps at random to fill in eHaps
        ## cols_to_replace <- sample(1:ncol(reference_haps), K)
        ## choose some haps using sampling from PCA approach
        for(s in 1:S) {
            ## get new cols each time? inefficient! revisit?
            cols_to_replace <- sample_haps_to_use(
                rhb = rhb,
                ref_alleleCount_at_L = ref_alleleCount[rh_in_L, ],
                N_haps = N_haps,
                nRefSNPs = nRefSNPs,
                K = K,
                max_snps = 1000,
                max_haps_to_build = 2000,
                max_haps_to_project = 20000
            )
            ## need to add real noise, otherwise, can cause problem
            eHapsCurrent_tc[, reference_panel_SNPs, s] <-
                0.95 * inflate_fhb_t(rhb_t = rhb_t, haps_to_get = cols_to_replace - 1, nSNPs = nRefSNPs) +
                0.05 * array(runif(prod(K * nRefSNPs)), c(K, nRefSNPs))
            ## do not do e and 1-e here
            if (sum(!reference_panel_SNPs) > 0) {
                eHapsCurrent_tc[, !reference_panel_SNPs, s] <- 0.5
            }
        }

        if (reference_iterations > 0) {

            print_message(paste0("Running reference EM with ", N_haps, " reference haplotypes"))
            out <- run_reference_EM(
                eHapsCurrent_tc = eHapsCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                hapSumCurrent_tc = hapSumCurrent_tc,
                sigmaCurrent_m = sigmaCurrent_m,
                priorCurrent_m = priorCurrent_m,
                reference_iterations = reference_iterations,
                rhb = rhb,
                rh_in_L = rh_in_L,
                N_haps = N_haps, nCores = nCores,
                tempdir = tempdir,
                regionName = regionName, nSNPs = nSNPs, nGrids = nGrids, L = L, nGen = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, minRate = minRate, maxRate = maxRate, pseudoHaploidModel = pseudoHaploidModel, reference_phred = reference_phred, grid_distances = grid_distances, reference_shuffleHaplotypeIterations = reference_shuffleHaplotypeIterations, L_grid = L_grid, grid = grid, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, shuffle_bin_radius = shuffle_bin_radius, snps_in_grid_1_based = snps_in_grid_1_based, outputdir = outputdir, plotHapSumDuringIterations = plotHapSumDuringIterations)
            eHapsCurrent_tc <- out$eHapsCurrent_tc
            alphaMatCurrent_tc <- out$alphaMatCurrent_tc
            hapSumCurrent_tc <- out$hapSumCurrent_tc
            sigmaCurrent_m <- out$sigmaCurrent_m
            priorCurrent_m <- out$priorCurrent_m
        }

        ## now - recall - reference_panel_SNPs defines SNPs that exist
        ## if there are some that do not, add in noise
        ## need this noise! important
        n1 <- sum(!reference_panel_SNPs)
        if (sum(n1) > 0) {
            ## add in noise at other positions
            for(s in 1:S) {
                eHapsCurrent_tc[, !reference_panel_SNPs, s] <- matrix(
                    runif(n1 * K, min = 0.4, max = 0.6),
                    nrow = K,
                    ncol = n1
                )
            }
        }
        
    }
    
    print_message("Done initializing paramters using reference haplotypes")


    return(
        list(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            hapSumCurrent_tc = hapSumCurrent_tc,
            sigmaCurrent_m = sigmaCurrent_m,
            priorCurrent_m = priorCurrent_m,
            reference_panel_SNPs = reference_panel_SNPs,
            ref_alleleCount = ref_alleleCount
        )
    )


}


## central function controlling the loading of reference haplotypes

#' @export
get_haplotypes_from_reference <- function(
    reference_haplotype_file,
    reference_legend_file,
    reference_sample_file,
    reference_populations,
    pos,
    tempdir = tempdir(),
    regionName = "test",
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    chr,
    niterations = 40,
    extraction_method = "all", ## both for
    load_rhb_method = "R"
) {

    print_message("Begin get haplotypes from reference")
    print_message("Load and validate reference legend header")
    
    legend_header <- as.character(unlist(read.table(reference_legend_file, nrow = 1, sep = " ")))
    validate_legend_header(legend_header)

    print_message("Load and validate reference legend")
    colClasses <- get_reference_colClasses(
        reference_sample_file,
        reference_populations,
        chr
    )
    out <- load_reference_legend(
        legend_header,
        reference_legend_file,
        regionStart,
        regionEnd,
        buffer
    )
    legend <- out$legend
    legend_position <- out$legend_position
    position_in_haps_file <- out$position_in_haps_file

    validate_reference_legend(legend, reference_legend_file)

    out <- extract_validate_and_load_haplotypes(
        legend = legend,
        pos = pos,
        reference_haplotype_file = reference_haplotype_file,
        position_in_haps_file = position_in_haps_file,
        regionName = regionName,
        tempdir = tempdir,
        colClasses = colClasses,
        niterations = niterations,
        extraction_method = extraction_method,
        load_rhb_method = load_rhb_method
    )

    print_message(paste0("Succesfully extracted ", ncol(out[["rhb3"]]), " haplotypes from reference data"))
    print_message("End get haplotypes from reference")

    return(out)

}

validate_reference_legend <- function(
    legend,
    reference_legend_file
) {
    if (nrow(legend) == 0) {
        stop(paste0("After extraction, the legend file had no rows. Please check the intersection between the regionStart, regionEnd, buffer and the reference_legend_file:", reference_legend_file))
    }
    legend_snps <- paste(
        legend[, "position"], legend[,"a0"], legend[, "a1"],
        sep = "-"
    )
    position <- as.numeric(as.character(legend[, "position"]))
    x <- position <= 0
    if (sum(x) > 0) {
        stop(paste0("There are variants with position <= 0 in the reference legend file ", reference_legend_file, ". One such example is ", legend[which.max(x), "position"], " which occurs at entry ", which.max(x)))
    }
    ## re-use previous functionality
    chr <- "dummy"
    decoy_legend <- data.frame(chr, position, legend[,"a0"], legend[, "a1"])
    validate_pos(decoy_legend, chr = chr, stop_file_name = paste0("The reference legend file ", reference_legend_file, " "))
    ## t <- table(legend_snps)
    ## if (sum(t > 1) > 0) {
    ##     ## argh R - get character not factor
    ##     m <- match(names(t[t>1])[1], legend_snps)
    ##     example <- sapply(
    ##         legend[m, c("position", "a0", "a1")],
    ##         as.character
    ##     )
    ##     stop(
    ##         paste0(
    ##             "There are ", sum(t > 1), " duplicate row ids. ",
    ##             "One such example is ",
    ##             paste0(example, collapse = " ")
    ##         )
    ##     )
    ## }
    return(NULL)
}


validate_haps <- function(
    haps,
    lines_to_get
) {
    if (nrow(haps) != length(lines_to_get)) 
        stop ("Something has gone wrong during processing the reference haplotype and legend files. Please double check the legend file and haplotype file have the same number of entries, other than headers, and there are SNPs in the target region in the reference haplotype file")
    if (sum(is.na(haps)) > 0)
        stop ("The reference_haplotype_file contains entries other than 0 or 1")
    if (sum(haps != 0 & haps != 1) > 0)
        stop ("The reference_haplotype_file contains entries other than 0 or 1")
    if (ncol(haps) <= 1)
        stop ("The reference_haplotype_file contains less than 2 haplotypes. Please use more at least 2 reference haplotypes. If you think this message is an error, please check the reference haplotype file is space separated")
    return(NULL)
}


validate_legend_header <- function(
    legend_header
) {
    if (length(legend_header) < 3)
        stop("The reference_legend_file has fewer than 3 columns. Perhaps it is not a space separated file?")
    if (is.na(match("position", legend_header)))
        stop ("Cannot find 'position' column in reference_legend_file")
    if (is.na(match("a0", legend_header)))
        stop ("Cannot find reference allele 'a0' column in reference_legend_file")
    if (is.na(match("a1", legend_header)))
        stop ("Cannot find alternate allele 'a1' column in reference_legend_file")
    return(NULL)
}



compare_reference_haps_against_alleleCount <- function(
    alleleCount,
    ref_alleleCount,
    outputdir,
    regionName
) {

    all_cor <- suppressWarnings(
        cor(alleleCount[, 3], ref_alleleCount[, 3], use = "pairwise.complete.obs")
    )
    low_maf_cor <- NA
    high_maf_cor <- NA

    w <- alleleCount[, 3] > 0.05 & alleleCount[, 3] < 0.95
    if (sum(w) > 1) {
        high_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], ref_alleleCount[w, 3], use = "pairwise.complete.obs")
        )
    }

    w <- (alleleCount[,3] < 0.05 | alleleCount[,3] > 0.95)
    if (sum(w) > 1) {
        low_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], ref_alleleCount[w, 3], use = "pairwise.complete.obs")
        )
    }

    print_message(paste0("The following correlations are observed between allele frequencies estimated from the sequencing pileup and the reference haplotype counts:"))
    print_message(paste0(round(all_cor, 3), " for all SNPs"))
    print_message(paste0(round(high_maf_cor, 3), " for > 5% MAF SNPs"))
    print_message(paste0(round(low_maf_cor, 3), " for < 5% MAF SNPs"))

    out_plot <- file.path(outputdir, "plots", paste0("alleleFrequency_pileup_vs_reference_haplotypes.", regionName, ".png"))
    print_message(paste0("A plot of allele frequencies from sequencing pileup vs reference haplotype counts is at:", out_plot))
    png(out_plot, height = 500, width = 500)
    plot(alleleCount[, 3], ref_alleleCount[, 3], xlab = "Allele frequency from pileup", ylab = "Allele frequency from reference haplotypes")
    dev.off()

    return(NULL)

}




run_reference_EM <- function(
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    hapSumCurrent_tc,
    sigmaCurrent_m,
    priorCurrent_m,
    reference_iterations,
    rhb,
    rh_in_L,
    N_haps,
    nCores,
    tempdir,
    regionName,
    nSNPs,
    nGrids,
    L,
    nGen,
    emissionThreshold,
    alphaMatThreshold,
    expRate,
    minRate,
    maxRate,
    pseudoHaploidModel,
    reference_phred,
    grid_distances,
    reference_shuffleHaplotypeIterations,
    L_grid,
    grid,
    plot_shuffle_haplotype_attempts,
    shuffle_bin_radius,
    snps_in_grid_1_based,
    outputdir,
    plotHapSumDuringIterations
) {

    ## note - for haplotype shuffling
    ## after   (iteration - 1), find the spots to examine
    ## during  (iteration - 0), examine the shuffling
    ## after   (iteration - 0), do the shuffling
    print_message("Begin EM using reference haplotypes")
    S <- dim(eHapsCurrent_tc)[3]
    K <- dim(eHapsCurrent_tc)[1]
    iteration <- 1
    
    for(iteration in 1:reference_iterations) {
        
        print_message(paste0("Reference iteration = ", iteration))

        ## this loads the break information if the proper iteration (does check)
        out <- get_nbreaks(iteration = iteration, tempdir = tempdir, regionName = regionName, shuffleHaplotypeIterations = reference_shuffleHaplotypeIterations, nGrids = nGrids, shuffle_bin_nSNPs = shuffle_bin_nSNPs, shuffle_bin_radius = shuffle_bin_radius, S = S)
        nbreaks <- out$nbreaks
        list_of_break_results <- out$list_of_break_results

        ## here
        out <- single_reference_iteration(
            eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, sigmaCurrent_m = sigmaCurrent_m, priorCurrent_m = priorCurrent_m, N_haps = N_haps, nCores = nCores, tempdir = tempdir, regionName = regionName, L = L, grid = grid, grid_distances = grid_distances, nGen = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, maxRate = maxRate, minRate = minRate, reference_phred = reference_phred, nbreaks = nbreaks, list_of_break_results = list_of_break_results, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, iteration = iteration,
            rhb = rhb,
            rh_in_L = rh_in_L
        )

        eHapsCurrent_tc <- out$gammaSum_tc
        alphaMatCurrent_tc <- out$alphaMatSum_tc
        hapSumCurrent_tc <- out$hapSum_tc
        sigmaCurrent_m <- out$sigmaSum_m
        sigmaSum_m_unnormalized <- out$sigmaSum_m_unnormalized
        priorCurrent_m <- out$priorSum_m
        list_of_fromMat <- out$list_of_fromMat

        if (iteration %in% (reference_shuffleHaplotypeIterations - 1)) {
            ## in before iteration, define the breaks
            define_and_save_breaks_to_consider(tempdir = tempdir, regionName = regionName, sigmaSum_m_unnormalized = sigmaSum_m_unnormalized, L_grid = L_grid, grid_distances = grid_distances, nGrids = nGrids, nGen = nGen, minRate = minRate, maxRate = maxRate, iteration = iteration, shuffle_bin_radius = shuffle_bin_radius, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts)
        }
        if (sum(nbreaks) > 0) {
            out <- apply_better_switches_if_appropriate(list_of_fromMat = list_of_fromMat, nbreaks = nbreaks, list_of_break_results = list_of_break_results, eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, grid = grid, iteration = iteration, snps_in_grid_1_based = snps_in_grid_1_based, tempdir = tempdir, regionName = regionName, grid_distances = grid_distances, L_grid = L_grid, outputdir = outputdir, is_reference = TRUE, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts)
            eHapsCurrent_tc <- out$eHapsCurrent_tc
            alphaMatCurrent_tc <- out$alphaMatCurrent_tc
        }

        if (plotHapSumDuringIterations) {
            ## these are the new values
            interim_plotter(
                outputdir = outputdir,
                regionName = regionName,
                iteration = iteration,
                L_grid = L_grid,
                hapSumCurrent_tc = hapSumCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                sigmaCurrent_m = sigmaCurrent_m,
                N = N_haps,
                is_reference = TRUE
            )
        }

    }

    print_message("Done EM using reference haplotypes")

    return(
        list(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            hapSumCurrent_tc = hapSumCurrent_tc,
            sigmaCurrent_m = sigmaCurrent_m,
            priorCurrent_m = priorCurrent_m
        )
    )

}





single_reference_iteration <- function(
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    sigmaCurrent_m,
    priorCurrent_m,
    N_haps,
    nCores,
    tempdir,
    regionName,
    L,
    grid,
    grid_distances,
    nGen,
    rhb,
    rh_in_L,
    emissionThreshold = 1e-4,
    alphaMatThreshold = 1e-4,
    expRate = 0.5,
    maxRate = 100,
    minRate = 0.1,
    reference_phred = 20,
    nbreaks = NA,
    list_of_break_results = NULL,
    plot_shuffle_haplotype_attempts = FALSE,
    iteration = 1,
    maxDifferenceBetweenReads = 1000,
    maxEmissionMatrixDifference = 1e10
) {

    sampleRanges <- getSampleRange(N_haps, nCores)
    transMatRate_tc_H <- get_transMatRate_m(method = "pseudoHaploid", sigmaCurrent_m)

    K <- nrow(eHapsCurrent_tc)
    nSNPs <- ncol(eHapsCurrent_tc)
    S <- dim(eHapsCurrent_tc)[[3]]
    nGrids <- ncol(alphaMatCurrent_tc) + 1

    if (is.na(nbreaks[1])) {
        nbreaks <- rep(0, S)
    }

    ## unclear how much bounding might make a difference in the future
    ## leave this in for now, possibly revisit
    if (tail(grid, 1) == (length(grid) - 1)) {
        rescale_eMatGrid_t <- FALSE
        bound_eMatGrid_t <- FALSE
    } else {
        rescale_eMatGrid_t <- TRUE
        bound_eMatGrid_t <- TRUE
    }
    
    out2 <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = new_subset_of_single_reference_iteration,
        tempdir = tempdir,
        regionName = regionName,
        eHapsCurrent_tc = eHapsCurrent_tc,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        priorCurrent_m = priorCurrent_m,
        grid = grid,
        list_of_break_results = list_of_break_results,
        nbreaks = nbreaks,        
        reference_phred = reference_phred,
        sampleRanges = sampleRanges,
        plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts,
        iteration = iteration,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        maxEmissionMatrixDifference = maxEmissionMatrixDifference,
        rescale_eMatGrid_t = rescale_eMatGrid_t,
        bound_eMatGrid_t = bound_eMatGrid_t,
        rhb = rhb,
        rh_in_L = rh_in_L
    )

    check_mclapply_OK(out2, "There has been an error generating the input. Please see error message above")

    out <- calculate_updates(
        out2 = out2, sampleRanges = sampleRanges, N = N_haps,
        nGen = nGen, expRate = expRate, minRate = minRate, maxRate = maxRate,
        emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, L = L, grid_distances = grid_distances
    )

    sigmaSum_m <- out$sigmaSum_m
    sigmaSum_m_unnormalized <- out$sigmaSum_m_unnormalized
    priorSum_m <- out$priorSum_m
    alphaMatSum_tc <- out$alphaMatSum_tc
    gammaSum_tc <- out$gammaSum_tc
    hapSum_tc <- out$hapSum_tc

    if (sum(nbreaks) > 0) {
        list_of_fromMat <- lapply(1:S, function(s) {
            nbreak <- nbreaks[s]
            fromMat <- array(0, c(nbreak, K, K))            
            for(i in 1:length(sampleRanges)) {
                fromMat <- fromMat + out2[[i]][["list_of_fromMat"]][[s]]
            }
            return(fromMat)
        })
    } else {
        list_of_fromMat <- NULL
    }
    
    return(
        list(
            sigmaSum_m = sigmaSum_m,
            sigmaSum_m_unnormalized = sigmaSum_m_unnormalized,
            priorSum_m= priorSum_m,
            alphaMatSum_tc = alphaMatSum_tc,
            gammaSum_tc = gammaSum_tc,
            hapSum_tc = hapSum_tc,
            list_of_fromMat = list_of_fromMat
        )
    )

}

new_subset_of_single_reference_iteration <- function(
    sampleRange,
    tempdir,
    regionName,
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    transMatRate_tc_H,
    priorCurrent_m,
    grid,
    list_of_break_results,
    nbreaks,
    reference_phred,
    sampleRanges,
    plot_shuffle_haplotype_attempts,
    iteration,
    rhb,
    rh_in_L,
    rescale_eMatGrid_t = FALSE,
    bound_eMatGrid_t = FALSE,
    maxDifferenceBetweenReads = 1000,
    maxEmissionMatrixDifference = 1e10,
    suppressOutput = 1
) {

    K <- nrow(eHapsCurrent_tc)
    nSNPs <- ncol(eHapsCurrent_tc)
    S <- dim(eHapsCurrent_tc)[[3]]
    nGrids <- ncol(alphaMatCurrent_tc) + 1

    ## no, make like before
    ## maxDifferenceBetweenReads <- 1 / (10 **(-(reference_phred / 10)))
    ## maxEmissionMatrixDifference <- 1 / (10 **(-(reference_phred / 10)))
    
    alphaHat_t <- array(0, c(K, nGrids))
    betaHat_t <- array(0, c(K, nGrids))
    gamma_t <- array(0, c(K, nGrids))
    eMatGrid_t <- array(0, c(K, nGrids))

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

    alphaMatCurrentX_tc <- array(0, dim(alphaMatCurrent_tc))
    rcpp_make_alphaMatSumX_tc(alphaMatCurrent_tc = alphaMatCurrent_tc, alphaMatCurrentX_tc = alphaMatCurrentX_tc, transMatRate_tc_H = transMatRate_tc_H)
    
    ## for updating
    gammaSum0_tc <- array(0, c(K, nSNPs, S))
    gammaSum1_tc <- array(0, c(K, nSNPs, S))
    alphaMatSum_tc <- array(0, c(K, nGrids - 1, S))
    ## alphaMatSum_tc <- array(0, c(K, nGrids, S)) ## make 1 bigger for now!
    hapSum_tc <- array(0, c(K, nGrids, S))
    priorSum_m <- array(0, c(K, S))
    
    if (sum(nbreaks) > 0) {
        list_of_fromMat <- lapply(1:S, function(s) {
            return(array(0, c(nbreaks[s], K, K)))
        })
        return_gammaK <- TRUE
    } else {
        list_of_fromMat <- NULL
        return_gammaK <- FALSE
    }

    ## this is dinky - just make anyway
    pshaM <- min(20, sampleRange[2]) ## how many to plot                    
    list_of_fbd_store <- lapply(1:S, function(s) {
        return(as.list(1:pshaM))
    })
    i_core <- match(sampleRange[1], sapply(sampleRanges, function(x) x[[1]]))
    if ((sum(nbreaks) > 0) & plot_shuffle_haplotype_attempts & i_core == 1) {
        save_fbd_store <- TRUE
    } else {
        save_fbd_store <- FALSE
    }

    for(iSample1 in sampleRange[1]:sampleRange[2]) {

        ## note - makes almost no difference to push the entire below into c++
        ## if (iSample1 == 10) {
        ##     print("remove me")
        ##     suppressOutput <- 0
        ## } else {
        ##     suppressOutput <- 1
        ## }

        ## storage is more efficient at least!
        ## do raw later!
        ## or something more efficient...
        reference_hap <- inflate_fhb(rhb, haps_to_get = iSample1 - 1, nSNPs = nSNPs)
        fbsoL <- reference_fbh(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrentX_tc = alphaMatCurrentX_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            priorCurrent_m = priorCurrent_m,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            suppressOutput = suppressOutput,
            reference_hap = reference_hap,
            rh_in_L = rh_in_L,
            gammaSum0_tc = gammaSum0_tc,
            gammaSum1_tc = gammaSum1_tc,
            alphaMatSum_tc = alphaMatSum_tc,
            hapSum_tc = hapSum_tc,
            priorSum_m = priorSum_m,
            list_of_break_results = list_of_break_results,
            list_of_fromMat = list_of_fromMat,
            nbreaks = nbreaks,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t,
            gamma_t = gamma_t,
            eMatGrid_t = eMatGrid_t,
            ehh_h1_S = ehh_h1_S,
            ehh_h1_D = ehh_h1_D,
            ehh_h0_S = ehh_h0_S,
            ehh_h0_D = ehh_h0_D,
            grid = grid,
            return_gammaK = FALSE,
            return_extra = FALSE,
            reference_phred = reference_phred,
            list_of_fbd_store = list_of_fbd_store,
            save_fbd_store = save_fbd_store,
            iSample1 = iSample1,
            pshaM = pshaM
        )

    }

    if (save_fbd_store) {
        save(list_of_fbd_store, file = file_fbdStore(tempdir, regionName, iteration))
    }
    

    ## remove one entry
    ## alphaMatSum_tc <- alphaMatSum_tc[, -1, ]
    ## in c++, for normal mode, avoided some multiplications
    ## here, add them back in
    rcpp_finalize_alphaMatSum_tc(
        alphaMatSum_tc = alphaMatSum_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaMatCurrentX_tc = alphaMatCurrentX_tc
    )

    return(
        list(
            gammaSum0_tc = gammaSum0_tc,
            gammaSum1_tc = gammaSum1_tc,
            alphaMatSum_tc = alphaMatSum_tc,
            hapSum_tc = hapSum_tc,
            priorSum_m = priorSum_m,
            list_of_fromMat = list_of_fromMat
        )
    )

}










get_reference_colClasses <- function(
    reference_sample_file,
    reference_populations,
    chr
) {
    colClasses <- NA
    if (reference_sample_file != "" & !is.na(reference_populations[1])) {
        reference_samples <- read.table(reference_sample_file, header = TRUE)
        validate_reference_sample_file(reference_samples, reference_sample_file)
        keep_samples <- which(
            is.na(
                match(
                    rep(as.character(reference_samples[, "POP"]), each = 2),
                    reference_populations
                )
            ) == FALSE
        )
        print_message(paste0("Keeping ", length(keep_samples) / 2, " out of ", nrow(reference_samples), " reference samples from populations:", paste0(reference_populations, collapse = ",")))
        colClasses <- rep("NULL", nrow(reference_samples) * 2)
        colClasses[keep_samples] <- "integer"
    }
    return(colClasses = colClasses)
}






load_reference_legend <- function(
    legend_header,
    reference_legend_file,
    regionStart,
    regionEnd,
    buffer
) {
    ## extract reference
    r <- match("position", legend_header)
    if (is.na(regionStart) == FALSE) {
        extract_condition <- paste0(
            "$", r, " >= ", regionStart - buffer, " && ",
            "$", r, "<=", regionEnd + buffer
        )
    } else {
        extract_condition <- "NR>1" ## do not load the header
    }
    legend_and_range <- system(
        paste0(
            "gunzip -c ", shQuote(reference_legend_file), " | ", "awk '{if(", extract_condition,  ") {print NR\" \"$0}}'"
        ), intern = TRUE
    )
    a <- strsplit(legend_and_range, " ")
    if (length(a) == 0) {
        stop(paste0("No SNPs were loaded from the reference_legend_file:", reference_legend_file, " for region with regionStart=", regionStart, ", regionEnd=", regionEnd, " and buffer=", buffer, ". Perhaps the reference legend file does not span the region requested?"))
    }
    col <- as.integer(sapply(a, function(x) x[1]))
    legend <- t(sapply(a, function(x) x[-1]))
    colnames(legend) <- legend_header
    legend_position <- as.integer(legend[, "position"])
    position_in_haps_file <- as.integer(sapply(a, function(x) x[1])) - 1
    return(
        list(
            legend = legend,
            legend_position = legend_position,
            position_in_haps_file = position_in_haps_file
        )
    )
}



## extract the relevant rows from the reference haplotypes file
## do this using awk to avoid loading way too many sites
extract_validate_and_load_haplotypes <- function(
    legend,
    pos,
    reference_haplotype_file,
    position_in_haps_file,
    regionName,
    tempdir,
    colClasses,
    niterations,
    extraction_method = "all",
    load_rhb_method = "R"
) {

    print_message("Extract reference haplotypes")
    temp_haps_file <- file.path(tempdir, paste0("haps.", regionName, ".txt.gz"))

    pos_snps <- paste(pos[,2], pos[,3], pos[,4], sep = "-")
    legend_snps <- paste(
        legend[, "position"], legend[,"a0"], legend[, "a1"],
        sep = "-"
    )

    both_snps <- intersect(legend_snps, pos_snps)

    validate_pos_and_legend_snps_for_niterations_equals_1(legend_snps, pos_snps, niterations)
    print_and_validate_reference_snp_stats(pos_snps, legend_snps, both_snps)

    ## only extract those we need
    ## note - legend_snps validated - so no duplicate legend SNPs
    hap_snps_to_extract <- is.na(match(legend_snps, both_snps)) == FALSE
    hap_snps_position_in_pos <- is.na(match(pos_snps, both_snps)) == FALSE
    lines_to_get <- position_in_haps_file[hap_snps_to_extract]
    rh_in_L <- which(hap_snps_position_in_pos)    

    nSNPs <- nrow(pos)
    n_haps_snps <- max(lines_to_get) ## stop at this, anyway
    reference_haps <- NA
    rhb1 <- NA
    rhb2 <- NA
    rhb3 <- NA
    ref_alleleCount1 <- NA
    ref_alleleCount2 <- NA
    ref_alleleCount3 <- NA    
    if ((extraction_method == "all") | (extraction_method == "hap_v1")) {
        ##
        reference_haps <- load_haps_at_positions_no_NAs(
            reference_haplotype_file = reference_haplotype_file,
            lines_to_get = lines_to_get,
            colClasses = colClasses,
            nSNPs = nSNPs,
            tempdir = tempdir
        )
        rhb1 <- make_rhb_from_rhi(reference_haps)
        ref_alleleCount <- array(NA, c(nSNPs, 3))
        ref_alleleCount[rh_in_L, 1] <- rowSums(reference_haps)
        ref_alleleCount[rh_in_L, 2] <- ncol(reference_haps)
        ref_alleleCount[, 3] <- ref_alleleCount[, 1] / ref_alleleCount[, 2]
        ref_alleleCount1 <- ref_alleleCount
    }
    if ((extraction_method == "all") | (extraction_method == "hap_v2")) {
        ##
        out <- load_rhb_at_positions_no_NAs(
            reference_haplotype_file = reference_haplotype_file,
            lines_to_get = lines_to_get,
            colClasses = colClasses,
            nSNPs = nSNPs,
            tempdir = tempdir,
            n_haps_snps = n_haps_snps,
            rh_in_L = rh_in_L,            
            load_rhb_method = "R"
        )
        rhb2 <- out[["rhb"]]
        ref_alleleCount2 <- out[["ref_alleleCount"]]
    }
    if ((extraction_method == "all") | (extraction_method == "hap_v3")) {
        ##
        out <- load_rhb_at_positions_no_NAs(
            reference_haplotype_file = reference_haplotype_file,
            lines_to_get = lines_to_get,
            colClasses = colClasses,
            nSNPs = nSNPs,
            tempdir = tempdir,
            n_haps_snps = n_haps_snps,
            rh_in_L = rh_in_L,
            load_rhb_method = "Rcpp"
        )
        rhb3 <- out[["rhb"]]
        ref_alleleCount3 <- out[["ref_alleleCount"]]
    }

    return(
        list(
            reference_haps = reference_haps,
            rh_in_L = rh_in_L,
            rhb1 = rhb1,
            rhb2 = rhb2,
            rhb3 = rhb3,
            ref_alleleCount1 = ref_alleleCount1,
            ref_alleleCount2 = ref_alleleCount2,
            ref_alleleCount3 = ref_alleleCount3            
        )
    )
}


remove_NA_columns_from_haps <- function(haps) {
    na_sum <- colSums(is.na(haps))
    x <- na_sum != 0 & na_sum != nrow(haps)
    if (sum(x) > 0)
        stop(
            "The missing haplotype entry '-' for a reference male sample on chromosome X does not make up the entirety of a sample haplotype for at least one sample for the region of interest"
        )
    if (sum(na_sum == 0) == 0)
        stop("There are no viable haplotypes from the reference samples for the region of interest")
    if (sum(na_sum == nrow(na_sum)) > 0)
        print_message(paste0("Removing ", X, " out of ", Y, " male haplotypes from the reference haplotypes"))
    if (sum(na_sum != 0) > 0) {
        haps <- haps[, na_sum == 0, drop = FALSE]
    }
    return(haps)
}


## lines to get is 1-based (right) what to get from haps file
## hap_snps_position_in_pos is 1-based, where this goes wrt big file
load_haps_at_positions_no_NAs <- function(
    reference_haplotype_file,
    lines_to_get,
    colClasses,
    nSNPs,
    tempdir = tempdir()
) {

    temp_extract_file <- tempfile(fileext = ".txt", tmpdir = tempdir)
    cat(
        lines_to_get,
        file = temp_extract_file, sep = "\n"
    )
    awk_command <- paste0(
        "awk '{if (NR==FNR) {a[$1]=$1; next} ",
        "if (FNR in a) {print $0}}' "
    )
    command <- paste0(
        "gunzip -c ", shQuote(reference_haplotype_file), " | ",
        awk_command, shQuote(temp_extract_file), " - "
    )
    ## 
    haps <- data.table::fread(
        cmd = command,
        colClasses = colClasses,
        sep = " ",
        na.strings = "-",
        data.table = FALSE
    )
    colnames(haps) <- NULL
    unlink(temp_extract_file)
    ##
    ## argh - data.table doing weird things
    ## 
    haps2 <- matrix(0L, nrow = nrow(haps), ncol = ncol(haps))
    haps2[] <- as.matrix(haps)
    haps <- haps2
    ## haps <- base::as.matrix(haps) ## as.matrix
    haps <- remove_NA_columns_from_haps(haps)
    validate_haps(haps, lines_to_get)
    ## convert to rhb around here!
    return(haps)
    ##
    ## no longer "expand" haps here. now make the size of whatever was loaded
    ## control later with rh_in_L, saying how it compares to pos
    ## 
    ## haps_both <- array(NA, c(nSNPs, ncol(haps)))
    ## haps_both[hap_snps_position_in_pos, ] <- haps
    ## 
}




## was originally
## cols_to_replace <- sample(1:ncol(reference_haps), K)
## but this fails on low K, trivial examples
##
## note - on toy data - this works very well
##        ## how often is this working properly
##        sapply(1:1000, function(i) {
##           set.seed(i)
##            cols_to_replace <- sample_haps_to_use(reference_haps, K)
##            a <- sort(colSums(reference_haps[, cols_to_replace]))
##            sum(c(a[1] == 1, a[2] == 3, a[3] == 6))
##        })
##
sample_haps_to_use <- function(
    rhb,
    ref_alleleCount_at_L,
    N_haps,
    nRefSNPs,
    K,
    max_snps = 1000,
    max_haps_to_build = 2000,
    max_haps_to_project = 20000
) {
    
    print_message("Begin determining haplotypes to use for initialization")
    if (max_haps_to_project < K) {
        max_haps_to_project <- K
    }

    if (nRefSNPs > max_snps) {
        
        ## make weight proportional to allele frequency, i.e. sample more frequent
        a <- ref_alleleCount_at_L[, 3]
        ## this can error, but I am uncertain how, as it is fairly well tested
        ## make it throw a better error if it does
        if (sum(is.na(a)) > 0) {
            print_message("BEGIN ERROR REPORTING")
            print_message(paste0("There are ", sum(is.na(a)), " NA reference frequency values out of ", length(a), " SNPs"))
            print_message(paste0("Here is the first one, occuring in row ", which.max(is.na(a))))
            print(ref_alleleCount_at_L[which.max(is.na(a)), ])
            stop("There was an error with sample_haps_to_use. At least one allele frequency for a reference SNP is NA when it should be a 0-1 double. Please report this and the message above")
        }
        a[a > 0.5] <- 1 - a[a > 0.5]
        prob <- a / sum(a)
        prob[is.na(prob)] <- 0
        keep_snps <- sort(sample(1:nRefSNPs, size = min(sum(prob > 0), max_snps), replace = FALSE, prob = prob))
    } else {
        keep_snps <- NA
        prob <- rep(1 / nRefSNPs, nRefSNPs)
    }
    ##
    if (N_haps > max_haps_to_build) {
        ## meh, do at random
        keep_samples <- sort(sample(1:N_haps, size = max_haps_to_build, replace = FALSE))
    } else {
        keep_samples <- 1:N_haps
    }
    
    ## inflate here!
    h <- inflate_fhb(rhb, haps_to_get = keep_samples - 1, nSNPs = nRefSNPs)
    if (!is.na(keep_snps[1])) {
        h <- h[keep_snps, ] ## meh
    }

    ##
    ## first, build the pca, using the building samples
    ## 
    c <- rowMeans(h) ## rows are SNPs
    h <- h - c
    h2 <- t(h) %*% (h)
    out <- eigen(h2, symmetric = TRUE)

    ## keep this many columns
    eigen_cols_to_keep <- min(ncol(out$vectors), max(K, which.max((cumsum (out$values / sum(out$values))) > 0.50) - 1))
    
    ## only do if appropriate number of haps
    if ((max_haps_to_build < N_haps) & (max_haps_to_project < N_haps)) {
        ##
        ## now, project "everyone" onto this
        ## don't do everyone as R kmeans not linear in N (although close!)
        ##
        vs <- out$vectors[, 1:eigen_cols_to_keep]
        ls <- out$values[1:eigen_cols_to_keep]
        ds <- diag(ls)
        ds_inv <- diag(1 / ls)
        ## next line - if many SNPs have been provided with AF 0, this could break. this fixes it
        keep_snps <- sort(sample(1:nRefSNPs, size = min(sum(prob > 0), max_snps), replace = FALSE, prob = prob))
        print_message("Perform projection")
        h_all <- inflate_fhb(rhb, haps_to_get = keep_samples - 1, nSNPs = nRefSNPs)
        if (!is.na(keep_snps[1])) {
            h_all <- h_all[keep_snps, ] ## meh
        }
        h2b <- t(h) %*% (h_all - c) 
        b <- (t(h2b) %*% vs %*% ds_inv) ## this is the potentially slow one!
        print_message("Done performing projection")
    } else {
        ## can do directly
        b <- out$vectors[, 1:eigen_cols_to_keep, drop = FALSE]
    }
    ## change values slightly
    for(icol in 1:eigen_cols_to_keep) {
        b[, icol] <- b[, icol] * out$values[icol]
    }
    ##
    b <- round(b, 3)
    b <- b[, colSums(b != 0) > 0, drop = FALSE]
    ## yuck - if K < number of unique, will fail
    local_K <- min(K, nrow(unique(b)))
    if (local_K == nrow(b)) {
        local_K <- local_K - 1
    }
    print_message("Perform K-means")
    kmeans_results <- suppressWarnings(kmeans(b, centers = local_K, iter.max = 100, nstart = 10))
    print_message("Done performing K-means")    
    ##
    cluster_membership <- kmeans_results$cluster
    ## note - keep_samples is fine here - importantly, it is how the projection is done
    cols_to_replace <- choose_cols_to_replace(cluster_membership, keep_samples, local_K, K)
    ##

    print_message("Done determining haplotypes to use for initialization")
    
    return(cols_to_replace)
    
}

choose_cols_to_replace <- function(cluster_membership, keep_samples, local_K, K) {
    cols_to_replace <- array(NA, K)
    for(k in 1:local_K) {
        ## god dammit R, "sample" changes definition on single valued input
        w <- which(cluster_membership == k)
        if (length(w) == 1) {
            cols_to_replace[k] <- w
        } else {
            cols_to_replace[k] <- sample(keep_samples[w], size = 1)
        }
    }
    if (local_K < K) {
        ## remaining, choose at random
        w <- setdiff(keep_samples, cols_to_replace[1:local_K])
        cols_to_replace[-c(1:local_K)] <- sort(sample(w, K - local_K, replace = FALSE))
    }
    return(sort(cols_to_replace))
}



make_sampleReads_from_hap <- function(rh_in_L, reference_phred, reference_hap) {
    sampleReads <- lapply(
        rh_in_L,
        function(x) {
        return(list(
            0, x - 1,
            reference_phred * (2 * reference_hap[x] - 1),
            x - 1
        ))
    }
    )
    return(sampleReads)
}
