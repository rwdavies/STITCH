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
    windowSNPs,
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
    snps_in_grid_1_based
) {

    print_message("Begin initializing paramters using reference haplotypes")
    ## get reference haplotypes matched to posfile
    ## NA's where there are no match
    ## only for populations of interest
    reference_haps <- get_haplotypes_from_reference(
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
        niterations = niterations
    )
    reference_panel_SNPs <- is.na(reference_haps[, 1]) == FALSE

    ## build ref_alleleCount
    ref_alleleCount <- array(NA, c(nrow(reference_haps), 3))
    ref_alleleCount[, 1] <- rowSums(reference_haps)
    ref_alleleCount[, 2] <- ncol(reference_haps)
    ref_alleleCount[, 3] <- ref_alleleCount[, 1] / ref_alleleCount[, 2]

    if (is.null(alleleCount) == FALSE) {
        compare_reference_haps_against_alleleCount(
            alleleCount = alleleCount,
            reference_haps = reference_haps,
            outputdir = outputdir,
            regionName = regionName
        )
    }

    if (K > ncol(reference_haps)) {

        ## fill in rest with noise
        print_message("You have set K to be more than the number of reference haplotypes. The rest of the K ancestral haplotypes will be filled with noise to start")
        h <- t(reference_haps[reference_panel_SNPs, ])
        for(s in 1:S) {
            eHapsCurrent_tc[1:ncol(reference_haps), reference_panel_SNPs, s] <- h
        }

    } else if (K == ncol(reference_haps)) {

        print_message("There are exactly as many reference haplotypes as K. Using these haplotypes directly as the initial estimate of the ancestral haplotypes")
        ## shrink them from 0 -> e and 1 -> (1-e)
        e <- 0.001
        reference_haps[reference_haps == 0] <- e
        reference_haps[reference_haps == 1] <- (1 - e)
        reference_haps_t <- t(reference_haps[reference_panel_SNPs, ])
        for(s in 1:S) {
            eHapsCurrent_tc[, reference_panel_SNPs, s] <- reference_haps_t
        }

    } else {

        N_haps <- ncol(reference_haps)
        reference_bundling_info <- get_bundling_position_information(
            N = N_haps,
            nCores = nCores,
            blockSize = inputBundleBlockSize
        )

        make_fake_sample_reads_from_haplotypes(
            reference_haps = reference_haps,
            nCores = nCores,
            reference_bundling_info = reference_bundling_info,
            tempdir = tempdir,
            N_haps = N_haps,
            reference_phred = reference_phred,
            regionName = regionName
        )

        ## snap to grid
        if (nSNPs > nGrids) {
            snap_reads_to_grid(
                N = N_haps,
                nCores = nCores,
                regionName = regionName,
                tempdir = tempdir,
                bundling_info = reference_bundling_info,
                grid = grid,
                whatKindOfReads = "referenceSampleReads",
                downsampleToCov = downsampleToCov,
                allSampleReads = NULL
            )
        }

        ## chose some haps at random to fill in eHaps
        ## cols_to_replace <- sample(1:ncol(reference_haps), K)
        ## choose some haps using sampling from PCA approach
        for(s in 1:S) {
            ## get new cols each time!
            cols_to_replace <- sample_haps_to_use(reference_haps, K)
            reference_haps_t <- t(reference_haps[, cols_to_replace])
            ## add some noise
            e <- 0.01
            reference_haps_t[reference_haps_t == 0] <- e
            reference_haps_t[reference_haps_t == 1] <- (1 - e)
            ##
            eHapsCurrent_tc[, reference_panel_SNPs, s] <- reference_haps_t[, reference_panel_SNPs]
            ## make the other positions have more noise
            n1 <- sum(!reference_panel_SNPs)
            ## these do not matter for now
            eHapsCurrent_tc[, !reference_panel_SNPs, s] <- 0.5
        }

        if (reference_iterations > 0) {
            out <- run_EM_on_reference_sample_reads(
                eHapsCurrent_tc = eHapsCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                hapSumCurrent_tc = hapSumCurrent_tc,
                sigmaCurrent_m = sigmaCurrent_m,
                priorCurrent_m = priorCurrent_m,
                reference_iterations = reference_iterations,
                N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, nSNPs = nSNPs, nGrids = nGrids, L = L, nGen = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, minRate = minRate, maxRate = maxRate, pseudoHaploidModel = pseudoHaploidModel, reference_phred = reference_phred, grid_distances = grid_distances, reference_shuffleHaplotypeIterations = reference_shuffleHaplotypeIterations, L_grid = L_grid, grid = grid, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, shuffle_bin_radius = shuffle_bin_radius, snps_in_grid_1_based = snps_in_grid_1_based, outputdir = outputdir)
            eHapsCurrent_tc <- out$eHapsCurrent_tc
            alphaMatCurrent_tc <- out$alphaMatCurrent_tc
            hapSumCurrent_tc <- out$hapSumCurrent_tc
            sigmaCurrent_m <- out$sigmaCurrent_m
            priorCurrent_m <- out$priorCurrent_m
        }

        ## now - recall - reference_panel_SNPs defines SNPs that exist
        ## if there are some that do not, add in noise
        ## choose some haps using sampling from PCA approach
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
    niterations = 40
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

    haps <- extract_validate_and_load_haplotypes(
        legend,
        pos,
        reference_haplotype_file,
        position_in_haps_file,
        regionName,
        tempdir,
        colClasses,
        niterations
    )

    print_message(paste0("Succesfully extracted ", ncol(haps), " haplotypes from reference data"))
    print_message("End get haplotypes from reference")

    return(haps)

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
    x <- as.numeric(as.character(legend[, "position"])) <= 0
    if (sum(x) > 0) {
        stop(paste0("There are variants with position <= 0 in the reference legend file ", reference_legend_file, ". One such example is ", legend[which.max(x), "position"], " which occurs at entry ", which.max(x)))
    }
    t <- table(legend_snps)
    if (sum(t > 1) > 0) {
        ## argh R - get character not factor
        m <- match(names(t[t>1])[1], legend_snps)
        example <- sapply(
            legend[m, c("position", "a0", "a1")],
            as.character
        )
        stop(
            paste0(
                "There are ", sum(t > 1), " duplicate row ids. ",
                "One such example is ",
                paste0(example, collapse = " ")
            )
        )
    }
    return(NULL)
}


validate_haps <- function(
    haps,
    both_snps
) {
    if (nrow(haps) != length(both_snps))
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
    reference_haps,
    outputdir,
    regionName
) {

    all_cor <- suppressWarnings(
        cor(alleleCount[, 3], rowSums(reference_haps), use = "pairwise.complete.obs")
    )
    low_maf_cor <- NA
    high_maf_cor <- NA

    w <- alleleCount[,3] > 0.05 & alleleCount[,3] < 0.95
    if (sum(w) > 1)
        high_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], rowSums(reference_haps[w, ]), use = "pairwise.complete.obs")
        )

    w <- alleleCount[,3] < 0.05 | alleleCount[,3] > 0.95
    if (sum(w) > 1)
        low_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], rowSums(reference_haps[w, ]), use = "pairwise.complete.obs")
        )

    print_message(paste0("The following correlations are observed between allele frequencies estimated from the sequencing pileup and the reference haplotype counts:"))
    print_message(paste0(round(all_cor, 3), " for all SNPs"))
    print_message(paste0(round(high_maf_cor, 3), " for > 5% MAF SNPs"))
    print_message(paste0(round(low_maf_cor, 3), " for < 5% MAF SNPs"))

    N_haps <- ncol(reference_haps)

    out_plot <- file.path(outputdir, "plots", paste0("alleleFrequency_pileup_vs_reference_haplotypes.", regionName, ".png"))
    print_message(paste0("A plot of allele frequencies from sequencing pileup vs reference haplotype counts is at:", out_plot))
    png(out_plot, height = 500, width = 500)
    plot(alleleCount[, 3], rowSums(reference_haps) / N_haps, xlab = "Allele frequency from pileup", ylab = "Allele frequency from reference haplotypes")
    dev.off()

    return(NULL)

}


make_fake_sample_reads_from_haplotypes <- function(
    reference_haps,
    nCores,
    reference_bundling_info,
    tempdir,
    N_haps,
    reference_phred,
    regionName
) {

    print_message(paste0("Begin convert reference haplotypes to internal input format"))

    sampleRanges <- getSampleRange(N_haps, nCores)
    non_NA_cols <- which(is.na(reference_haps[ , 1]) == FALSE)

    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        function(sampleRange) {
        for(iSample in sampleRange[1]:sampleRange[2]) {
            sampleReads <- rcpp_make_sampleReads_from_hap(
                non_NA_cols,
                reference_phred,
                reference_hap = reference_haps[, iSample]
            )
            save(
                sampleReads,
                file = file_referenceSampleReads(tempdir, iSample, regionName)
            )
            if (length(reference_bundling_info) > 0) {
                last_in_bundle <- reference_bundling_info$matrix[iSample, "last_in_bundle"]
                if (last_in_bundle == 1) {
                    bundle_inputs_after_generation(
                        bundling_info = reference_bundling_info,
                        iBam = iSample,
                        dir = tempdir,
                        regionName = regionName,
                        what = "referenceSampleReads"
                    )
                }
            }
        }
        return(NULL)
    })

    check_mclapply_OK(out, "There has been an error generating the input. Please see error message above")

    print_message("End convert reference haplotypes to internal input format")

}


run_EM_on_reference_sample_reads <- function(
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    hapSumCurrent_tc,
    sigmaCurrent_m,
    priorCurrent_m,
    reference_iterations,
    N_haps,
    nCores,
    reference_bundling_info,
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
    outputdir
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
        out <- single_reference_iteration(eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, sigmaCurrent_m = sigmaCurrent_m, priorCurrent_m = priorCurrent_m, N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, L = L, grid = grid, grid_distances = grid_distances, nGen = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, maxRate = maxRate, minRate = minRate, reference_phred = reference_phred, nbreaks = nbreaks, break_results = break_results, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, iteration = iteration)

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
            apply_better_switches_if_appropriate(list_of_fromMat = list_of_fromMat, nbreaks = nbreaks, list_of_break_results = list_of_break_results, eHapsCurrent_tc = eHapsCurrent_tc, alphaMatCurrent_tc = alphaMatCurrent_tc, grid = grid, iteration = iteration, snps_in_grid_1_based = snps_in_grid_1_based, tempdir = tempdir, regionName = regionName, grid_distances = grid_distances, L_grid = L_grid, outputdir = outputdir)
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
                N = N,
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
    reference_bundling_info,
    tempdir,
    regionName,
    L,
    grid,
    grid_distances,
    nGen,
    emissionThreshold = 1e-4,
    alphaMatThreshold = 1e-4,
    expRate = 0.5,
    maxRate = 100,
    minRate = 0.1,
    reference_phred = 20,
    nbreaks = NA,
    break_results = NULL,
    plot_shuffle_haplotype_attempts = FALSE,
    iteration = 1
) {

    sampleRanges <- getSampleRange(N_haps, nCores)
    transMatRate_tc_H <- get_transMatRate_m(method = "pseudoHaploid", sigmaCurrent_m)

    maxDifferenceBetweenReads <- 1 / (10 **(-(reference_phred / 10)))
    maxEmissionMatrixDifference <- 1 / (10 **(-(reference_phred / 10)))
    K <- nrow(eHapsCurrent_tc)
    nSNPs <- ncol(eHapsCurrent_tc)
    S <- dim(eHapsCurrent_tc)[[3]]
    nGrids <- ncol(alphaMatCurrent_tc) + 1

    if (is.na(nbreaks[1])) {
        nbreaks <- rep(0, S)
    }

    out2 <- mclapply(sampleRanges, mc.cores = nCores, function(sampleRange) {

        alphaHat_t <- array(0, c(K, nGrids))
        betaHat_t <- array(0, c(K, nGrids))
        gamma_t <- array(0, c(K, nGrids))
        eMatGrid_t <- array(0, c(K, nGrids))

        ## for updating
        gammaSum0_tc <- array(0, c(K, nSNPs, S))
        gammaSum1_tc <- array(0, c(K, nSNPs, S))
        alphaMatSum_tc <- array(0, c(K, nGrids - 1, S))
        hapSum_tc <- array(0, c(K, nGrids, S))
        priorSum_m <- array(0, c(K, S))

        if (sum(nbreaks) > 0) {
            list_of_fromMat <- lapply(nbreaks, function(nbreak) {
                array(0, c(nbreak, K, K))
            })
        } else {
            list_of_fromMat <- NULL
        }

        bundledSampleReads <- NULL
        for(iSample in sampleRange[1]:sampleRange[2]) {

            out <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = reference_bundling_info,
                bundledSampleReads = bundledSampleReads,
                what = "referenceSampleReads"
            )
            sampleReads <- out$sampleReads
            bundledSampleReads <- out$bundledSampleReads
            
            fbsoL <- forwardBackwardHaploid(
                sampleReads = sampleReads,
                eHapsCurrent_tc = eHapsCurrent_tc,
                alphaMatCurrent_tc = alphaMatCurrent_tc,
                transMatRate_tc_H = transMatRate_tc_H,
                priorCurrent_m = priorCurrent_m,
                pRgivenH1 = as.double(0),
                pRgivenH2 = as.double(0),
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                maxEmissionMatrixDifference = maxEmissionMatrixDifference,
                suppressOutput = 1,
                model = -1, ## irrelavent here
                run_pseudo_haploid = as.integer(0), ## just run haploid
                blocks_for_output = array(0, c(1, 1)),
                update_in_place = TRUE,
                gammaSum0_tc = gammaSum0_tc,
                gammaSum1_tc = gammaSum1_tc,
                alphaMatSum_tc = alphaMatSum_tc,
                hapSum_tc = hapSum_tc,
                priorSum_m = priorSum_m,
                pass_in_alphaBeta = TRUE,
                alphaHat_t = alphaHat_t,
                betaHat_t = betaHat_t,
                gamma_t = gamma_t,
                eMatGrid_t = eMatGrid_t,
                grid = grid,
                Jmax = 10,
                prev_list_of_alphaBetaBlocks = list()
            )

            for(s in 1:S) {
                if (nbreaks[s] > 0) {
                    for(iBreak in 1:nbreaks[s]) {
                        from <- break_results[iBreak, "left_grid_break_0_based"]
                        to <- break_results[iBreak, "right_grid_break_0_based"]
                        hp1 <- fbsoL$gamma_t[, from + 1]
                        hp2 <- fbsoL$gamma_t[, to + 1]
                        list_of_fromMat[[s]][iBreak, ] <-
                            list_of_fromMat[[s]][iBreak, ] + hp1 %*% t(hp2)
                    }
                    if (plot_shuffle_haplotype_attempts) {
                        i_core <- match(sampleRange[1], sapply(sampleRanges, function(x) x[[1]]))
                        if (i_core == 1) {
                            M <- min(20, sampleRange[2]) ## how many to plot
                            if (iSample == 1) {
                                fbd_store <- as.list(1:M)
                            }
                            if (iSample <= M) {
                                ## this needs to be R <- fbd_store[[iSample]]$gammaK_t
                                fbd_store[[iSample]] <- list(gammaK_t = fbsoL[["gamma_t"]])
                            }
                            if (iSample == M) {
                                save(fbd_store, file = file_fbdStore(tempdir, regionName, iteration))
                                fbd_store <- NULL
                            }
                        }
                    }
                }
            }
            
        }
        
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
    })

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
                fromMat_c <- fromMat_c + out2[[i]][["list_of_fromMat"]][[s]]
            }
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


make_sampleReads_from_hap <- function(non_NA_cols, reference_phred, reference_hap) {
    sampleReads <- lapply(
        non_NA_cols,
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





get_reference_colClasses <- function(
    reference_sample_file,
    reference_populations,
    chr
) {
    colClasses <- NA
    if (reference_sample_file != "") {
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
    niterations
) {

    print_message("Extract reference haplotypes")
    temp_haps_file <- file.path(tempdir, paste0("haps.", regionName, ".txt.gz"))
    temp_extract_file <- file.path(tempdir, paste0("extract.", regionName, ".txt"))

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

    ## extract rows from haplotypes file
    cat(
        position_in_haps_file[hap_snps_to_extract],
        file = temp_extract_file, sep = "\n"
    )

    awk_command <- paste0(
        "awk '{if (NR==FNR) {a[$1]=$1; next} ",
        "if (FNR in a) {print $0}}' "
    )
    command <- paste0(
        "gunzip -c ", shQuote(reference_haplotype_file), " | ",
        awk_command, shQuote(temp_extract_file), " - | ",
        "gzip > ", shQuote(temp_haps_file)
    )
    out <- system(command)

    print_message("Load reference haplotypes")
    haps <- read.table(
        temp_haps_file,
        colClasses = colClasses,
        sep = " ",
        na.strings = "-"
    )

    unlink(temp_haps_file)
    haps <- as.matrix(haps)

    haps <- remove_NA_columns_from_haps(haps)
    validate_haps(haps, both_snps)

    haps_both <- array(NA, c(nrow(pos), ncol(haps)))
    haps_both[hap_snps_position_in_pos, ] <- haps

    return(haps_both)
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
    haps <- haps[, na_sum == 0, drop = FALSE]
    return(haps)
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
sample_haps_to_use <- function(reference_haps, K, max_snps = 1000, max_samples = 2000) {
    reference_haps <- reference_haps[!is.na(reference_haps[, 1]), ]
    if (nrow(reference_haps) > max_snps) {
        ## make weight proportional to allele frequency
        a <- rowSums(haps) / ncol(haps)
        a[a > 0.5] <- 1 - a[a > 0.5]
        prob <- a / sum(a)
        keep <- sort(sample(1:nrow(haps), size = max_snps, replace = FALSE, prob = prob))
        reference_haps <- reference_haps[keep, ]
    }
    ##
    if (ncol(reference_haps) > max_samples) {
        ## meh, do at random
        keep_samples <- sort(sample(1:ncol(reference_haps), size = max_samples, replace = FALSE))
        reference_haps <- reference_haps[, keep_samples]
    } else {
        keep_samples <- 1:ncol(reference_haps)
    }
    ##
    h <- reference_haps
    c <- rowMeans(h)
    h <- h - c
    h <- t(h)
    h2 <- h %*% t(h)
    out <- eigen(h2)
    v <- out$vectors
    ## at least K, at most >50% fit
    k <- min(ncol(v), max(K, which.max((cumsum (out$values / sum(out$values))) > 0.50) - 1))
    ##
    b <- v[, 1:k, drop = FALSE]
    for(i in 1:k) {
        b[, i] <- b[, i] * out$values[i]
    }
    ##
    b <- round(b, 3)
    b <- b[, colSums(b != 0) > 0, drop = FALSE]
    ## yuck - if K < number of unique, will fail
    local_K <- min(K, nrow(unique(b)))
    out2 <- suppressWarnings(kmeans(b, centers = local_K, iter.max = 100, nstart = 10))
    ## remove useless columns
    ## k-nearest neighbours
    ## take average of members
    ## eHapsCurrent <- array(NA, c(nrow(reference_haps), K))
    cols_to_replace <- array(NA, K)
    for(k in 1:K) {
        if (k <= local_K) {
            w <- which(out2$cluster == k)
        } else {
            ## any remaining sample
            w <- setdiff(keep_samples, cols_to_replace[1:(k - 1)])
        }
        ## yuck - just sample
        cols_to_replace[k] <- sample(keep_samples[w], 1)
        ## eHapsCurrent[, k] <- rowSums(reference_haps[, keep_samples[w]]) / length(w)
        ## eHapsCurrent[, k] <- reference_haps[, sample(keep_samples[w], 1)]
    }
    return(cols_to_replace)
    ##
    ## WHAT THE FUCK KMEANS
    ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ## plot(v[, 1] + rnorm(n = 60) / 100, v[, 2] + rnorm(n = 60) / 100, col = cbPalette[out2$cluster])
    ## ## but proportional to values (multiply by that value)
    ##
    ## pdf("~/temp.pdf", height = 8, width = 20)
    ## par(mfrow = c(2, 5))
    ## col <-
    ## for(i in 1:10) {
    ##     plot(v[, 1], v[, i], col = col)
    ## }
    ## dev.off()
}


