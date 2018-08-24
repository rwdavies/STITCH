## in a simple refill, look every 100 SNP. in each one with an average below a threshold, sample from other haplotypes to re-copy
## across all haplotypes, when this happen, inject noise
refillSimple <- function(
    hapSum_t,
    nGrids,
    K,
    gammaSum_t,
    N,
    L_grid,
    grid,
    distance_between_check = 5000 ## check every 5000 bp
) {
    nSNPs <- ncol(gammaSum_t)
    a <- floor(L_grid / distance_between_check)
    a <- c(match(unique(a), a))
    if (length(a) == 1) {
        region_starts <- 1
        region_ends <- nGrids
    } else {
        region_starts <- a
        region_ends <- c(a[-1] - 1, length(L_grid))
    }
    avHapSumInBlock <- array(0, c(length(region_starts), K))
    for(i in 1:length(region_starts)) {
        avHapSumInBlock[i, ] <- rowSums(hapSum_t[, region_starts[i]:region_ends[i], drop = FALSE])         
    }
    hapFreqMin <- N / K ## less than average - aggressive!
    replaceBlock <- avHapSumInBlock < (hapFreqMin)
    ## 
    ## for each k, for each region, fill in
    ## within each continuous region, fill in with respect to frequencies of all other haplotypes
    ##
    k_to_replace <- which(colSums(replaceBlock) > 0)
    ever_changed <- array(FALSE, ncol(gammaSum_t))
    for (k in k_to_replace) {
        ## change into intervals
        z1 <- replaceBlock[, k]
        x <- (1:(length(z1)-1))[diff(z1)!=0]
        start <- c(1,x+1)
        end <- c(x,length(z1))
        ## only take "bad" regions
        start <- start[z1[start]==TRUE]
        end <- end[z1[end]==TRUE    ]
        ## OK!
        if (length(start) > 0) {
            for(iR in 1:length(start)) {
                ## now - in each region, get sums of haplotypes
                p <- colSums(matrix(avHapSumInBlock[start[iR]:end[iR],],ncol=K))
                p[k] <- 0 ## do not re-sample!
                p <- p / sum(p)
                replacement <- sample(K, 1, prob = p)
                r1 <- region_starts[start[iR]]
                r2 <- region_ends[end[iR]]
                snps_to_replace <- ((r1 - 1) <= grid) & (grid <= (r2 - 1))
                ## now need to get grid as well
                gammaSum_t[k, snps_to_replace] <- gammaSum_t[replacement, snps_to_replace]
                ever_changed[snps_to_replace] <- TRUE
            }
        }
    }
    print_message(
        paste0(
            "Refill infrequently used haplotypes - ",
            round(100 * sum(replaceBlock == TRUE )/ prod(dim(replaceBlock)), 1),
            "% of regions replaced"
    ))
    return(
        list(
            gammaSum_t = gammaSum_t,
            ever_changed = ever_changed
        )
    )
}



## breaks is X
## nbreaks is number of breaks
## rStart is X
## rEnd is X
get_nbreaks <- function(
    iteration,
    tempdir,
    regionName,
    shuffleHaplotypeIterations,
    nGrids,
    shuffle_bin_nSNPs = NULL,
    shuffle_bin_radius = 5000
) {
    break_results <- NULL ## matrix otherwise
    if(is.na(match(iteration, shuffleHaplotypeIterations))==FALSE) {
        if (is.null(shuffle_bin_radius) == FALSE) {
            ## have previously (previous iteration) ended with define_breaks_to_consider
            ## now, load break_results
            ## NOTE - can be NULL, but unlikely
            load(file_break_results(tempdir, regionName))
        } else {
            ## mimic format here as well
            start <- round(shuffle_bin_nSNPs / 2) + 1
            if (nGrids >= (start + shuffle_bin_nSNPs)) {
                ##
                left_break <- seq(start, nGrids - shuffle_bin_nSNPs, by = shuffle_bin_nSNPs)
                nbreaks <- length(left_break) - 1
                break_results <- cbind(
                    left_grid_break_0_based = left_break - 1,
                    left_grid_focal_0_based = left_break + (shuffle_bin_nSNPs - start),
                    right_grid_focal_0_based = left_break + (shuffle_bin_nSNPs - start + 1),
                    right_grid_break_0_based = left_break + shuffle_bin_nSNPs - 1
                )
            }
        }
    }
    if (is.null(break_results)) {
        nbreaks <- 0
    } else {
        nbreaks <- nrow(break_results)
    }
    return(
        list(
            break_results = break_results,
            nbreaks = nbreaks
        )
    )
}


determine_switch_order <- function(fromMat, nbreaks, K) {
    switchOrder <- t(sapply(1:nbreaks,fromMat=fromMat,K=K,function(ib,fromMat,K) {
        fromMatL <- fromMat[ib, , ]
        ## columns are to
        ## rows are from
        ## ie entry fromMat[10,2,3] is from breaks[10] hap 2to breaks[10+1] hap 3
        orderL=array(0, K)
        fromL=1:K
        toL=1:K
        fromR <- array(0, K)
        toR <- array(0, K)
        for(j in 1:(K - 1)) {
            ## 
            x <- which(fromMatL == max(fromMatL), arr.ind = TRUE)[1, ]
            f <- x[1]
            t <- x[2]
            fromR[j] <- fromL[f]
            toR[j] <- toL[t]
            ## now - remove from consideration
            fromMatL <- fromMatL[-f,-t]
            fromL <- fromL[-f]
            toL <- toL[-t]
        }
        fromR[K] <- fromL
        toR[K] <- toL
        ## re-order by fromR
        return(toR[order(fromR)])
    }))
    return(switchOrder)
}


getBetterSwitchesSimple <- function(
    fromMat,
    nbreaks,
    break_results,
    K,
    eHapsFuture_t,
    alphaMatFuture_t,
    grid,
    snps_in_grid_1_based,
    iteration = 1
) {
    ## greedy version
    ## choose, in order, one which means fewest moves
    switchOrder <- determine_switch_order(fromMat, nbreaks, K) 
    whichIsBest <- as.integer(apply(switchOrder,1,function(x) sum(x==(1:K)))==K)
    ##
    ## determine the order of subsequent regions
    ##
    ## now - do the unravelling
    tempMat <- array(0, c(nbreaks + 1,K))
    tempMat[1, ] <- 1:K
    currentState <- 1:K
    nextState <- 1:K
    for (iBreak in 1:(nbreaks)) {
        ## how this works - first, start in your current state
        ## next, choose new state from switchOrder. copy that state, then remember where to go next
        currentState <- nextState
        for(k in 1:K) {
            ## first - the current value is wherever you are
            tempMat[iBreak + 1, k] <- switchOrder[iBreak, currentState[k]]
            nextState[k] <- tempMat[iBreak + 1, k]
        }
    }
    ##
    ## do the shuffling
    ##
    print_message(paste0(
        "Shuffle haplotypes - Iteration ", iteration, " - change ", sum(whichIsBest!=1), " intervals out of ", nbreaks, " considered"
    ))
    ## 0-based start and end of the grids with the switches
    grid_starts <- c(0, break_results[, "left_grid_focal_0_based"] + 1)
    grid_ends <- c(break_results[, "left_grid_focal_0_based"], ncol(alphaMatFuture_t) - 1)
    ## do switches around focal points
    for(iBreak in 1:(nbreaks + 1)) {
        ## grids, 0-based, can just be a number
        which_grids <- grid_starts[iBreak]:grid_ends[iBreak] 
        ## snps, 1-based
        which_snps <- c(
            snps_in_grid_1_based[grid_starts[iBreak] + 1, "snps_start"]:
            snps_in_grid_1_based[grid_ends[iBreak] + 1, "snps_end"]
        )
        ## 
        permL <- tempMat[iBreak, ]
        eHapsFuture_t[, which_snps] <- eHapsFuture_t[permL, which_snps]
        alphaMatFuture_t[, which_grids] <- alphaMatFuture_t[permL, which_grids]
    }
    ## add in noise around the breaks
    for(iBreak in 1:nbreaks) {
        if (whichIsBest[iBreak] == 0) {
            s <- break_results[iBreak, "left_grid_break_0_based"]
            e <- (break_results[iBreak, "right_grid_break_0_based"] - 1) ## do not actually include this one
            which_grids <- s:e
            ## snps, 1-based
            which_snps <- c(
                snps_in_grid_1_based[s + 1, "snps_start"]:
                snps_in_grid_1_based[e + 1, "snps_end"]
            )
            ##
            norm_component <- abs(1 - seq(0, 2, length = length(which_snps)))
            ## fill in with more noise closer to the break
            for(ii in 1:length(which_snps)) {
                eHapsFuture_t[, which_snps[ii]] <-
                    norm_component[ii] * eHapsFuture_t[, which_snps[ii]] +
                    (1 - norm_component[ii]) * runif(K)
            }
            ## just reset these, can be re-determined pretty quickly
            alphaMatFuture_t[, which_grids + 1] <- matrix(1 / K, nrow = K, ncol = length(which_grids))
        }
    }
    return(
        list(
            eHapsFuture_t = eHapsFuture_t,
            alphaMatFuture_t = alphaMatFuture_t,
            whichIsBest = whichIsBest
        )
    )
}


## this function attempts to use unnormalized recombination rate
## to find regions with elevated rate beyond bound,
## and then to return a series of peaks, start and end points
## so that in the subsequent iteration,
## we can examine how making breaks on those changes fitting
##
##
## what is returned:
## break from_snp = what SNP to start to look at
##       central SNP = first in pair of SNPs to central
##       to_snp = end of block
## so if shuffle_bin_radius is 2000, and looking between SNPs at positions a and b
## then smooth rate over those SNPs is between round((a + b) / 2) minus 2000 and plus 2000
define_and_save_breaks_to_consider <- function(
    tempdir,
    regionName,
    sigmaSum_unnormalized,
    L_grid,
    grid_distances,
    nGrids,
    nGen,
    minRate,
    maxRate,
    iteration = NULL,
    shuffle_bin_radius = 2000,
    plot_shuffle_haplotype_attempts = FALSE
) {
    ## if too few SNPs, do not bother
    ## EM algorithm should be fine, won't get stuck this way
    if ((nGrids < 5)) {
        break_results <- NULL
        save(break_results, file = file_break_results(tempdir, regionName))
    } else {
        ## generate smoothed rate
        smoothed_rate <- rcpp_make_smoothed_rate(
            sigmaSum = sigmaSum_unnormalized,
            sigma_rate = -log(sigmaSum_unnormalized) / grid_distances,
            L_grid = L_grid,
            grid_distances = grid_distances, ## when no grid, this is diff(L)
            nGrids = nGrids, ## when no grid, this equals nSNPs 
            shuffle_bin_radius = shuffle_bin_radius
        )
        ## normalize against max rate
        smoothed_rate <- smoothed_rate / 
            ((nGen * maxRate * 1/100/1000000) * shuffle_bin_radius * 2)
        ##
        out <- choose_points_to_break(
            smoothed_rate = smoothed_rate,
            nGrids = nGrids,
            L_grid = L_grid,
            shuffle_bin_radius = shuffle_bin_radius,
            grid_distances = grid_distances
        )
        break_results <- out$results
        break_thresh <- out$thresh
        ## NULL return unlikely but probably probable
        if (plot_shuffle_haplotype_attempts == FALSE) {
            save(break_results, file = file_break_results(tempdir, regionName))
        } else {
            ## add in useful things for plot on subsequent iteration
            ##binSize <- 10000
            ##rate <- array(0, max(L)) ## this could be huge?
            ##for(i in 1:(nSNPs - 1)) {
            ##    rate[L[i]:L[i + 1]] <- (-log(sigmaSum_unnormalized[i]) / grid_distances[i])
            ## }
            x1 <- exp(-nGen * minRate * grid_distances/100/1000000) # lower
            x2 <- exp(-nGen * maxRate * grid_distances/100/1000000) # upper
            ## make matrix with useful things
            realized_rate <- (-log(sigmaSum_unnormalized) / grid_distances)
            realized_rate_no_NA <- realized_rate
            realized_rate_no_NA[is.na(realized_rate)] <- 0
            cumu_rate <- cumsum(realized_rate_no_NA) - realized_rate_no_NA[1]
            recomb_usage <- cbind(
                grid_distances = grid_distances,
                min_prob = x1,
                realized_prob = sigmaSum_unnormalized,
                max_prob = x2,
                min_rate = (-log(x1) / grid_distances),
                realized_rate = realized_rate,
                max_rate = (-log(x2) / grid_distances),
                cumu_rate = cumu_rate
            )
            save(break_thresh, smoothed_rate, recomb_usage, break_results, file = file_break_results(tempdir, regionName))
            ## save extra copy - useful for debugging etc
            save(break_thresh, smoothed_rate, recomb_usage, break_results, file = file_break_results(tempdir, regionName, iteration))
        }
    }
    return(NULL)
}


##
## returns one entry of average rate between pairs of SNPs
## value is average for middle point +/- shuffle_bin_radius
## 
make_smoothed_rate <- function(
    sigmaSum_unnormalized,
    sigma_rate,
    L_grid,
    grid_distances,
    nGrids,
    shuffle_bin_radius = 5000
) {
    smoothed_rate <- array(0, nGrids - 1) ## between SNPs
    ##
    min_L_grid <- min(L_grid)
    max_L_grid <- max(L_grid)
    for(iSNP in 1:(nGrids - 1)) {
        focal_point <- floor((L_grid[iSNP] + L_grid[iSNP + 1]) / 2)        
        if (
            min_L_grid <= (focal_point - shuffle_bin_radius) &
            (focal_point + shuffle_bin_radius) <= max_L_grid
        ) {
            ## left
            ## so re-call
            ## start with focal SNP, defined by position average of iSNP and iSNP + 1
            ## first region is iSNP_left
            iSNP_left <- iSNP
            bp_remaining <- shuffle_bin_radius
            bp_prev <- focal_point
            while(0 < bp_remaining) {
                bp_to_add <- (bp_prev - L_grid[iSNP_left] + 1)
                if ((bp_remaining - bp_to_add) < 0) {
                    bp_to_add <- bp_remaining
                    bp_remaining <- 0
                } else {
                    bp_remaining <- bp_remaining - bp_to_add
                    bp_prev <- L_grid[iSNP_left]
                }
                ## add bit
                smoothed_rate[iSNP] <- smoothed_rate[iSNP] + bp_to_add * sigma_rate[iSNP_left]
                iSNP_left <- iSNP_left - 1                
            }
            ## right
            ## so recall. start with iSNP, then iSNP_right starts one to the right
            ## at first, difference between focal point (mid-way between iSNP and iSNP + 1)
            ## and iSNP + 1. then rate for that region is for iSNP (iSNP_right -1)
            iSNP_right <- iSNP + 1
            bp_remaining <- shuffle_bin_radius
            bp_prev <- focal_point
            while(0 < bp_remaining) {
                bp_to_add <- (L_grid[iSNP_right] - bp_prev + 1)
                if ((bp_remaining - bp_to_add) < 0) {
                    bp_to_add <- bp_remaining
                    bp_remaining <- 0
                } else {
                    bp_remaining <- bp_remaining - bp_to_add
                    bp_prev <- L_grid[iSNP_right] ## reset
                }
                ## add bit
                smoothed_rate[iSNP] <- smoothed_rate[iSNP] +
                    bp_to_add * sigma_rate[iSNP_right - 1]
                iSNP_right <- iSNP_right + 1
            }
        } else {
            smoothed_rate[iSNP] <- NA
        }
    }
    return(smoothed_rate)
}



## 
choose_points_to_break <- function(
    smoothed_rate,
    nGrids,
    L_grid,
    shuffle_bin_radius,
    grid_distances,
    min_breaks = 10,
    max_breaks = NULL
) {
    if (is.null(max_breaks)) {
        ## not ideal
        max_breaks <- max(1000, round(nGrids / 100))
    }
    smoothed_rate[1] <- NA
    smoothed_rate[length(smoothed_rate)] <- NA    
    ## now, those > 1 are great candidates
    ## alternatively, if nGen not quite right (or rates not changed!)
    ## anything really "peaky"
    ## remember, nothing forces it to make the change later, so can be really liberal here
    ## only downside is later calculation is slow-ish
    ## so choose as threshold lower of 1, or 50th percentile
    thresh <-  min(1, quantile(smoothed_rate[smoothed_rate != 0], probs = 0.50, na.rm = TRUE))
    ideal <- smoothed_rate > thresh
    ideal[is.na(ideal)] <- FALSE
    available <- array(TRUE, length(smoothed_rate))
    available[is.na(smoothed_rate)] <- FALSE
    best <- order(smoothed_rate, decreasing = TRUE)
    ##
    results <- NULL
    ## if nothing really available, return
    if (sum(available) <= 3) {
        ## return results as NULL
        return(list(results = results, thresh = thresh))
    }
    ## keep going to get
    iBest <- 1 ## ordered on "best" ordered vector
    while(iBest < nGrids) {
        ## now find start, end for this
        ## this is between grids "snp_best" and grid "snp_best" + 1
        snp_best <- best[iBest]
        ## only consider if have >1 Grid either side
        if (sum(available[snp_best + -1:1]) == 3) {
            ## no matter what, one further away
            ## return 1-based as well
            snp_left <- determine_where_to_stop(
                smoothed_rate,
                available,
                snp_best,
                thresh,
                nGrids,
                side = "left"
            )
            ##
            snp_right <- determine_where_to_stop(
                smoothed_rate,
                available,
                snp_best,
                thresh,
                nGrids,
                side = "right"
            )
            ## store
            ## everything here 1_based coordinates
            ## left_break_1_based is left grid before elevated rate switching
            ## left_focal_1_based is left grid point around highest local rate
            ## right_focal_1_based is right grid point around highest local rate
            ## right_break_1_based is right grid after elevated rate switching
            colnames <- c(
                "left_grid_break_0_based",
                "left_grid_focal_0_based",
                "right_grid_focal_0_based",
                "right_grid_break_0_based"
            )
            results <- rbind(
                results,
                c(snp_left - 1, snp_best - 1, snp_best, snp_right - 1)
            )
            ## nuke out (do not consider for future) at least
            ## anything within 2 * shuffle_bin_radius of focal SNP or within decreasing region
            to_nuke <- get_snps_to_nuke(
                grid_distances, nGrids, shuffle_bin_radius, snp_best, snp_left, snp_right
            )
            ## for available - remove the outside SNPs. so can be outside for another?
            available[min(to_nuke[1] + 1, snp_best):max(snp_best, to_nuke[2] - 1)] <- FALSE
            ideal[to_nuke[1]:to_nuke[2]] <- FALSE
        } else {
            ## otherwise, either SNP upstream or downstream cannot be considered
            ## remove consideration
            available[snp_best + -1:1] <- FALSE
        }
        if (sum(available, na.rm = TRUE) == 0) {
            iBest <- nGrids
        } else if ((length(results) / 4) == max_breaks) {
            iBest <- nGrids
        } else {
            if ((sum(ideal) == 0) & ((length(results) / 4) >= min_breaks)) {
                    iBest <- nGrids
            } else {
                ## move forward in iBest
                iBest <- iBest + 1
                while((available[best[iBest]] == FALSE) & (iBest < (nGrids - 1))) {
                    iBest <- iBest + 1
                }
            }
        }
    }
    if ((length(results) / 4) > 0) {
        colnames(results) <- colnames
        results <- results[order(results[, "left_grid_focal_0_based"]), , drop = FALSE]        
    }
    return(
        list(
            results = results,
            thresh = thresh
        )
    )
}


get_snps_to_nuke <- function(grid_distances, nGrids, shuffle_bin_radius, snp_best, snp_left, snp_right) {
    ## nuke those within 2 * shuffle_bin_radius of central SNP
    ## or those within snp_left, snp_right
    dist_left <- grid_distances[snp_best] / 2
    ## 
    snp_nuke_left <- snp_best
    while((1 < snp_nuke_left) & (dist_left < (2 * shuffle_bin_radius))) {
        snp_nuke_left <- snp_nuke_left - 1
        dist_left <- dist_left + grid_distances[snp_nuke_left]
    }
    if ((2 * shuffle_bin_radius) <= dist_left) {
        snp_nuke_left <- snp_nuke_left + 1
    }
    ##
    snp_nuke_right <- snp_best
    dist_right <- grid_distances[snp_best] / 2
    while((snp_nuke_right < (nGrids - 1)) & (dist_right < (2 * shuffle_bin_radius))) {
        snp_nuke_right <- snp_nuke_right + 1
        dist_right <- dist_right + grid_distances[snp_nuke_right]
    }
    if ((2 * shuffle_bin_radius) <= dist_right) {
        snp_nuke_right <- snp_nuke_right - 1
    }
    ##
    return(c(
        nuke_left = min(snp_left, snp_nuke_left),
        nuke_right = max(snp_right, snp_nuke_right)
    ))
}


determine_where_to_stop <- function(
    smoothed_rate,
    available,
    snp_best,
    thresh,
    nGrids,
    side = "left"
) {
    if (side == "left") {
        mult <- 1
    } else {
        mult <- -1
    }
    snp_consider <- snp_best
    val_cur <- smoothed_rate[snp_consider]
    val_prev <- smoothed_rate[snp_best]
    snp_min <- snp_consider
    val_min <- smoothed_rate[snp_min]
    c <- 1
    done <- 0
    while(done == 0) {
        snp_consider <- snp_consider + (-1) * mult
        val_cur <- smoothed_rate[snp_consider]
        if (c >= 5) {
            val_prev <- smoothed_rate[snp_consider + 5 * mult]
        }
        c <- c + 1
        if (val_cur < val_min) {
            snp_min <- snp_consider
            val_min <- val_cur
        }
        ## do not continue if would go out of bound
        if ((snp_consider <= 2) | ((nGrids - 3) <= snp_consider)) {
            done <- 1
        }
        ## do not continue if not available i.e. in previous region        
        if (available[snp_consider + (-1) * mult] == FALSE) {
            done <- 1
        }
        ## do not continue if 300% of minimum
        ## also, re-set to minimum value
        if ((3 * val_min) < val_cur) {
            done <- 1
        }
        ## stop if below threshold AND not still decreasing over 5 SNP
        if ((val_cur < thresh) & (val_prev < val_cur)) {
            done <- 1
        }
    }
    return(snp_min)
}



plot_attempt_to_find_shuffles <- function(
    grid_distances,
    L_grid,
    fbd_store,
    tempdir,
    outputdir,
    regionName,
    iteration,
    whichIsBest = NULL
) {
    add_grey_background <- function(L_grid) {
        from <- floor(min(L_grid) / 1e6)
        to <- ceiling(max(L_grid) / 1e6)
        for(f in from:to) {
            abline(v = 1e6 * f, col = "grey")
        }
    }
    ##
    ## the commented out bit at the bottom assumes starting at 0
    ## would be nice but need to re-write to make general
    ##
    ## 0-5Mbp in 2000 bins (1-2000, 2001-4000, etc)
    ## rate binned
    ##rate_binned <- array(0, 5000000 / binSize)
    ##for(iBin in 1:length(rate_binned)) {
    ##    w <- (1 + binSize * (iBin - 1)):(binSize * (iBin - 0))
    ##    rate_binned[iBin] <- sum(rate[w], na.rm = TRUE) / sum(is.na(rate[w]) == FALSE)
    ## }
    ##
    load(file_break_results(tempdir, regionName))
    if (is.null(break_results)) {
        print_message("Insufficient SNPs to try to find shuffled haplotypes")
        return(NULL)
    }
    xlim <- c(head(L_grid)[1], tail(L_grid)[1])    
    ## do plot here
    outname <- file.path(outputdir, "plots", paste0("shuffleHaplotypes.", iteration, ".",regionName,".png"))
    ## make a 5 Mbp segment 60 wide. then bound up and down at 20 and 200
    width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)
    png(outname, height = 30, width = width, res = 100, units = "in")
    par(mfrow = c(3, 1))
    x <- L_grid[-1] - grid_distances    
    xleft <- L_grid[-length(L_grid)]
    xright <- L_grid[-1]
    ## for fbdstore
    midpoints <- L_grid[-1] - grid_distances / 2    
    xleft2 <- c(L_grid[1], midpoints)
    xright2 <- c(midpoints, L_grid[length(L_grid)])
    ##
    ## 0) plot fine and "coarse" scale recombination
    ##
    r <- c(
        range(recomb_usage[, "min_rate"]),        
        range(recomb_usage[, "realized_rate"]),
        range(recomb_usage[, "max_rate"])
    )
    ylim <- log10(range(unlist(r), na.rm = TRUE))
    z <- recomb_usage[, "cumu_rate"]
    z <- z * (1 / max(z)) ## 0-1 scaled
    cumu_rate_scaled <- ylim[2] + z * (ylim[1] - ylim[2])
    plot(x = 0, y = 0, col = "white", ylim = ylim, xlim = xlim, main = "Various recombination rates. Blue is expected min, red is expected max, green is unnormalized desired. Values above max suggest shuffling.\nOrange is cumultative decreasing arbitrarily scaled from 0-1 on same axis")
    add_grey_background(L_grid)
    lines(x = x, y = log10(recomb_usage[, "min_rate"]), ylim = ylim, col = "blue")
    lines(x = x, y = log10(recomb_usage[, "realized_rate"]), ylim = ylim, col = "green")
    ## lines(x = seq(binSize / 2 + 1, 5000000, binSize), y = log10(rate_binned), ylim = ylim, col = "purple")
    lines(x = x, y = log10(recomb_usage[, "max_rate"]), ylim = ylim, col = "red")
    lines(x = x, y = cumu_rate_scaled, ylim = ylim, col = "orange")
    ## then add all the switches?
    ##
    ## 1) find peaks to check
    ##
    ylim <- c(0, max(break_thresh, max(smoothed_rate, na.rm = TRUE)))
    plot(x = 0, y = 0, xlab = "Physical position", ylab = "Rate", main = "Location of shuffles to check", ylim = ylim, xlim = xlim)
    add_grey_background(L_grid)
    lines(x = x, y = smoothed_rate, lwd = 2)
    for(iBreak in 1:nrow(break_results)) {
        abline(v = midpoints[break_results[iBreak, "left_grid_break_0_based"] + 1], col = "red")
        abline(v = midpoints[break_results[iBreak, "left_grid_focal_0_based"] + 1], col = "purple")
        abline(v = midpoints[break_results[iBreak, "right_grid_break_0_based"] + 1], col = "red")
        if (is.null(whichIsBest) == FALSE) {
            col <- c("green", "red")[as.integer(whichIsBest[iBreak]) + 1] ## switch: red = no, green = yes
        } else {
            col <- "orange"
        }
        rect(
            xleft = midpoints[break_results[iBreak, "left_grid_break_0_based"] + 1],
            xright = midpoints[break_results[iBreak, "right_grid_break_0_based"] + 1],
            ybottom = -1, ytop = 0,
            col = col
        ) ##whichIsBest 0 = no switch
    }
    abline(h = break_thresh, col = "black", lwd = 2)
    ##
    ## 2) plot individuals with their switches
    ##
    plot_fbd_store(fbd_store = fbd_store, xleft2 = xleft2, xright2 = xright2, xlim = xlim)
    ##
    for(iBreak in 1:nrow(break_results)) {
        abline(v = midpoints[break_results[iBreak, "left_grid_break_0_based"] + 1], col = "red")
        abline(v = midpoints[break_results[iBreak, "left_grid_focal_0_based"] + 1], col = "purple")
        abline(v = midpoints[break_results[iBreak, "right_grid_break_0_based"] + 1], col = "red")
    }
    add_grey_background(L_grid)
    K <- nrow(fbd_store[[1]]$gammaK_t)
    cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)    
    legend("topright", paste0("anc_hap=", 1:K), col = cbPalette[1:K], lwd = 2)
    dev.off()
}


plot_fbd_store <- function(fbd_store, xleft, xright, xlim, main = "Haplotype usage per-sample", mbp_xlab = FALSE, xleft2 = NULL, xright2 = NULL) {
    if (is.null(xleft2) == FALSE) {
        xleft <- xleft2
        xright <- xright2
    }
    if (mbp_xlab) {
        xlab <- "Physical position (Mbp)"
        xleft <- xleft / 1e6
        xright <- xright / 1e6
        xlim <- xlim / 1e6
        d <- 1e1 / 1e6
    } else {
        xlab <- "Physical position"
        d <- 0
    }
    NN <- length(fbd_store)
    ylim <- c(1, NN + 1)
    plot(x = 0, y = 0, ylim = ylim, xlim = xlim, xlab = xlab, ylab = "Sample", main = main, col = "white")
    cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)
    nr <- length(xleft)
    for(iSample in 1:NN) {
        ## plot each of them, as rectangle?
        R <- fbd_store[[iSample]]$gammaK_t
        ## should sum to 1, plot rectangles of each
        ybottom <- iSample + array(0, nr)
        ytop <- iSample + array(0, nr)
        for(k in 1:nrow(R)) {
            if (ncol(R) != nr) {
                ytop <- ytop + R[k, -ncol(R)]
            } else {
                ytop <- ytop + R[k, ]
            }
            rect(
                xleft = xleft - d,
                xright = xright + d, ## should not be necessary, R artefact I think
                ybottom = ybottom,
                ytop = ytop,
                col = cbPalette[k],
                border = NA
            )
            ybottom <- ytop
        }
    }
}


alpha_col <- function(col, alpha) {
    x <- col2rgb(col) / 255
    return(rgb(x["red", 1], x["green", 1], x["blue", 1], alpha = alpha)    )
}
