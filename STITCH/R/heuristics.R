## in a simple refill, look every 100 SNP.
## in each one with an average below a threshold,
## see who was in that haplotype before. then extend out their choice 
## sample from other haplotypes to re-copy
## across all haplotypes, when this happen, inject noise
refillSimple <- function(
    gammaSum_tc,
    hapSum_tc,
    N,
    L_grid,
    grid,
    distance_between_check = 5000 ## check every 5000 bp
) {
    ## 
    K <- nrow(gammaSum_tc)
    S <- dim(gammaSum_tc)[3]
    nGrids <- ncol(hapSum_tc)
    ##
    a <- floor(L_grid / distance_between_check)
    a <- c(match(unique(a), a))
    if (length(a) == 1) {
        region_starts <- 1
        region_ends <- nGrids
    } else {
        region_starts <- a
        region_ends <- c(a[-1] - 1, length(L_grid))
    }
    avHapSumInBlock <- array(0, c(length(region_starts), K, S))
    ever_changed <- array(FALSE, c(ncol(gammaSum_tc), S))    
    for(s in 1:S) {
        for(i in 1:length(region_starts)) {
            n <- region_ends[i] - region_starts[i] + 1
            avHapSumInBlock[i, , s] <- rowSums(hapSum_tc[, region_starts[i]:region_ends[i], s, drop = FALSE]) / n
        }
        ## really target haplotypes that basically aren't used    
        hapFreqMin <- N * min(c(1 / (10 * K), 1 / 100))
        replaceBlock <- avHapSumInBlock[, , s, drop = FALSE] < (hapFreqMin)
        ## 
        ## for each k, for each region, fill in
        ## within each continuous region, fill in with respect to frequencies of all other haplotypes
        ##
        k_to_replace <- which(colSums(replaceBlock) > 0)
        for (k in k_to_replace) {
            ## change into intervals
            z1 <- replaceBlock[, k, 1, drop = FALSE] ## yes, 1, not s here
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
                    p <- colSums(matrix(avHapSumInBlock[start[iR]:end[iR], , s, drop = FALSE], ncol = K))
                    p[k] <- 0 ## do not re-sample!
                    p <- p / sum(p)
                    replacement <- sample(K, 1, prob = p)
                    r1 <- region_starts[start[iR]]
                    r2 <- region_ends[end[iR]]
                    snps_to_replace <- ((r1 - 1) <= grid) & (grid <= (r2 - 1))
                    ## now need to get grid as well
                    gammaSum_tc[k, snps_to_replace, s] <- gammaSum_tc[replacement, snps_to_replace, s]
                    ever_changed[snps_to_replace, s] <- TRUE
                }
            }
        }
    }
    print_message(
        paste0(
            "Refill infrequently used haplotypes - on average, ",
            round(100 * sum(replaceBlock)/ prod(dim(replaceBlock)), 1),
            "% of regions replaced"
    ))
    return(
        list(
            gammaSum_tc = gammaSum_tc,
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
    S = 1,
    shuffle_bin_nSNPs = NULL,
    shuffle_bin_radius = 5000
) {
    ## 
    list_of_break_results <- lapply(1:S, function(x) NULL)
    nbreaks <- rep(0, S)            
    ## 
    if (iteration %in% shuffleHaplotypeIterations) {
        for(s in 1:S) {
            if (!is.null(shuffle_bin_radius)) {
                ## have previously (previous iteration) ended with define_breaks_to_consider
                ## now, load break_results
                ## NOTE - can be NULL, but unlikely
                load(file_break_results(tempdir = tempdir, regionName = regionName, s = s))
                ## 
                if (length(break_results) > 0) {
                    list_of_break_results[[s]] <- break_results                    
                    nbreaks[s] <- nrow(break_results)
                }
            } else {
                ## mimic format here as well
                start <- round(shuffle_bin_nSNPs / 2) + 1
                if (nGrids >= (start + shuffle_bin_nSNPs)) {
                    ##
                    left_break <- seq(start, nGrids - shuffle_bin_nSNPs, by = shuffle_bin_nSNPs)
                    nbreaks[s] <- length(left_break) - 1
                    list_of_break_results[[s]] <- cbind(
                        left_grid_break_0_based = left_break - 1,
                        left_grid_focal_0_based = left_break + (shuffle_bin_nSNPs - start),
                        right_grid_focal_0_based = left_break + (shuffle_bin_nSNPs - start + 1),
                        right_grid_break_0_based = left_break + shuffle_bin_nSNPs - 1
                    )
                }
            }
        }
    }
    return(
        list(
            list_of_break_results = list_of_break_results,
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



getTempMat <- function(nbreak, K, switchOrder) {
    ## determine the order of subsequent regions
    ##
    ## now - do the unravelling
    tempMat <- array(0, c(nbreak + 1,K))
    tempMat[1, ] <- 1:K
    currentState <- 1:K
    nextState <- 1:K
    for (iBreak in 1:(nbreak)) {
        ## how this works - first, start in your current state
        ## next, choose new state from switchOrder. copy that state, then remember where to go next
        currentState <- nextState
        for(k in 1:K) {
            ## first - the current value is wherever you are
            tempMat[iBreak + 1, k] <- switchOrder[iBreak, currentState[k]]
            nextState[k] <- tempMat[iBreak + 1, k]
        }
    }
    return(tempMat)
}


getBetterSwitchesSimple <- function(
    fromMat,
    nbreak,
    break_results,
    eHapsFuture_t,
    alphaMatFuture_t,
    grid,
    snps_in_grid_1_based,
    iteration = 1
) {
    ## to_run <- which(nbreaks > 0)
    ## s <- 1
    ## S <- length(nbreaks)
    K <- dim(eHapsFuture_t)[1]
    ## to_store <- array(0, length(to_run))
    ## list_of_whichIsBest <- lapply(1:S, function(x) return(NULL))
    ##for(i_s in 1:length(to_run)) {
    ## s <- to_run[i_s]
    ## fromMat <- list_of_fromMat[[s]]
    ## nbreak <- nbreaks[s]
    ## choose, in order, one which means fewest moves
    switchOrder <- determine_switch_order(fromMat, nbreak, K) 
    whichIsBest <- as.integer(apply(switchOrder,1,function(x) sum(x==(1:K)))==K)
    ## to_store[i_s] <- sum(whichIsBest!=1)
    ## list_of_whichIsBest[[i_s]]<- whichIsBest
    tempMat <- getTempMat(nbreak, K, switchOrder)
    ##
    ## do the shuffling
    ##
    ## 0-based start and end of the grids with the switches
    grid_starts <- c(0, break_results[, "left_grid_focal_0_based"] + 1)
    grid_ends <- c(break_results[, "left_grid_focal_0_based"], ncol(alphaMatFuture_t) - 1)
    for(iBreak in 1:(nbreak + 1)) {
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
    for(iBreak in 1:nbreak) {
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
    ##
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
    sigmaSum_m_unnormalized,
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
    S <- ncol(sigmaSum_m_unnormalized)
    ## if too few SNPs, do not bother
    if (nGrids < 5) {
        break_results <- NULL    
        for(s in 1:S) {        
            save(break_results, file = file_break_results(tempdir = tempdir, regionName = regionName, s = s))
        }
    }
    ## EM algorithm should be fine, won't get stuck this way
    for(s in 1:S) {
        ## these are defined and fine to work on
        sigmaSum_unnormalized <- sigmaSum_m_unnormalized[, s]
        ## generate smoothed rate
        smoothed_rate <- rcpp_make_smoothed_rate(
            sigma_rate = -log(sigmaSum_unnormalized) / grid_distances,
            L_grid = L_grid,
            shuffle_bin_radius = shuffle_bin_radius
        )
        ## normalize against max rate
        smoothed_rate <- smoothed_rate / 
            ((nGen * maxRate * 1/100/1000000))
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
        if (!plot_shuffle_haplotype_attempts) {
            save(break_results, file = file_break_results(tempdir = tempdir, regionName = regionName, s = s))
        } else {
            ## add in useful things for plot on subsequent iteration
            ##binSize <- 10000
            ##rate <- array(0, max(L)) ## this could be huge?
            ##for(i in 1:(nSNPs - 1)) {
            ##    rate[L[i]:L[i + 1]] <- (-log(sigmaSum_unnormalized[i]) / grid_distances[i])
            ## }
            x1_pre_exp <- -nGen * minRate * grid_distances/100/1000000
            x1 <- exp(x1_pre_exp) # lower
            x2_pre_exp <- -nGen * maxRate * grid_distances/100/1000000
            x2 <- exp(x2_pre_exp) # upper
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
                min_rate = (-x1_pre_exp / grid_distances),
                realized_rate = realized_rate,
                max_rate = (-x2_pre_exp / grid_distances),
                cumu_rate = cumu_rate
            )
            save(break_thresh, smoothed_rate, recomb_usage, break_results, file = file_break_results(tempdir = tempdir, regionName = regionName, s = s))
            ## save extra copy - useful for debugging etc
            ## save(break_thresh, smoothed_rate, recomb_usage, break_results, file = file_break_results(tempdir, regionName, iteration))
        }
    }
    return(NULL)
}



## smooth AT the point
## not BETWEEN the points
make_smoothed_cM_rate <- function(
    cM,
    L_grid,
    shuffle_bin_radius = 5000
) {
    ## 
    new_L <- L_grid[-length(L_grid)] + diff(L_grid)/2
    n <- length(new_L)
    new_L <- c(
        new_L[1] - (new_L[1] - L_grid[1]) / 2,
        new_L,
        new_L[n] + (L_grid[length(L_grid)] - new_L[n]) / 2
    )
    cM_smoothed <- rcpp_make_smoothed_rate(
        sigma_rate = cM,
        L_grid = new_L,
        shuffle_bin_radius = shuffle_bin_radius
    )
    return(cM_smoothed)
}

##
## returns one entry of average rate between pairs of SNPs
## value is average for middle point +/- shuffle_bin_radius
## 
make_smoothed_rate <- function(
    sigma_rate,
    L_grid,
    shuffle_bin_radius = 5000,
    verbose = FALSE,
    range = NA
) {
    nGrids <- length(L_grid)
    if (is.na(range)) {
        range <- 1:(nGrids - 1)
    }
    smoothed_rate <- array(0, nGrids - 1) ## between SNPs
    ##
    min_L_grid <- min(L_grid)
    max_L_grid <- max(L_grid)
    for(iGrid in range) {
        focal_point <- ((L_grid[iGrid] + L_grid[iGrid + 1]) / 2)
        total_bp_added <- 0
        if (verbose) {
            print_message(paste0("This is iGrid ", iGrid))
            print_message(paste0("This is between SNPs ", iGrid, " and ", iGrid + 1))
            print_message(paste0("The position of those two SNPs are", L[iGrid], " and ", L[iGrid + 1]))
            print_message(paste0("The focal point is:", focal_point))
            print_message(paste0("So, arguably, we want from:", focal_point - shuffle_bin_radius, " to ", focal_point + shuffle_bin_radius))
        }
        ## 
        ## left
        ## so re-call
        ## start with focal SNP, defined by position average of iGrid and iGrid + 1
        ## first region is iGrid_left
        iGrid_left <- iGrid
        bp_remaining <- shuffle_bin_radius
        bp_prev <- focal_point
        while ((0 < bp_remaining) & (0 < iGrid_left)) {
            bp_to_add <- (bp_prev - L_grid[iGrid_left])
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add <- bp_remaining
                bp_remaining <- 0
            } else {
                bp_remaining <- bp_remaining - bp_to_add
            }
            ## add bit
            if (verbose) {
                print_message(paste0("This is iGrid left :", iGrid_left))
                print_message(paste0("This is effectively from, inclusively:", bp_prev - bp_to_add, " to ", bp_prev - 1))
                print_message(paste0("We are adding this many bp: ", bp_to_add))
                print_message(paste0("Using rate: ", sigma_rate[iGrid_left]))
            }
            smoothed_rate[iGrid] <- smoothed_rate[iGrid] + bp_to_add * sigma_rate[iGrid_left]
            total_bp_added <- total_bp_added + bp_to_add
            bp_prev <- L_grid[iGrid_left]
            iGrid_left <- iGrid_left - 1                
        }
        ## right
        ## so recall. start with iGrid, then iGrid_right starts one to the right
        ## at first, difference between focal point (mid-way between iGrid and iGrid + 1)
        ## and iGrid + 1. then rate for that region is for iGrid (iGrid_right -1)
        iGrid_right <- iGrid + 1
        bp_remaining <- shuffle_bin_radius
        bp_prev <- focal_point
        while ((0 < bp_remaining) & (iGrid_right <= (nGrids))) {
            bp_to_add <- (L_grid[iGrid_right] - bp_prev)
            ##
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add <- bp_remaining
                bp_remaining <- 0 ## will end loop
            } else {
                bp_remaining <- bp_remaining - bp_to_add
            }
            ## add bit
            if (verbose) {
                print_message(paste0("This is iGrid right :", iGrid_right))
                print_message(paste0("This is effectively from, inclusively:", bp_prev, " to ", bp_prev + bp_to_add - 1))
                print_message(paste0("We are adding this many bp: ", bp_to_add))            
                print_message(paste0("Using rate: ", sigma_rate[iGrid_right - 1]))
            }
            smoothed_rate[iGrid] <- smoothed_rate[iGrid] +
                bp_to_add * sigma_rate[iGrid_right - 1]
            total_bp_added <- total_bp_added + bp_to_add
            bp_prev <- L_grid[iGrid_right] ## reset
            iGrid_right <- iGrid_right + 1
        }
        smoothed_rate[iGrid] <- smoothed_rate[iGrid] / total_bp_added
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
    whichIsBest = NULL,
    is_reference = FALSE,
    s = 1,
    S = 1
) {
    add_grey_background <- function(L_grid) {
        from <- floor(min(L_grid) / 1e6)
        to <- ceiling(max(L_grid) / 1e6)
        for(f in from:to) {
            abline(v = 1e6 * f, col = "grey")
        }
    }
    if (S == 1) {
        suffix <- ".png"
    } else {
        suffix <- paste0(".s.", s, ".png")
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
    if (is_reference) {
        ref_text <- ".ref"
    } else {
        ref_text <- ""
    }
    outname <- file.path(outputdir, "plots", paste0("shuffleHaplotypes", ref_text, ".", iteration, ".",regionName, suffix))
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


## make interim plots of things that might be useful to understand performance
## hapSumCurrent_t, regular and log10
## alphaMatCurrent_t
interim_plotter <- function(
    outputdir,
    regionName,
    iteration,
    L_grid,
    hapSumCurrent_tc,
    alphaMatCurrent_tc,
    sigmaCurrent_m,
    N,
    final_iteration = FALSE,
    is_reference = FALSE
) {
    ## easiest - just break up by S
    nGrids <- ncol(alphaMatCurrent_tc) + 1
    K <- nrow(alphaMatCurrent_tc)
    S <- dim(hapSumCurrent_tc)[3]
    for(s in 1:S) {
        alphaMatCurrent_t <- array(0, c(K, nGrids - 1))
        alphaMatCurrent_t[] <- alphaMatCurrent_tc[, , s]
        hapSumCurrent_t <- array(0, c(K, nGrids))
        hapSumCurrent_t[] <- hapSumCurrent_tc[, , s]
        plotHapSumCurrent_t(
            L_grid = L_grid,
            K = K,
            hapSumCurrent_t = hapSumCurrent_t,
            nGrids = nGrids,
            N = N,
            outputdir = outputdir,
            iteration = iteration,
            regionName = regionName,
            s = s,
            S = S,
            final_iteration = final_iteration,
            is_reference = is_reference
        )
        plotHapSumCurrent_t_log(
            L_grid = L_grid,
            K = K,
            hapSumCurrent_t = hapSumCurrent_t,
            nGrids = nGrids,
            N = N,
            outputdir = outputdir,
            regionName = regionName,
            iteration = iteration,
            s = s,
            S = S,
            final_iteration = final_iteration,
            is_reference = is_reference            
        )
        plotAlphaMatCurrent_t(
            L_grid = L_grid,
            alphaMatCurrent_t = alphaMatCurrent_t,
            sigmaCurrent = sigmaCurrent_m[, s],
            outputdir = outputdir,
            iteration = iteration,
            regionName = regionName,
            s = s,
            S = S,
            final_iteration = final_iteration,
            is_reference = is_reference            
        )
    }
    return(NULL)
}


interim_plot_name <- function(name, regionName, outputdir, s, S, iteration, final_iteration, is_reference, what = "") {
    if (final_iteration) {
        it_name <- ""
    } else {
        it_name <- paste0(".iteration.", iteration)
    }
    if (S > 1) {
        suffix <- ".png"
    } else {
        suffix <- paste0(".s.", s, ".png")
    }
    if (is_reference) {
        name <- paste0("ref.", name)
    }
    if (nchar(what) > 0) {
        what <- paste0(".", what)
    }
    outname <- file.path(outputdir, "plots", paste0(name, ".", regionName, it_name, what, suffix))
    return(outname)
}

## plot hapSumCurrent along genome
plotHapSumCurrent_t <- function(
    outname,
    L_grid,
    K,
    hapSumCurrent_t,
    nGrids,
    N,
    outputdir,
    regionName,
    iteration,
    s,
    S,
    final_iteration,
    is_reference
) {
    ##
    outname <- interim_plot_name("hapSum", regionName, outputdir, s, S, iteration, final_iteration, is_reference)
    ## 
    width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)    
    png(outname, height = 10, width = width, res = 100, units = "in")    
    colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    nCols <- length(colStore)
    sum <- array(0, nGrids)
    xlim <- range(L_grid)
    ylim <- c(0, 1)
    ## OK so if there are grids, use the grid points
    plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE)
    x <- c(L_grid[1], L_grid, L_grid[length(L_grid):1])
    m <- array(0, c(nGrids, K + 1))
    for(i in 1:K) {
        m[, i + 1] <- m[, i] + hapSumCurrent_t[i, , drop = FALSE] / N
    }
    for(i in K:1) {
        polygon(
            x = x, y = c(m[1, i], m[, i + 1], m[nGrids:1, i]),
            xlim = xlim, ylim = ylim, col = colStore[(i %% nCols) + 1]
        )
    }
    dev.off()
}



## plot hapSumCurrent along genome
plotHapSumCurrent_t_log <- function(
    L_grid,
    K,
    hapSumCurrent_t,
    nGrids,
    N,
    outputdir,
    regionName,
    iteration,
    s,
    S,
    final_iteration,
    is_reference
) {
    outname <- interim_plot_name("hapSum_log", regionName, outputdir, s, S, iteration, final_iteration, is_reference)
    ## 
    width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)
    png(outname, height = 10, width = width, res = 100, units = "in")
    main <- "log10 of average haplotype usage vs physical position"
    colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    nCols <- length(colStore)
    sum <- array(0, nGrids)
    xlim <- range(L_grid) / 1e6    
    ylim <- c(log10(max(1, min(hapSumCurrent_t))), log10(max(hapSumCurrent_t)))
    if (sum(hapSumCurrent_t) == 0) {
        stop("Something has done wrong and an ampty hapSumCurrent_t has been passed to plotHapSumCurrent_t_log")
    }
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = main, xlab = "Physical position (Mbp)", ylab = "log10 Haplotype usage")
    axis(1)
    axis(2)
    x <- L_grid / 1e6
    for(k in 1:K) {
        points(x = x, y = log10(hapSumCurrent_t[k, ]), col = colStore[(k %% nCols) + 1], type = "l")
    }
    for(i in floor(ylim[1]):ceiling(ylim[2])) {
        abline(h = i, col = "grey")
    }
    abline(h = log10(N / K), col = "red")
    dev.off()
    ##    system(paste0("rsync -av '", outname, "' florence:~/"))
}



plotAlphaMatCurrent_t <- function(L_grid, alphaMatCurrent_t, sigmaCurrent, outputdir, iteration, regionName, s, S, final_iteration, is_reference) {
    if (S > 1) {
        suffix <- ".png"
    } else {
        suffix <- paste0(".s.", s, ".png")
    }
    nGrids <- ncol(alphaMatCurrent_t)
    K <- nrow(alphaMatCurrent_t)
    ## two views - proportional, total?
    xleft <- (L_grid)[-length(L_grid)]
    xright <- (L_grid)[-1]
    xlim <- range(L_grid)
    ##
    for(i_what in 1:2) {
        ## make barplot type thing
        if (i_what == 2) {
            ## normalize by sigmaCurrent
            main <- "alphaMatCurrent_t P(q_t = k, I_t = 1)"
            x <- alphaMatCurrent_t
            for(k in 1:K) {
                x[k, ] <- alphaMatCurrent_t[k, ] * (1 - sigmaCurrent)
            }
            fbd_store <- list(list(gammaK_t = x))
            what <- "normalized"
        } else {
            main <- "alphaMatCurrent_t P(q_t | I_t = 1)"
            fbd_store <- list(list(gammaK_t = alphaMatCurrent_t))
            what <- "all"
        }
        outname <- interim_plot_name("alphaMat", regionName, outputdir, s, S, iteration, final_iteration, is_reference, what)
        width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)
        png(outname, height = 10, width = width, res = 100, units = "in")
        plot_fbd_store(fbd_store, xleft, xright, xlim, main = main)
        dev.off()
    }
}



apply_better_switches_if_appropriate <- function(
    list_of_fromMat,
    nbreaks,
    list_of_break_results,
    eHapsCurrent_tc,
    alphaMatCurrent_tc,
    grid,
    iteration,
    snps_in_grid_1_based,
    tempdir,
    regionName,
    grid_distances,
    L_grid,
    outputdir,
    is_reference,
    plot_shuffle_haplotype_attempts
) {
    ## only do for s that suit it
    K <- dim(eHapsCurrent_tc)[1]
    S <- dim(eHapsCurrent_tc)[3]
    ## 
    list_of_whichIsBest <- lapply(1:S, function(s) return(NA))
    for(s in 1:S) {
        if (nbreaks[s] > 0) {
            out <- getBetterSwitchesSimple(
                fromMat = list_of_fromMat[[s]],
                nbreak = nbreaks[s],
                break_results = list_of_break_results[[s]],
                eHapsFuture_t = eHapsCurrent_tc[, , s],
                alphaMatFuture_t = alphaMatCurrent_tc[, , s],
                grid = grid,
                iteration = iteration,
                snps_in_grid_1_based = snps_in_grid_1_based
            )
            eHapsCurrent_tc[, , s] <- out$eHapsFuture_t
            alphaMatCurrent_tc[, , s] <- out$alphaMatFuture_t
            list_of_whichIsBest[[s]] <- out$whichIsBest
            rm(out)
        }
    }
    if (sum(nbreaks) > 0) {
        to_store <- sapply(list_of_whichIsBest, function(whichIsBest) sum(whichIsBest!=1, na.rm = TRUE))
        print_message(paste0(
            "Shuffle haplotypes - Iteration ", iteration, " - ", 
            "change on average ", round(sum(to_store) / S, 1), " intervals out of ", sum(nbreaks) / S, " considered"
        ))
    }
    if (plot_shuffle_haplotype_attempts && sum(nbreaks) > 0) {
        ## a single run, everything inside
        load(file = file_fbdStore(tempdir, regionName, iteration)) ## to load fbdStore
        for(s in 1:S) {
            if (nbreaks[s] > 0) {
                plot_attempt_to_find_shuffles(grid_distances = grid_distances, L_grid = L_grid, fbd_store = list_of_fbd_store[[s]], tempdir = tempdir, outputdir = outputdir, regionName = regionName, iteration = iteration, whichIsBest = list_of_whichIsBest[[s]], is_reference = is_reference, s = s, S = S)
            }
        }
    }
    return(
        list(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrent_tc = alphaMatCurrent_tc            
        )
    )
}


findRecombinedReadsPerSample <- function(
    gammaK_t,
    eHapsCurrent_t,
    K,
    L,
    iSample,
    sampleReads,
    tempdir,
    regionName,
    grid,
    verbose = TRUE,
    method = "diploid",
    pRgivenH1_m = NULL,
    pRgivenH2_m = NULL,
    srp = NULL
) {
    K <- dim(eHapsCurrent_t)[1]
    ## needs a full run
    ## only do for some - need at least 3 SNPs to consider
    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    w <- w[w != 1 & w != length(w)]
    count <- 0
    if (length(w) > 0) {
        for (w1 in w) {
            out <- split_a_read(
                sampleReads = sampleReads,
                read_to_split = w1,
                gammaK_t = gammaK_t,
                L = L,
                eHapsCurrent_t = eHapsCurrent_t,
                K = K,
                grid = grid,
                method = method,
                pRgivenH1_m = pRgivenH1_m,
                pRgivenH2_m = pRgivenH2_m,
                srp = srp
            )
            sampleReads <- out$sampleReads
            pRgivenH1_m <- out$pRgivenH1_m
            pRgivenH2_m <- out$pRgivenH2_m
            srp <- srp
            count <- count + as.integer(out$did_split)
        } # end of loop on reads
        new_order <- order(unlist(lapply(sampleReads,function(x) x[[2]])))        
        sampleReads <- sampleReads[new_order]
        save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName), compress = FALSE)
        if (verbose) {
            print_message(paste0(
                "sample ", iSample, " readsSplit ", count, " readsTotal ", length(sampleReads)
            ))
        }
        if (method == "pseudoHaploid") {
            ## randomize those of split reads
            srp <- srp[new_order]
            pRgivenH1_m <- pRgivenH1_m[new_order, , drop = FALSE]
            pRgivenH2_m <- pRgivenH2_m[new_order, , drop = FALSE]
            save(
                srp, pRgivenH1_m, pRgivenH2_m,
                file = file_sampleProbs(tempdir, iSample, regionName)
            )
        }        
    }
    return(
        list(
            readsSplit = count,
            readsTotal = length(sampleReads)
        )
    )
}


split_a_read <- function(
    sampleReads,
    read_to_split,
    gammaK_t,
    L,
    eHapsCurrent_t,
    K,
    grid,
    method,
    pRgivenH1_m,
    pRgivenH2_m,
    srp
) {

    did_split <- FALSE
    sampleRead <- sampleReads[[read_to_split]]

    ## get the likely for and after states
    ## by checking, a few reads up and downstream, what changes
    ## note - probably want to change this to per-base, not per-read!
    startRead <- max(1, read_to_split - 10)
    endRead <- min(length(sampleReads), read_to_split + 10)
    ##
    startGammaGrid <- sampleReads[[startRead]][[2]] + 1
    endGammaGrid <- sampleReads[[endRead]][[2]] + 1

    change <- apply(gammaK_t[, c(startGammaGrid, endGammaGrid)], 1, diff)
    ## from is the one that drops the most
    ## to is the one that gains the most
    from <- which.min(change)
    to <- which.max(change)

    ## try a break between all SNPs - take the best one
    y <- eHapsCurrent_t[, sampleRead[[4]] + 1]

    bqProbs <- convertScaledBQtoProbs(sampleRead[[3]])
    g1 <- bqProbs[, 2] * y[from, ] +
          bqProbs[, 1] * (1 - y[from, ])
    g2 <- bqProbs[, 2] * y[to, ] +
          bqProbs[, 1] * (1 - y[to, ])

    ## sum probabilities
    ## best split is between SNPs
    z <- sapply(1:(sampleRead[[1]]),function(a) {
        x1 <- sum(g1[1:a])
        x2 <- sum(g2[(a+1):(sampleRead[[1]] + 1)])
        sum(x1 + x2)
    })
    best_split <- which.max(z) ## between this and + 1

    u <- log(
        split_function(eHapsCurrent_t, sampleRead, K)[ K + 1]
    )

    if (max(z)>u) {
        new_read_1 <- get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, 1, best_split, L, grid
        )
        new_read_2 <- get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, best_split + 1, sampleRead[[1]] + 1, L, grid
        )
        sampleReads[[read_to_split]] <- new_read_1
        sampleReads <- append(
            sampleReads,
            list(new_read_2)
        )
        did_split <- TRUE
        if (method == "pseudoHaploid") {
            pRgivenH1_m[read_to_split, ] <- runif(ncol(pRgivenH1_m))
            pRgivenH2_m[read_to_split, ] <- runif(ncol(pRgivenH2_m))
            srp[read_to_split] <- sampleReads[[read_to_split]][[2]]
            pRgivenH1_m <- rbind(pRgivenH1_m, runif(ncol(pRgivenH1_m)))
            pRgivenH2_m <- rbind(pRgivenH2_m, runif(ncol(pRgivenH2_m)))
            srp <- c(srp, sampleReads[[length(sampleReads)]][[2]])
        }        
    }

    return(
        list(
            sampleReads = sampleReads,
            did_split = did_split,
            pRgivenH1_m = pRgivenH1_m,
            pRgivenH2_m = pRgivenH2_m,
            srp = srp            
        )
    )
}
