## in a simple refill, look every 100 SNP. in each one with an average below a threshold, sample from other haplotypes to re-copy
## across all haplotypes, when this happen, inject noise
refillSimple <- function(
    hapSum,
    nGrids,
    K,
    gammaSum,
    N,
    L_grid,
    grid,
    distance_between_check = 5000 ## check every 5000 bp
) {
    nSNPs <- nrow(gammaSum)
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
        avHapSumInBlock[i, ] <- colSums(hapSum[region_starts[i]:region_ends[i], , drop = FALSE]) 
    }
    hapFreqMin <- N / K ## less than average - aggressive!
    replaceBlock <- avHapSumInBlock < (hapFreqMin)
    ## 
    ## for each k, for each region, fill in
    ## within each continuous region, fill in with respect to frequencies of all other haplotypes
    ##
    k_to_replace <- which(colSums(replaceBlock) > 0)
    ever_changed <- array(FALSE, nrow(gammaSum))
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
                gammaSum[snps_to_replace, k] <- gammaSum[snps_to_replace, replacement]
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
    return(list(gammaSum = gammaSum, ever_changed = ever_changed))
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
                nbreaks <- length(left_breaks) - 1
                break_results <- cbind(
                    left_break = left_break,
                    left_focal_SNP = left_break + (shuffle_bin_nSNPs - start),
                    right_focal_SNP = left_break + (shuffle_bin_nSNPs - start + 1),
                    right_break = left_break + shuffle_bin_nSNPs - 1
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
    eHapsFuture,
    alphaMatFuture,
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
    switchCur <- 1:K
    rStart <- c(1, break_results[, "left_focal"])
    rEnd <- c(break_results[, "right_focal"], nrow(eHapsFuture))
    ## add in 
    for(iBreak in 1:(nbreaks + 1)) {
        w <- rStart[iBreak]:rEnd[iBreak]
        permL <- tempMat[iBreak, ]
        eHapsFuture[w, ] <- eHapsFuture[w, permL]
        if (iBreak == (nbreaks + 1)) {
            w <- rStart[iBreak]:(rEnd[iBreak] - 1)
        }
        alphaMatFuture[w, ] <- alphaMatFuture[w, permL]
    }
    ## add in noise around the breaks
    for(iBreak in 1:nbreaks) {
        if(whichIsBest[iBreak]==0) {
            ## add in some noise around break            
            ## progressively more in the middle
            s <- break_results[iBreak, "left_break"]
            e <- break_results[iBreak, "right_break"]
            norm_component <- abs(1 - seq(0, 2, length = (e - s) + 1))
            w <- s:e
            ## fill in with more noise closer to the break
            for(ii in 1:length(w)) {
                i <- (w)[ii]
                eHapsFuture[i, ] <-
                    norm_component[ii] * eHapsFuture[i, ] +
                    (1 - norm_component[ii]) * runif(K)
            }
            ## just reset these, can be re-determined pretty quickly
            alphaMatFuture[w, ] <- matrix(1 / K, ncol = K, nrow = length(w))
        }
    }
    return(
        list(
            eHapsFuture = eHapsFuture,
            alphaMatFuture = alphaMatFuture
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
## then smooth rate over those SNPs is betwene round((a + b) / 2) minus 2000 and plus 2000
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
    shuffle_bin_radius = 2000,
    plot_shuffle_haplotype_attempts = FALSE
) {
    ## if too few SNPs, do not bother
    ## EM algorithm should be fine, won't get stuck this way
    if ((nGrids < 100)) {
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
            nGrids = nGrids
        )
        break_results <- out$results
        break_thresh <- out$thresh
        ## NULL return unlikely but probably probable
        if (plot_shuffle_haplotype_attempts == FALSE) {
            save(break_results, file = file_break_results(tempdir, regionName))
        } else {
            ## add in useful things for plot on subsequent iteration
            binSize <- 10000
            rate <- array(0, 5000000)
            for(i in 1:(nSNPs - 1)) {
                rate[L[i]:L[i + 1]] <- (-log(sigmaSum_unnormalized[i]) / grid_distances[i])
            }
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
        }
    }
    return(NULL)
}

make_smoothed_rate <- function(sigmaSum_unnormalized, sigma_rate, L_grid, grid_distances, nGrids, shuffle_bin_radius = 5000) {
    smoothed_rate <- array(0, nGrids - 1) ## between SNPs
    ##
    min_L_grid <- min(L_grid)
    max_L_grid <- max(L_grid)
    for(iSNP in 1:(nGrids - 1)) {
        focal_point <- floor((L_grid[iSNP] + L_grid[iSNP + 1]) / 2)        
        if (
            min_L_grid < (focal_point - shuffle_bin_radius) &
            (focal_point + shuffle_bin_radius) < max_L_grid
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
    min_breaks = 10,
    max_breaks = Inf
) {
    smoothed_rate[1] <- NA
    smoothed_rate[length(smoothed_rate)] <- NA    
    ## now, those > 1 are great candidates
    ## alternatively, if nGen not quite right (or rates not changed!)
    ## anything really "peaky"
    ## remember, nothing forces it to make the change later, so can be really liberal here
    ## so choose as threshold lower of 1, or 50th percentile
    ## boundaries where < thresh or start rising 5+ SNPs
    thresh <-  min(1, quantile(smoothed_rate[smoothed_rate != 0], probs = 0.95, na.rm = TRUE))
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
        snp_best <- best[iBest]
        ## only consider if have >1 SNP either side
        if (sum(available[snp_best + -1:1]) == 3) {
            ##
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
            results <- rbind(
                results,
                c(snp_left, snp_best, snp_best + 1, snp_right)
            )
            ## nuke out (do not consider) at least
            ##   those SNPs under consideration
            ##   50 SNPs either way
            nuke_left <- max(1, min(snp_left, snp_best - 50))
            nuke_right <- min(nGrids - 2, max(snp_right, snp_best + 50))
            available[nuke_left:nuke_right] <- FALSE
            ideal[nuke_left:nuke_right] <- FALSE
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
        colnames(results) <- c("left_break", "left_focal", "right_focal", "right_break")
        results <- results[order(results[, "left_focal"]), ]        
    }
    return(
        list(
            results = results,
            thresh = thresh
        )
    )
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
        ## do not continue if not avaialble i.e. in previous region        
        if (available[snp_consider + (-1) * mult] == FALSE) {
            done <- 1
        }
        ## do not continue if 300% of minimum
        ## also, re-set to minimum value
        if ((3 * val_min) < val_cur) {
            done <- 1
            snp_consider <- snp_min
        }
        ## stop if below threshold AND not still decreasing over 5 SNP
        if ((val_cur < thresh) & (val_prev < val_cur)) {
            done <- 1
        }
        ## stop if not decreasing over 10 SNPs
    }
    return(snp_consider)
}



plot_attempt_to_find_shuffles <- function(
    grid_distances,
    L,
    fbd_store,
    tempdir,
    outputdir,
    regionName,
    iteration
) {
    add_grey_background <- function(L) {
        from <- floor(min(L) / 1e6)
        to <- ceiling(max(L) / 1e6)
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
    xlim <- c(head(L)[1], tail(L)[1])    
    ## do plot here
    outname <- file.path(outputdir, "plots", paste0("shuffleHaplotypes.", iteration, ".",regionName,".png"))
    ## make a 5 Mbp segment 60 wide. then bound up and down at 20 and 200
    width <- min(max(20, (L[length(L)] - L[1]) / 1e6 * 12), 200)
    png(outname, height = 30, width = width, res = 100, units = "in")
    par(mfrow = c(3, 1))    
    x <- L[-1] - grid_distances ## useful
    xleft <- L[-length(L)]
    xright <- L[-1]
    ##
    ## 0) find peaks to check
    ##
    ylim <- c(0, max(break_thresh, max(smoothed_rate, na.rm = TRUE)))
    plot(x = 0, y = 0, xlab = "Physical position", ylab = "Rate", main = "Location of shuffles to check", ylim = ylim, xlim = xlim)
    add_grey_background(L)
    lines(x = x, y = smoothed_rate, lwd = 2)
    for(iBreak in 1:nrow(break_results)) {
        abline(v = L[break_results[iBreak, "left_break"]], col = "red")
        abline(v = L[break_results[iBreak, "left_focal"]], col = "purple")
        abline(v = L[break_results[iBreak, "right_break"]], col = "red")
        rect(
            xleft = L[break_results[iBreak, "left_break"]],
            xright = L[break_results[iBreak, "right_break"]],
            ybottom = -1, ytop = 0,
            col = "green"
        )
    }
    abline(h = break_thresh, col = "black", lwd = 2)
    ##
    ## 1) plot fine and "coarse" scale recombination
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
    add_grey_background(L)
    lines(x = x, y = log10(recomb_usage[, "min_rate"]), ylim = ylim, col = "blue")
    lines(x = x, y = log10(recomb_usage[, "realized_rate"]), ylim = ylim, col = "green")
    ## lines(x = seq(binSize / 2 + 1, 5000000, binSize), y = log10(rate_binned), ylim = ylim, col = "purple")
    lines(x = x, y = log10(recomb_usage[, "max_rate"]), ylim = ylim, col = "red")
    lines(x = x, y = cumu_rate_scaled, ylim = ylim, col = "orange")
    ## then add all the switches?
    ##
    ## 2) plot individuals with their switches
    ##
    NN <- length(fbd_store)
    ylim <- c(1, NN + 1)
    plot(x = 0, y = 0, ylim = ylim, xlim = xlim, xlab = "Physical position", ylab = "Sample", main = "Haplotype usage per-sample", col = "white")
    cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)
    for(iSample in 1:NN) {
        ## plot each of them, as rectangle?
        R <- fbd_store[[iSample]]$gammaK_t
        ## should sum to 1, plot rectangles of each
        ybottom <- iSample + array(0, ncol(R) - 1)
        ytop <- iSample + array(0, ncol(R) - 1)
        for(k in 1:nrow(R)) {
            ytop <- ytop + R[k, -ncol(R)]
            rect(
                xleft = xleft,
                xright = xright,
                ybottom = ybottom,
                ytop = ytop,
                col = cbPalette[k],
                border = NA
            )
            ybottom <- ytop
        }
    }
    add_grey_background(L)
    legend("topright", paste0("anc_hap=", 1:nrow(R)), col = cbPalette[1:nrow(R)], lwd = 2)
    dev.off()
}




alpha_col <- function(col, alpha) {
    x <- col2rgb(col) / 255
    return(rgb(x["red", 1], x["green", 1], x["blue", 1], alpha = alpha)    )
}
