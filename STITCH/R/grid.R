## assign each of the positions to a grid
## input is numbers, e.g. 3, 5, 10, 15
## and a windowSize, like 5
## output is 1-based on grid coordinates, like
## 1-5 -> 0, 6-10 -> 1, etc
## remove holes, sigma will be made able to handle it with bounding

## method = "normal",
##    smoothed_rate_cM = NULL,
##    desired_gridWindowSize = NULL

#' @export
assign_positions_to_grid <- function(
    L,
    gridWindowSize = NA,
    cM = NULL,
    grid32 = FALSE
) {
    cM_grid <- NULL
    if (grid32) {
        nSNPs <- length(L)
        grid <- rep(0:(nSNPs - 1), each = 32, length.out = length(L))
        ## this is average in each
        L_grid <- array(NA, max(grid) + 1)
        for(i in min(grid):max(grid)) {
            s <- 32 * i + 1
            e <- min(32 * (i + 1), nSNPs)
            L_grid[i + 1] <- round(mean(L[s:e]))
        }
        nGrids <- length(L_grid)
        grid_distances <- diff(L_grid)
        ## L_grid is average within the grid
        snps_in_grid_1_based <- cbind(
            snps_start = match(unique(grid), grid),
            snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
        )
    } else if (is.na(gridWindowSize)) {
        ## simple case - no gridding
        grid <- 0:(length(L) - 1)
        grid_distances <- diff(L)
        L_grid <- L
        nGrids <- length(L)
        snps_in_grid_1_based <- cbind(
            snps_start = 1:nGrids,
            snps_end = 1:nGrids
        )
        cM_grid <- cM
    } else {
        ## if no genetic map, make based on physical distance only
        if (is.null(cM)) {
            grid <- ceiling(L / gridWindowSize)
            grid <- grid - min(grid)
            ## for L_grid, get first mid-point
            L_grid_start <- gridWindowSize * (ceiling(L[1] / gridWindowSize) - 0.5)
            grid_distances <- diff(unique(grid)) * gridWindowSize
            L_grid <- L_grid_start + c(0, cumsum(grid_distances))
            grid <- match(grid, unique(grid)) - 1
            nGrids <- length(grid_distances) + 1
            if ((length(L_grid) - length(grid_distances)) != 1) {
                stop("An error has been made assigning SNP positions to grid. Please report this")
            }
            ## this allows one to go from grid to start and end of SNPs in that grid
            ## so e.g. the fifth grid with 0-based index 4
            ## has 1-based SNPs starting from
            ## snps_in_grid_1_based[4 + 1, "grid_starts"]
            ## to
            ## snps_in_grid_1_based[4 + 1, "grid_ends"]        
            snps_in_grid_1_based <- cbind(
                snps_start = match(unique(grid), grid),
                snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
            )
        } else {
            ## here we have a genetic map!
            out <- attempt_to_better_grid(
                cM = cM,
                gridWindowSize = gridWindowSize, ## why the name change
                L = L
            )
            return(out)
        }
    }
    return(
        list(
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based,
            cM_grid = cM_grid
        )
    )
    return(to_out)
}


## alternative approach - approximate!
attempt_to_better_grid <- function(
    gridWindowSize,
    L,
    smoothed_rate_cM = NULL,
    cM = NULL
) {
    if (is.null(smoothed_rate_cM)) {
        ## convert cM to smoothed rate
        smoothed_rate_cM <- 1e6 * diff(cM) / diff(L)
    }
    ## note - smoothed_rate_cM is the genetic distance between points
    smoothed_rate_cM[is.na(smoothed_rate_cM)] <- 0
    dist <- L[length(L)] - L[1]    
    total_smoothed_rate_cM <- sum(smoothed_rate_cM)
    nSNP <- length(smoothed_rate_cM) + 1
    nGrids <- round(dist / gridWindowSize) ## approximate number of grids we want
    average_smoothed_rate_cM_per_grid <- total_smoothed_rate_cM / nGrids
    ## so want to divide to average this! how to do...
    running_total <- 0
    iGrid <- 0
    grid <- array(NA, nSNP)
    grid[1] <- 0 ## 0-based
    L_grid <- array(NA, nSNP) ## will be shortened!
    first_snp_in_this_grid <- 1
    for(iSNP in 2:nSNP) {
        running_total <- running_total + smoothed_rate_cM[iSNP - 1] ##
        ## if it exceeds - break
        if ((running_total > average_smoothed_rate_cM_per_grid)) {
            ## this is now the first SNP of it's own grid
            running_total <- 0
            w <- first_snp_in_this_grid:(iSNP - 1)
            L_grid[iGrid + 1] <- mean(L[w])
            if (iSNP != nSNP) {
                iGrid <- iGrid + 1
            }
            ## Now - set L_grid is average of those in the grid
            first_snp_in_this_grid <- iSNP
        } else if (iSNP == nSNP) {
            ## if it is just done, end in its own way
            ## for this grid, include this particular SNP
            w <- first_snp_in_this_grid:(iSNP)
            L_grid[iGrid + 1] <- mean(L[w])
        }
        grid[iSNP] <- iGrid
    }
    L_grid <- L_grid[1:(iGrid + 1)]
    nGrids <- length(L_grid)
    grid_distances <- diff(L_grid)
    ## L_grid is average within the grid
    snps_in_grid_1_based <- cbind(
        snps_start = match(unique(grid), grid),
        snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
    )
    ##
    ## now - give back gridded genetic map, to use for sigmaCurrent later
    ##
    temp_genetic_map <- array(NA, c(length(L), 3))
    colnames(temp_genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")    
    temp_genetic_map[, "position"] <- L
    temp_genetic_map[, "Genetic_Map.cM."] <- cM
    cM_grid <- match_genetic_map_to_L(
        genetic_map = temp_genetic_map,
        L = L_grid,
        expRate = NA
    )  ## expRate should not matter - is in interior
    return(
        list(
            grid = grid,
            L_grid = L_grid,
            nGrids = nGrids,
            grid_distances = grid_distances,
            snps_in_grid_1_based = snps_in_grid_1_based,
            cM_grid = cM_grid
        )
    )
}



## L_grid is distance between average grid poins
assign_grid_variables <- function(grid, L) {
    ## distances
    L_grid <- array(NA, length(L)) ## will be shortened!
    iGrid <- 0
    first_snp_in_this_grid <- 1
    for(iSNP in 2:length(L)) {
        if ((grid[iSNP] > iGrid) | (iSNP == length(L))) {
            print(iSNP)
            w <- first_snp_in_this_grid:(iSNP - 1)
            L_grid[iGrid] <- mean(L[w])
            iGrid <- iGrid + 1
            first_snp_in_this_grid <- iSNP + 1
        }
    }
    ## final one
    L_grid <- L_grid[1:(iGrid - 1)]
    if (sum(is.na(L_grid)) > 0) {
        print(sum(is.na(L_grid)))
        stop("Something has gone wrong with STITCH internals. Please report this")
    }
    grid_distances <- diff(L_grid)    
    nGrids <- length(grid_distances) + 1
    ## this allows one to go from grid to start and end of SNPs in that grid
    ## so e.g. the fifth grid with 0-based index 4
    ## has 1-based SNPs starting from
    ## snps_in_grid_1_based[4 + 1, "grid_starts"]
    ## to
    ## snps_in_grid_1_based[4 + 1, "grid_ends"]        
    snps_in_grid_1_based <- cbind(
        snps_start = match(unique(grid), grid),
        snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
    )
    return(
        list(
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based
        )
    )
}




get_new_grid <- function(nGen, sigmaCurrent, L, desired_gridWindowSize, shuffle_bin_radius = 5000) {
    grid_distances <- diff(L)
    recombRate_cM_per_Mb <- (-(log(sigmaCurrent) / grid_distances) * 1e6 * 100) / nGen
    recombRate_cM <- (-log(sigmaCurrent) * 100) / nGen
    ##
    nGrids <- length(L)
    smoothed_rate_cM_per_Mb <- rcpp_make_smoothed_rate(
        sigmaSum = NA, ## this doesn't matter!
        sigma_rate = recombRate_cM_per_Mb,
        L_grid = L,
        grid_distances = grid_distances, ## when no grid, this is diff(L)
        nGrids = nGrids, ## when no grid, this equals nSNPs
        shuffle_bin_radius = shuffle_bin_radius
    ) / (shuffle_bin_radius * 2)
    ##
    smoothed_rate_cM <- smoothed_rate_cM_per_Mb * grid_distances / 1e6
    ##  
    out <- attempt_to_better_grid(
        smoothed_rate_cM = smoothed_rate_cM,
        desired_gridWindowSize = desired_gridWindowSize,
        L = L
    )
    out2 <- make_new_sigmaCurrent_and_alphaMatCurrent_from_new_grid(
        snps_in_grid_1_based = out$snps_in_grid_1_based,
        sigmaCurrent = sigmaCurrent,
        alphaMatCurrent = alphaMatCurrent,
        nGen = nGen,
        recombRate_cM = recombRate_cM
    )
    out$new_sigmaCurrent <- out2$new_sigmaCurrent
    out$new_alphaMatCurrent <- out2$new_alphaMatCurrent
    return(out)
}


make_new_sigmaCurrent_and_alphaMatCurrent_from_new_grid <- function(
    snps_in_grid_1_based,
    sigmaCurrent,
    alphaMatCurrent,
    nGen,
    recombRate_cM
) {
    ##
    nGrids <- nrow(snps_in_grid_1_based)
    ##
    new_recombRate_cM <- array(NA, nGrids - 1)
    new_alphaMatCurrent <- array(NA, c(nGrids - 1, ncol(alphaMatCurrent)))
    iGrid <- 1
    s1 <- snps_in_grid_1_based[iGrid, "snps_start"]
    e1 <- snps_in_grid_1_based[iGrid, "snps_end"]
    ## here, recombRate_cM is the difference in cM between points
    ## so e.g. if the first grid is SNPs 1 -> 5
    ## we want the average cM in that group
    ## so take differences vs start of 0
    ## i.e. sum(c(0, then 1->4 on grid))
    if (e1 == 1) {
        prevRate <- 0
        a <- ncol(alphaMatCurrent) ## might be K or K2 etc
        prevAlpha <- rep(1 / a, a)
    } else {
        prevRate <- sum(c(0, recombRate_cM[s1:(e1 - 1)]))
        prevAlpha <- colSums(alphaMatCurrent[s1:(e1 - 1), , drop = FALSE]) / (e1 - s1 + 1)
    }
    ##
    for(iGrid in 2:(nGrids)) {
        ## get difference in cM between the two
        s2 <- snps_in_grid_1_based[iGrid, "snps_start"] - 1
        e2 <- snps_in_grid_1_based[iGrid, "snps_end"] - 1
        ## want half of previous plus half of this?
        curRate <- sum(recombRate_cM[s2:e2])
        new_recombRate_cM[iGrid - 1] <- (curRate + prevRate) / 2
        prevRate <- curRate
        ## 
        curAlpha <- colSums(alphaMatCurrent[s2:e2, , drop = FALSE]) / (e2 - s2 + 1)
        new_alphaMatCurrent[iGrid - 1, ] <- (curAlpha + prevAlpha) / 2
        prevAlpha <- curAlpha
    }
    ## now convert back to sigmaCurrent units!
    new_sigmaCurrent <- exp(-nGen * new_recombRate_cM / 100)
    return(
        list(
            new_sigmaCurrent = new_sigmaCurrent,
            new_alphaMatCurrent = new_alphaMatCurrent
        )
    )
}


