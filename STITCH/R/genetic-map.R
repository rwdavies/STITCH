fill_in_genetic_map_cm_column <- function(genetic_map) {
    genetic_map[1, "Genetic_Map.cM."] <- 0
    for(iRow in 2:nrow(genetic_map)) {
        distance <-
            genetic_map[iRow, "position"] -
            genetic_map[iRow - 1, "position"]
        rate <- genetic_map[iRow - 1, "COMBINED_rate.cM.Mb."]
        genetic_map[iRow, "Genetic_Map.cM."] <-
            genetic_map[iRow - 1, "Genetic_Map.cM."] + 
            (distance / 1e6) * rate
    }
    return(genetic_map)
}


validate_genetic_map <- function(genetic_map, verbose = TRUE) {
    if (ncol(genetic_map) != 3) {
        stop(paste0("Provided reference genetic map has ", ncol(genetic_map), " columns, where 3 are expected, with columns that provide position (in bp), combined genetic rate map (in cM/Mbp), and genetic map (in cM)"))
    }
    original_colnames_genetic_map <- colnames(genetic_map)
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    ## do not allow any NAs or non-numeric characters
    if (sum(is.na(genetic_map) | is.na(is.na(genetic_map * 2))) > 0) {
        stop(paste0("Provided reference genetic map contains non-numeric entries or NA entries. Please ensure all entries are numeric and non-empty"))
    }
    ## 
    check_col_3 <- fill_in_genetic_map_cm_column(genetic_map)[, 3]
    x <- (genetic_map[, 3] - check_col_3) > 1e-8
    if (sum(x) > 0) {
        if (verbose) {
            print(genetic_map[c(max(1, x - 1):min(x + 1, nrow(genetic_map))) , ])
        }
        stop(paste0("Error interpreting genetic map around row ", which.max(x), ". Provided genetic map does not satisfy expectations. See above print. Let g be the genetic map with 3 columns. In 1-based coordinates, recall that g[iRow, 3] must equal g[iRow - 1, 3] + ( (g[iRow, 1] - g[iRow - 1, 1]) * g[iRow - 1, 2])"))
    }
    return(NULL)
}

get_and_validate_genetic_map <- function(genetic_map_file) {
    ## 
    genetic_map <- read.table(reference_genetic_map_file, header = TRUE)
    validate_genetic_map(genetic_map)
    return(genetic_map)
}

match_genetic_map_to_L <- function(genetic_map, L, expRate = 0.5, minRate = 0.1, maxRate = 100) {

    ## create joint set of variants in genetic map
    ## then fill in
    ## then map over
    joint_L <- sort(union(L, genetic_map[, "position"]))
    new_genetic_map <- array(NA, c(length(joint_L), 3))
    colnames(new_genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")    
    new_genetic_map[, "position"] <- joint_L
    new_genetic_map[match(genetic_map[, "position"], new_genetic_map[, "position"]), "Genetic_Map.cM."] <-
        genetic_map[, "Genetic_Map.cM."]
    
    ## now - fill in top / bottom / remaining part of new genetic map
    cM <- new_genetic_map[, "Genetic_Map.cM."]
    miss <- which(is.na(cM))
    ##
    ## do the starts and the ends
    ##
    if (1 %in% miss) {
        ## first the first one that needs fixing, then proceed
        first_one_to_fix <- which.max(!is.na(cM)) - 1
        for(iSNP in first_one_to_fix:1) {
            ##
            distance <- joint_L[iSNP + 1] - joint_L[iSNP]
            cM[iSNP] <- cM[iSNP + 1] - (expRate * distance / 1e6)
        }
    }
    if (length(cM) %in% miss) {
        ## again, find the first one that needs fixing, then proceed
        first_one_to_fix <- (length(cM):1)[which.max(!is.na(cM[length(cM):1])) - 1]
        for(iSNP in first_one_to_fix:length(cM)) {
            ##
            distance <- joint_L[iSNP] - joint_L[iSNP - 1]
            cM[iSNP] <- cM[iSNP - 1] + (expRate * distance / 1e6)
        }
    }
    ##
    ## now - fill in the ones in the middle
    ##
    if (sum(is.na(cM)) > 0) {
        
        ## fill in runs of NA here
        r <- rle(is.na(cM))
        ## get the runs that are true - i.e. there are NAs
        starts <- cumsum(c(1, r[["lengths"]]))
        ends <- starts[-1] - 1
        starts <- starts[-length(starts)]
        ##
        starts <- starts[r[["values"]]]
        ends <- ends[r[["values"]]]
        ##         
        for(iRegion in 1:length(starts)) {
            ## these are the undefined SNPs
            start <- starts[iRegion]
            end <- ends[iRegion]
            ## so take the value before and after!
            ## and interpolate
            rate <- (cM[end + 1] - cM[start - 1]) / ((joint_L[end + 1] - joint_L[start - 1]) / 1e6)
            for(iSNP in start:end) {
                ## go forward
                distance <- joint_L[iSNP] - joint_L[iSNP - 1]
                cM[iSNP] <- cM[iSNP - 1] + (rate * distance / 1e6)
            }
        }
    }
    ##
    ## finally, get rates from the genetic map
    ##

    ## 
    ## AM HERE - finish this tomorrow
    ## convert cM rate to sigmaCurrent used rate
    ## 
    
}
