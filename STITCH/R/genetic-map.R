#' @export
load_validate_and_match_genetic_map <- function(
    genetic_map_file,
    L,
    expRate
) {
    ## 
    if (genetic_map_file != "") {
        ##
        print_message("Getting and validating genetic map")
        genetic_map <- get_and_validate_genetic_map(genetic_map_file)
        ## cM is the genetic map at sites from L
        print_message("Match genetic map to desired imputation positions")
        cM <- match_genetic_map_to_L(
            genetic_map = genetic_map,
            L = L,
            expRate = expRate
        )
        print_message("Done with getting genetic map")
        ## now - convert to sigmaCurrent_m
    } else {
        cM <- NULL
    }
    return(cM)
}


initialize_sigmaCurrent_m <- function(
    cM_grid,
    nGen,
    nGrids,
    S,
    dl,
    expRate
) {
    if (!is.null(cM_grid)) {
        ## sigmaCurrent is in morgans, hence first 100 division
        sigmaCurrent_m <- array(exp(-nGen * diff(cM_grid) / 100), c(nGrids - 1, S))
    } else {
        sigmaCurrent_m <- array(exp(-nGen * expRate / 100 / 1000000 * dl), c(nGrids - 1, S))
    }
    return(sigmaCurrent_m)
}


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


fill_in_genetic_map_rate_column <- function(genetic_map) {
    genetic_map[, "COMBINED_rate.cM.Mb."] <- 0
    for(iRow in 1:(nrow(genetic_map) - 1)) {
        distance <-
            genetic_map[iRow + 1, "position"] -             
            genetic_map[iRow, "position"]
        difference <-
            genetic_map[iRow + 1, "Genetic_Map.cM."] -
            genetic_map[iRow, "Genetic_Map.cM."]
        genetic_map[iRow, "COMBINED_rate.cM.Mb."] <-
            difference / (distance / 1e6)
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
    x <- (genetic_map[, 3] - check_col_3) > 1e-6
    if (sum(x) > 0) {
        if (verbose) {
            print(genetic_map[c(max(1, x - 1):min(x + 1, nrow(genetic_map))) , ])
        }
        stop(paste0("Error interpreting genetic map around row ", which.max(x), ". Provided genetic map does not satisfy expectations. See above print. Let g be the genetic map with 3 columns. In 1-based coordinates, recall that g[iRow, 3] must equal g[iRow - 1, 3] + ( (g[iRow, 1] - g[iRow - 1, 1]) * g[iRow - 1, 2])"))
    }
    return(NULL)
}

#' @export
get_and_validate_genetic_map <- function(genetic_map_file) {
    ## 
    genetic_map <- read.table(genetic_map_file, header = TRUE)
    validate_genetic_map(genetic_map)
    return(genetic_map)
}

fill_in_start_and_end_of_cM <- function(cM, joint_L, expRate) {
    ##
    ## do the starts and the ends
    ##
    miss <- which(is.na(cM))
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
    return(cM)
}

fill_in_middle_cM <- function(cM, joint_L) {
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
    return(cM)
}


#' @export
match_genetic_map_to_L <- function(genetic_map, L, expRate = 0.5) {
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
    ## 
    cM <- fill_in_start_and_end_of_cM(cM, joint_L, expRate)     
    cM <- fill_in_middle_cM(cM, joint_L)
    specific_cM <- cM[match(L, joint_L)]
    return(specific_cM)
    ## checks
    ##new_genetic_map[, "Genetic_Map.cM."] <- cM
    ##new_genetic_map <- fill_in_genetic_map_rate_column(new_genetic_map)     
}


#' @title Create a genetic map file from VCF assuming uniform recombination rate
#' @param outfile path of output genetic map file
#' @param vcffile a VCF/BCF file with single chromosome ideally
#' @param region chromosome name if the VCF has multiple chromosome
#' @export
make_genetic_map_file_from_vcf <- function(outfile, vcffile, region = "", expRate = 0.5) {
  print_message("take only a single chromsome from VCF")
  res <- get_pos_from_reference_vcf(vcffile, region)
  if(res[["n_skipped"]] > 0)
    print_message(paste("there are",res[["n_skipped"]],"sites that are ignored, because they are non-SNPs or there is missing genotype or genotypes are not phased in VCF" ))
  dat <- make_genetic_map_file(res[["pos"]], length(res[["pos"]] ), expRate)
  write.table(dat, file = gzfile(outfile), sep = ' ',  quote = FALSE,  row.names = FALSE )
}



