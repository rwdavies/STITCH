## for working with structure rhb_t
## where rows = K
## columns = sets of 32 SNPs


## where hap is a binary vector
int_contract <- function(hap, check = TRUE) {
    if (check) {
        if (class(hap) != "integer") {
            stop("hap not integer")
        }
    }
    nSNPs <- length(hap)
    nbSNPs <- ceiling(nSNPs / 32)    
    if ((nSNPs %% 32) != 0) {
        ## need to pad out final one
        bs <- floor(nSNPs / 32)        
        x <- integer(32)
        x[1:(nSNPs - 32 * bs)] <- tail(hap, nSNPs - 32 * bs)
        hapc <- c(
            packBits(head(hap, 32 * bs), type = "integer"),
            packBits(x, type = "integer")
        )
    } else {
        hapc <- packBits(hap, type = "integer")
    }
    return(hapc)
}

int_contract_manual <- function(hap, check = TRUE) {
    if (length(hap) != 32) {
        stop("manual check just for length 32")
    }
    if (sum(c(rep(0L, 31), 1L) == hap) == 32) {
        return(as.integer(NA))
    }
    a <- 2 ** (0:30)
    if (hap[32] == 1L ) {
        -1 - sum((1L - hap[-32]) * a)
    } else {
        sum((hap[-32]) * a)
    }
}



#' @export
int_expand <- function(hapc, nSNPs = NULL) {
    nbSNPs <- length(hapc)    
    if (is.null(nSNPs)) {
        nSNPs <- nbSNPs * 32
    }
    hap <- integer(nSNPs)
    for(bs in 0:(nbSNPs - 1)) {
        if (bs < (nbSNPs - 1)) {
            hap[32 * bs + 1:32] <- as.integer(intToBits(hapc[bs + 1]))
        } else {
            hap[32 * bs + 1:(nSNPs - 32 * bs)] <- as.integer(intToBits(hapc[bs + 1]))[1:(nSNPs - 32 * bs)]
        }
    }
    return(hap)
}


#' @export
int_determine_rspo <- function(scrambled_vals) {
    m <- sapply(scrambled_vals, function(x) {
        sum(int_expand(x) * 2 ** (1:32))
    })
    scrambled_vals[order(m)]
}

int_expand_manual <- function(hapc, nSNPs = NULL) {
    nbSNPs <- length(hapc)    
    if (is.null(nSNPs)) {
        nSNPs <- nbSNPs * 32
    }
    hap <- integer(nSNPs)
    for(bs in 0:(nbSNPs - 1)) {
        if (bs < (nbSNPs - 1)) {
            hap[32 * bs + 1:32] <- as.integer(intToBits(hapc[bs + 1]))
        } else {
            hap[32 * bs + 1:(nSNPs - 32 * bs)] <- as.integer(intToBits(hapc[bs + 1]))[1:(nSNPs - 32 * bs)]
        }
    }
    return(hap)
}

#' @export
make_rhb_t_from_rhi_t <- function(rhi_t) {
    K <- nrow(rhi_t)
    nSNPs <- ncol(rhi_t)
    nbSNPs <- ceiling(nSNPs / 32)
    rhb_t <- array(as.integer(0), c(K, nbSNPs))
    for(bs in 0:(nbSNPs - 1)) {
        w <- seq(32 * bs + 1, min(nSNPs, 32 * bs + 32))        
        if (bs < (nbSNPs - 1)) {
            ## do not worry about padding
            for(k in 1:K) {
                rhb_t[k, bs + 1] <- rcpp_int_contract(rhi_t[k, w])
            }
        } else {
            x <- integer(32)
            w2 <- 1:length(w)
            for(k in 1:K) {
                x[w2] <- rhi_t[k, w]
                rhb_t[k, bs + 1] <- rcpp_int_contract(x)
            }
        }
    }
    return(rhb_t)
}

#' @export
make_rhb_t_equality <- function(
    rhb_t,
    nSNPs,
    ref_error,
    nMaxDH = NA,
    verbose = TRUE,
    use_hapMatcherR = FALSE,
    zilong = FALSE
) {
    zilong = FALSE
    ## this overrides everything else
    if (is.na(nMaxDH)) {
        if (!use_hapMatcherR) {
            nMaxDH_default <- 2 ** 10 - 1
        } else {
            nMaxDH_default <- 2 ** 8 - 1
        }
        infer_nMaxDH <- TRUE
    } else {
        nMaxDH_default <- nMaxDH
        infer_nMaxDH <- FALSE
    }
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    if (!use_hapMatcherR) {
        ## --- hapMatcher
        ## matrix K x nGrids
        ## 0 = no match
        ## i is match to ith haplotype in distinctHaps i.e. i
        hapMatcher <- array(0L, c(K, nGrids))
        hapMatcherR <- array(as.raw(0), c(K, 1))
    } else {
        ## --- hapMatcherR
        ## same as above, but raw, and therefore 0 through 255, so use 255
        hapMatcher <- array(0L, c(K, 1))
        hapMatcherR <- array(as.raw(0), c(K, nGrids))
    }
    if (infer_nMaxDH) {
        temp_counter <- array(0L, c(nMaxDH_default, nGrids))
    }
    ## --- distinctHapsB
    ## matrix with nMaxDH x nGrids
    ## matrix with the distinct haplotypes
    distinctHapsB <- array(as.integer(0), c(nMaxDH_default, nGrids)) ## store encoded binary
    ## --- all_symbols
    ## list with nGrid entries
    ## each entry is a matrix with each row containing the ID of a symbol, and the 1-based number of entries
    all_symbols <- list(1:nGrids)
    for(iGrid in 1:nGrids) {
        ## can safely ignore the end, it will be zeros what is not captured
        a <- table(rhb_t[, iGrid], useNA = "always")
        if (zilong) {
            a <- a[order(-a)]
        } else {
            ## from stitch
            a <- a[match(int_determine_rspo(names(a)), names(a))]
        }
        a <- a[a > 0]
        ## flip order to get conventional binary order
        if (infer_nMaxDH) {
            if (length(a) > nMaxDH_default) {
                temp_counter[, iGrid] <- a[1:nMaxDH_default]
            } else {
                temp_counter[1:length(a), iGrid] <- a
            }
        }
        names_a <- as.integer(names(a))
        w <- names_a[1:min(length(names_a), nMaxDH_default)]
        distinctHapsB[1:length(w), iGrid] <- w
        ## match against
        if (!use_hapMatcherR) {
            hapMatcher[, iGrid] <- as.integer(match(rhb_t[, iGrid], distinctHapsB[, iGrid]))
            hapMatcher[which(is.na(hapMatcher[, iGrid])), iGrid] <- 0L
        } else {
            m <- match(rhb_t[, iGrid], distinctHapsB[, iGrid])
            hapMatcherR[is.na(m), iGrid] <- as.raw(0)
            hapMatcherR[!is.na(m), iGrid] <- as.raw(m[!is.na(m)])
        }
        ##
        a <- cbind(as.integer(names_a), as.integer(a))
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")
        all_symbols[[iGrid]] <- a
    }
    ##
    ## now, if we're inferring this, choose appropriate re-value downwards
    ##
    if (infer_nMaxDH) {
        running_count <- cumsum(as.numeric(rowSums(temp_counter)) / (as.numeric(nrow(rhb_t)) * as.numeric(nGrids)))
        ## really need to tune this better
        ## basically, larger K, important to set large
        if (K > 50000) {
            thresh <- 0.9999
        } else if (K > 10000) {
            thresh <- 0.9995
        } else if (K > 1000) {
            thresh <- 0.999
        } else {
            thresh <- 0.99
        }
        ## really want to almost never need this, within reason, for large K
        if (sum(running_count > thresh) == 0) {
            suggested_value <- length(running_count)
        } else {
            suggested_value <- which.max(running_count > thresh)
        }
        nMaxDH <- min(
            max(c(2 ** 4 - 1, suggested_value)),
            nMaxDH_default
        )
        if (use_hapMatcherR) {
            if (nMaxDH > 255) {
                nMaxDH <- 255
            }
        }
        if (verbose) {
            print_message(paste0("Using nMaxDH = ", nMaxDH))
        }
        distinctHapsB <- distinctHapsB[1:nMaxDH, , drop = FALSE]
        if (!use_hapMatcherR) {
            hapMatcher[hapMatcher > (nMaxDH)] <- 0L
        } else {
            hapMatcherR[hapMatcherR > (nMaxDH)] <- as.raw(0)
        }
    }
    ##
    ## inflate them too, they're pretty small
    ##
    distinctHapsIE <- array(0L, c(nMaxDH, nSNPs)) ## inflated, with ref_error
    for(iGrid in 0:(nGrids - 1)) {
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        for(k in 1:nMaxDH) {
            distinctHapsIE[k, s:e] <- rcpp_int_expand(distinctHapsB[k, iGrid + 1], nSNPsLocal)
        }
    }
    ##
    distinctHapsIE[distinctHapsIE == 0] <- ref_error
    distinctHapsIE[distinctHapsIE == 1] <- 1 - ref_error
    ##
    ## also, look specifically at the 0 matches
    ##
    if (!use_hapMatcherR) {
        which_hapMatcher_0 <- which(hapMatcher == 0, arr.ind = TRUE) - 1
    } else {
        which_hapMatcher_0 <- which(hapMatcherR == 0, arr.ind = TRUE) - 1
    }
    special_grids <- unique(which_hapMatcher_0[, 2]) + 1 ## this-is-1-based
    eMatDH_special_grid_which <- integer(nGrids)
    eMatDH_special_grid_which[special_grids] <- as.integer(1:length(special_grids))
    if (nrow(which_hapMatcher_0) > 0) {
        ## now build list with them
        x <- which_hapMatcher_0[, 2]
        y <- which((x[-1] - x[-length(x)]) > 0) ## last entry that is OK
        starts <- c(1, y + 1)
        ends <- c(y, length(x))
        ##
        ## eMatDH_special_values
        ##   list of length the number of special grids
        ##   entries are which ones to re-do, and where they are in rhb_t
        ##   entries inside this are 0-based
        eMatDH_special_values_list <- lapply(1:length(starts), function(i) {
            return(as.integer(which_hapMatcher_0[starts[i]:ends[i], 1]))
        })
        ##
        ## eMatDH_special_symbols
        ##   not great name, but these are the actual symbols
        ##
        eMatDH_special_symbols_list <- lapply(1:length(starts), function(i) {
            rhb_t[
                as.integer(which_hapMatcher_0[starts[i]:ends[i], 1]) + 1,
                which_hapMatcher_0[starts[i], 2] + 1
            ]
        })
        ##
        ## make a new matrix version that doesn't need to be converted (ARGH!)
        ## and an index into it
        ##
        eMatDH_special_matrix <- cbind(
            unlist(eMatDH_special_values_list),
            unlist(eMatDH_special_symbols_list)
        )
        eMatDH_special_matrix_helper <- array(as.integer(NA), c(length(eMatDH_special_grid_which > 0), 2))
        eMatDH_special_matrix_helper[eMatDH_special_grid_which > 0, ] <- cbind(as.integer(starts),as.integer(ends))
        ##
        ## fix all_symbols
        ##
        for(iGrid in 1:nGrids) {
            a <- all_symbols[[iGrid]]
            if (nrow(a) > nMaxDH) {
                ##
                ## old behaviour - cap this to be the max
                ## all_symbols[[iGrid]] <- a[1:nMaxDH, ]
                ## 
                ## new behaviour, if more than the max, make final entry the original, with number of missing
                ## 
                a_temp <- a[1:nMaxDH, ]
                n_non_missing <- sum(a_temp[, 2])
                n_missing <- K - n_non_missing
                a_temp <- rbind(a_temp, c(a[1, 1], n_missing))
                ## a_temp <- a_temp[order(-a_temp[, 2]), ]
                all_symbols[[iGrid]] <- a_temp
            }
        }
    } else {
        eMatDH_special_values_list <- list()
        eMatDH_special_matrix <- matrix()
        eMatDH_special_matrix_helper <- matrix()
    }
    nrow_which_hapMatcher_0 <- nrow(which_hapMatcher_0) ## for testing
    ##
    return(
        list(
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_matrix = eMatDH_special_matrix,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            nrow_which_hapMatcher_0 = nrow_which_hapMatcher_0,
            all_symbols = all_symbols
        )
    )
}

## where rows = SNPs, cols = K / hap
make_rhb_from_rhi <- function(rhi) {
    nSNPs <- nrow(rhi)
    K <- ncol(rhi)
    nbSNPs <- ceiling(nSNPs / 32)
    rhb <- array(as.integer(0), c(nbSNPs, K))
    for(bs in 0:(nbSNPs - 1)) {
        w <- seq(32 * bs + 1, min(nSNPs, 32 * bs + 32))
        if (bs < (nbSNPs - 1)) {
            ## do not worry about padding
            for(k in 1:K) {
                rhb[bs + 1, k] <- rcpp_int_contract(rhi[w, k])
            }
        } else {
            x <- integer(32)
            w2 <- 1:length(w)
            for(k in 1:K) {
                x[w2] <- rhi[w, k]
                rhb[bs + 1, k] <- rcpp_int_contract(x)
            }
        }
    }
    return(rhb)
}


## lines to get is 1-based (right) what to get from haps file
load_rhb_at_positions_no_NAs <- function(
    reference_haplotype_file,
    lines_to_get,
    colClasses,
    nSNPs,
    n_haps_snps,
    rh_in_L,
    tempdir = tempdir(),
    load_rhb_method = "R"
) {
    ##
    ## first initialization
    ##
    ## want non-NULL values
    ## determine columns to get
    h_row_1 <- read.table(reference_haplotype_file, nrow = 1)
    if (!is.na(colClasses[1])) {
        haps_to_get <- which((colClasses == "integer") & (h_row_1 != "-")) - 1
    } else {
        haps_to_get <- which(h_row_1 != "-") - 1
    }
    haps_to_get <- as.integer(haps_to_get)
    n_haps <- length(haps_to_get)
    nbSNPs <- ceiling(length(lines_to_get) / 32)    
    ## rhb <- array(as.integer(0), c(nbSNPs, n_haps))
    rhb <- matrix(0L, nrow = nbSNPs, ncol = n_haps)
    ## also, same time, ref_alleleCount
    ref_alleleCount <- array(-1, c(nSNPs, 3)) ## nSNPs here is the pos nSNPs
    ## now - need to build this - argh!
    ## inflate each, add?
    ##
    ## initialize variables needed for processing
    ##
    ## recal -- integer initializes at 0. 1 means length
    bs <- integer(1) ## keep 0-based, also, make vector, for proper pass by reference
    ihold <- integer(1) ## keep 0-based
    hold <- array(0L, c(32, n_haps)) ## hold results here until needed
    final_snp_to_get <- tail(lines_to_get, 1) - 1 ## 0-based
    start_snp <- 1
    binary_get_line <- array(FALSE, n_haps_snps)
    binary_get_line[lines_to_get] <- TRUE
    iSNP <- -1
    final_snp_gotten <- logical(1) ## defaults to false as desired
    ##
    ##
    ## validate_haps(haps, lines_to_get)
    ## convert to rhb around here!
    if (load_rhb_method == "Rcpp") {
        ##
        in.con <- gzfile(reference_haplotype_file, "r")
        while(TRUE) {
            chunk <- readLines(in.con, n = 32)
            if (length(chunk) == 0 | final_snp_gotten[1]) {
                break;
            }
            chunk_length <- length(chunk)
            end_snp <- start_snp + chunk_length - 1
            ## cpp
            Rcpp_rhb_reader_chunk_process(
                rhb = rhb,
                hold = hold,
                chunk = chunk,
                chunk_length = chunk_length,
                start_snp = start_snp,
                end_snp = end_snp,
                bs = bs,
                ihold = ihold,
                haps_to_get = haps_to_get,
                final_snp_to_get = final_snp_to_get,
                n_haps = n_haps,
                binary_get_line = binary_get_line,
                ref_alleleCount = ref_alleleCount,
                rh_in_L = rh_in_L,
                final_snp_gotten = final_snp_gotten
            )
            ##
            start_snp <- start_snp + 32
        }
        close(in.con)
        ## 
    } else if (load_rhb_method == "R") {
        in.con <- gzfile(reference_haplotype_file, "r")
        while(TRUE) {
            chunk <- readLines(in.con, n = 32)
            if (length(chunk) == 0 | final_snp_gotten) {
                break;
            }
            chunk_length <- length(chunk)
            end_snp <- start_snp + chunk_length - 1
            ## 
            ## re-cast this bit in cpp
            ##
            iiSNP <- 0
            while(iiSNP < chunk_length) {
                iiSNP <- iiSNP + 1
                iSNP <- start_snp - 1 + iiSNP ## 1-based
                if (binary_get_line[iSNP]) {
                    ## how to do this in cpp
                    x <- suppressWarnings(as.integer(strsplit(chunk[iiSNP], " ")[[1]]))
                    ## 
                    hold[ihold + 1, ] <- x[haps_to_get + 1]
                    w <- rh_in_L[32 * bs + ihold + 1] ## which SNP is this, 1-based
                    ref_alleleCount[w, 1] <- sum(hold[ihold + 1, ])
                    ref_alleleCount[w, 2] <- n_haps
                    ref_alleleCount[w, 3] <- ref_alleleCount[w, 1] / ref_alleleCount[w, 2]
                    ## if full or the final SNP to get, reset bs / ihold
                    if (((iSNP - 1) == final_snp_to_get) | (ihold == 31)) {
                        ## overflow or end. note - with reset, should be fine not resetting this
                        for(k in 1:n_haps) {
                            rhb[bs + 1, k] <- packBits(hold[, k], type = "integer")
                        }
                        hold[] <- 0L ## reset
                        bs <- bs + 1
                        ihold <- 0
                        if ((iSNP - 1) == final_snp_to_get) {
                            final_snp_gotten <- TRUE
                            iiSNP <- 100
                        }
                    } else {
                        ihold <- ihold + 1
                    }
                }
            }
            ##
            start_snp <- start_snp + 32
        }
        close(in.con)
    }
    ref_alleleCount[ref_alleleCount == -1] <- NA
    return(
        list(
            rhb = rhb,
            ref_alleleCount = ref_alleleCount
        )
    )
    ## 
}


