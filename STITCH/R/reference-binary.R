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


