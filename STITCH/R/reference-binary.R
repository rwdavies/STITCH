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
                rhb_t[k, bs + 1] <- int_contract(rhi_t[k, w])
            }
        } else {
            x <- integer(32)
            w2 <- 1:length(w)
            for(k in 1:K) {
                x[w2] <- rhi_t[k, w]
                rhb_t[k, bs + 1] <- int_contract(x)
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
                rhb[bs + 1, k] <- int_contract(rhi[w, k])
            }
        } else {
            x <- integer(32)
            w2 <- 1:length(w)
            for(k in 1:K) {
                x[w2] <- rhi[w, k]
                rhb[bs + 1, k] <- int_contract(x)
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
    tempdir = tempdir(),
    load_rhb_method = "R"
) {
    ##
    ## first initialization
    ##
    ## want non-NULL values
    ## determine columns to get
    h_row_1 <- read.table(reference_haplotype_file, nrow = 1)
    haps_to_get <- which((colClasses == "integer") & (h_row_1 != "-")) - 1
    n_haps <- length(haps_to_get)
    nSNPs <- length(lines_to_get)
    nbSNPs <- ceiling(nSNPs / 32)    
    rhb <- array(as.integer(0), c(nbSNPs, n_haps))
    ##
    ## initialize variables needed for processing
    ##
    bs <- 0 ## keep 0-based
    ihold <- 0 ## keep 0-based
    hold <- array(0L, c(32, n_haps)) ## hold results here until needed
    final_snp_to_get <- tail(lines_to_get, 1) - 1 ## 0-based
    start_snp <- 1
    binary_get_line <- array(FALSE, n_haps_snps)
    binary_get_line[lines_to_get] <- TRUE
    iSNP <- -1
    ##
    ##
    ## validate_haps(haps, lines_to_get)
    ## convert to rhb around here!
    if (load_rhb_method == "Rcpp") {
        ##
        in.con <- gzfile(reference_haplotype_file, "r")
        while(TRUE) {
            chunk <- readLines(in.con, n = 32)
            if (length(chunk) == 0 | final_snp_to_get < (iSNP - 1)) {
                break;
            }
            chunk_length <- length(chunk)
            end_snp <- start_snp + chunk_length - 1
            ## cpp
            Rcpp_rhb_reader_chunk_process(
                rhb,
                hold,
                chunk,
                chunk_length,
                start_snp,
                end_snp,
                bs,
                ihold,
                haps_to_get,
                final_snp_to_get,
                n_haps,
                binary_get_line
            )
            ##
        }
        ## 
    } else if (load_rhb_method == "R") {
        in.con <- gzfile(reference_haplotype_file, "r")
        while(TRUE) {
            chunk <- readLines(in.con, n = 32)
            if (length(chunk) == 0 | final_snp_to_get < (iSNP - 1)) {
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
    return(rhb)
    ## 
}

