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
