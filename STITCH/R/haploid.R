R_run_forward_haploid <- function(
    alphaHat_t,
    c,
    eMatGrid_t,
    alphaMat_t,    
    transMatRate_t_H,
    T,
    K,
    pi,
    alphaStart = 0,
    run_fb_subset = FALSE
) {
    nSNPs <- T
    ## 
    if (!run_fb_subset) {
        for(k1 in 0:(K - 1)) {
            alphaHat_t[k1 + 1,0 + 1] <- pi[k1 + 1] * eMatGrid_t[k1 + 1, 0 + 1]
        }
    } else {
        for(k1 in 0:(K - 1)) {        
            alphaHat_t[k1 + 1, 1] <- alphaStart[k1 + 1]
        }
    }
    c[1] <- 1 / sum(alphaHat_t[, 1])
    alphaHat_t[, 1] <- alphaHat_t[, 1] * c[1]
    ## 
    ## 
    alphaConst <- 0
    for(t_0_based in 1:(T - 1)) {
        alphaConst <-
            transMatRate_t_H[1 + 1, t_0_based - 1 + 1] *
            sum(alphaHat_t[, t_0_based - 1 + 1])
        ## 
        alphaHat_t[, t_0_based + 1] <-
            eMatGrid_t[, t_0_based + 1] * (
                transMatRate_t_H[0 + 1, t_0_based - 1 + 1] *
                alphaHat_t[, t_0_based - 1 + 1] + 
                alphaConst *
                alphaMat_t[, t_0_based - 1 + 1]
            )
        ## 
        c[t_0_based + 1] <- 1 / sum(alphaHat_t[, t_0_based + 1])
        alphaHat_t[, t_0_based + 1] <- alphaHat_t[, t_0_based + 1] * c[t_0_based + 1]
    }
    return(
        list(
            alphaHat_t = alphaHat_t,
            c = c
        )
    )
}


R_run_backward_haploid <- function(
    betaHat_t,
    c,
    eMatGrid_t,
    alphaMat_t,    
    transMatRate_t_H
) {
    K <- nrow(alphaMat_t)
    T <- ncol(alphaMat_t) + 1
    for(t_0_based in (T - 2):0) {
        t <- t_0_based + 1
        e_times_b <- eMatGrid_t[, t + 1] * betaHat_t[, t + 1]
        x <- transMatRate_t_H[1 + 1, t] * sum(alphaMat_t[, t] * e_times_b)
        betaHat_t[, t] <- c[t] * (x + transMatRate_t_H[0 + 1, t] * e_times_b)
    }
    return(list(betaHat_t = betaHat_t))
}
