## return an integer matrix with nSNPs rows and 2 columns
## with 0 and 1 entries
## with column1 = hap1 and col2 = hap2
phase_one_sample <- function(
    alphaMatCurrent_tc,
    transMatRate_tc_H,
    alphaHat_t,
    betaHat_t,
    gamma_t,
    eMatGrid_t,
    genProbs,
    c,
    phasing_method,
    phasing_n_votes,
    grid,
    eHapsCurrent_tc
) {

## save(
##     alphaMatCurrent_tc,
##     transMatRate_tc_H,
##     alphaHat_t,
##     betaHat_t,
##     gamma_t,
##     eMatGrid_t,
##     genProbs,
##     c,
##     phasing_method,
##     phasing_n_votes,
##     grid,
##     eHapsCurrent_tc,
##     file = "~/temp.RData")
##     stop("WER")
##     load("~/temp.RData")
    
    if (dim(eHapsCurrent_tc)[[3]] > 1) {
        stop("Phasing does not support S > 1 at this time")
    }
    
    K <- sqrt(nrow(alphaHat_t))
        
    if (phasing_method == 3) {

        paths <- phase_sample_paths_method_3(
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            transMatRate_tc_H = transMatRate_tc_H,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t,
            gamma_t = gamma_t,
            eMatGrid_t = eMatGrid_t,
            c = c,
            phasing_n_votes = phasing_n_votes
        )

        ## identify het sites
        is_homalt <- (genProbs[3, ] > genProbs[2, ])  & (genProbs[3, ] > genProbs[1, ])        
        is_het <- (genProbs[2, ] > genProbs[1, ])  & (genProbs[2, ] > genProbs[3, ])

        ## if not hets - done!
        if (sum(is_het) == 0) {

            final_hap1 <- integer(ncol(genProbs))
            final_hap2 <- integer(ncol(genProbs))
            ##
            final_hap1[is_homalt] <- 1L
            final_hap2[is_homalt] <- 1L
            
            ## split into per-block
            return(cbind(final_hap1, final_hap2))
            
        }

        ## look at paths 
        m <- t(paths[, grid[is_het] + 1])

        ## want 1 if eHapsCurrent_tc for first hap greater than second
        k1 <- (m - 1) %% K  + 1
        k2 <- ceiling(m / K)

        ## for each path
        ## for het SNPs
        ## code 1 if hap1 > hap2 (i.e. 1|0)
        ## then return diff of this
        ## e.g. 1|0 -> 1|0 is a 0 i.e. no swap
        ## e.g. 0|1 -> 0|1 is same thing
        ## but 0|1 >- 1|0 is a 1
        results <- sapply(1:phasing_n_votes, function(i) {
            ##
            kk <- paths[i, grid[is_het] + 1]
            k1 <- (kk - 1)%% K  + 1
            k2 <- ceiling(kk / K)
            ##
            hap1 <- eHapsCurrent_tc[cbind(k1, which(is_het), 1)]
            hap2 <- eHapsCurrent_tc[cbind(k2, which(is_het), 1)]
            ##
            abs(diff(as.integer(hap1 > hap2)))
        })

        ## majority votes
        majority <- rowSums(results) > (phasing_n_votes / 2)
        overall_hap1 <- cumsum(c(1L, as.integer(majority))) %% 2L
        overall_hap2 <- 1L - overall_hap1

        ## build full thing
        final_hap1 <- integer(ncol(genProbs))
        final_hap2 <- integer(ncol(genProbs))
        ##
        final_hap1[is_het] <- overall_hap1
        final_hap2[is_het] <- overall_hap2
        ##
        final_hap1[is_homalt] <- 1L
        final_hap2[is_homalt] <- 1L

        ## split into per-block
        return(cbind(final_hap1, final_hap2))
        
    } else {
        stop("Bad phasing method selected")
    }

}




## 
phase_sample_paths_method_3 <- function(
    alphaMatCurrent_tc,
    transMatRate_tc_H,
    alphaHat_t,
    betaHat_t,
    gamma_t,
    eMatGrid_t,
    c,
    phasing_n_votes
){
    K <- sqrt(nrow(alphaHat_t))
    nGrids <- ncol(alphaHat_t)
    paths <- array(0, c(phasing_n_votes, nGrids))
    for(i in 1:phasing_n_votes){
        paths[i,1] <- sample(1:K**2, prob = gamma_t[, 1], size = 1)
        t2 <- 2
        newprobsum <- numeric(0)
        newprobsum[1] <- 1
        for(t2 in 2:nGrids) {
            U <- runif(1)
            t1 <- t2 - 1
            subtrans <- transMatRate_tc_H[, t1, 1]
            amc <- alphaMatCurrent_tc[, t1, ]
            K1 <- paths[i,t1]
            kt1_1 <- (K1 - 1) %% K  + 1
            kt1_2 <- ceiling(K1 / K)
            transprobstay <- (subtrans[2] * amc[kt1_1] + subtrans[1])*(subtrans[2] * amc[kt1_2] + subtrans[1])
            pstay <- transprobstay*eMatGrid_t[K1,t2]*betaHat_t[K1, t2]*c[,t1]/betaHat_t[K1,t1]
            if(U<pstay){
                paths[i,t2] <- K1
                newprobsum[t2] <- 1
            }
            else{
                newprobs <- numeric(K*K)
                for(kt2_1 in 1:K) {
                    for(kt2_2 in 1:K) {        
                        K2 <- kt2_1 + K * (kt2_2 - 1)
                        jump1 <- subtrans[2] * amc[kt2_1]
                        if (kt1_1 == kt2_1) {
                            jump1 <- jump1 + subtrans[1]
                        }
                        jump2 <- subtrans[2] * amc[kt2_2]
                        if (kt1_2 == kt2_2){
                            jump2 <- jump2 + subtrans[1]
                        }
                        trans_value <- jump1 * jump2
                        newprobs[K2] <- trans_value*eMatGrid_t[K2, t2]*betaHat_t[K2, t2]*c[,t1]/betaHat_t[K1, t1]
                    }
                }
                newprobsum[t2] <- sum(newprobs)
                newprobs <- newprobs[-K1]
                newprobs <- newprobs / sum(newprobs)
                paths[i,t2] <- sample((1:K**2)[-K1], prob = newprobs, size = 1)
            }
        }
        k1old <- (paths[i,] - 1) %% K + 1
        k2old <- ceiling(paths[i,]/K)
        alternate <- k2old + K * (k1old - 1)
        dist <- sum(paths[1,]!=paths[i,])
        altdist <- sum(paths[1,]!=alternate)
        ##if(altdist < dist){
        ##    paths[i,] <- alternate
        ##}
    }
    paths
}


##
vote_at_het_sites <- function(
    paths

) {
    


}

