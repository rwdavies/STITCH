test_that("binary/integer storage format for reference haps makes sense", {

    set.seed(91)
    K <- 10
    
    ## try it just at, under, and over 32
    for(nSNPs in c(90, 95, 96, 97)) {
        
        ## 
        nbSNPs <- ceiling(nSNPs / 32)

        ## make some test data
        rhi_t <- array(as.integer(round(runif(K * nSNPs))), c(K, nSNPs))
        rhb_t <- make_rhb_t_from_rhi_t(rhi_t)

        ## check can expand / contract and get the same answer
        for(k in 1:K) {
            ## check
            expect_equal(rhi_t[k, ], int_expand(int_contract(rhi_t[k, ]), nSNPs))
            expect_equal(rhi_t[k, ], rcpp_int_expand(int_contract(rhi_t[k, ]), nSNPs))
        }
        
        ## specifically test edge cases
        for(i in 1:2) {
            hap <- integer(nSNPs)
            hap[] <- 0L
            if (i == 1) { hap[1] <- 1L}
            if (i == 2) { hap[length(hap)] <- 1L}
            if (i == 3) { hap[] <- 1L}
            expect_equal(hap, int_expand(int_contract(hap), nSNPs))
            expect_equal(hap, rcpp_int_expand(int_contract(hap), nSNPs))
        }


    }
    
})



test_that("can quickly calculate distance between binary haps and haplotype dosage vector", {

    ## am here! work on this
    K <- 10
    nSNPs <- 90
    
    ## try it just at, under, and over 32
    for(nSNPs in c(90, 95, 96, 97)) {
        
        ## 
        nbSNPs <- ceiling(nSNPs / 32)

        ## make some test data
        rhi_t <- array(as.integer(round(runif(K * nSNPs))), c(K, nSNPs))
        rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
        hap <- (rhi_t[1, ] + 1.5 * runif(nSNPs)) / 2.5

        ## calculate using R the easy way
        out <- integer(K)
        for(k in 1:K) {
            out[k] <- sum(abs(hap - rhi_t[k, ]))
        }

        ## calculate in cpp
        out_rcpp <- calc_dist_between_rhb_t_and_hap(rhb_t, hap, nSNPs)
        expect_equal(out, as.numeric(out_rcpp))

    }

})


test_that("can quickly subset rhb_t to a new eHaps like container", {

    ## am here! work on this
    K <- 10
    nSNPs <- 90
    
    ## try it just at, under, and over 32
    for(nSNPs in c(90, 95, 96, 97)) {
        
        ## 
        nbSNPs <- ceiling(nSNPs / 32)

        ## make some test data
        rhi_t <- array(as.integer(round(runif(K * nSNPs))), c(K, nSNPs))
        rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
        hap <- (rhi_t[1, ] + 1.5 * runif(nSNPs)) / 2.5

        haps_to_get <- c(4, 6, 8)
        rhi_t_subset_r <- rhi_t[haps_to_get, ]
        rhi_t_subset_rcpp <- inflate_fhb_t(rhb_t, haps_to_get - 1, nSNPs) 
        expect_equal(rhi_t_subset_r, rhi_t_subset_rcpp)

    }

})

