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
            expect_equal(rhi_t[k, ], int_expand(int_contract(rhi_t[k, ]), nSNPs))
            expect_equal(rhi_t[k, ], rcpp_int_expand(int_contract(rhi_t[k, ]), nSNPs))
            expect_equal(int_contract(rhi_t[k, ]), rcpp_int_contract(rhi_t[k, ]))
        }

        ## specifically test edge cases
        for(i in 1:2) {
            hap <- integer(nSNPs)
            hap[] <- 0L
            if (i == 1) { hap[1] <- 1L}
            if (i == 2) { hap[length(hap)] <- 1L}
            if (i == 3) { hap[] <- 1L}
            expect_equal(int_contract(hap), rcpp_int_contract(hap))
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
        rhi <- t(rhi_t)
        rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
        rhb <- t(rhb_t)
        hap <- (rhi_t[1, ] + 1.5 * runif(nSNPs)) / 2.5

        haps_to_get <- c(4, 6, 8)
        rhi_t_subset_r <- rhi_t[haps_to_get, ]
        rhi_t_subset_rcpp <- inflate_fhb_t(rhb_t, haps_to_get - 1, nSNPs) 
        expect_equal(rhi_t_subset_r, rhi_t_subset_rcpp)

        ## other way
        rhi_subset_r <- rhi[, haps_to_get]
        rhi_subset_rcpp <- inflate_fhb(rhb, haps_to_get - 1, nSNPs) 
        expect_equal(rhi_subset_r, rhi_subset_rcpp)
        

    }

})



test_that("math for subsetting for start makes sense", {

    nSNPs <- 100

    ## can crank these up!
    N_haps <- 1000
    max_haps_to_build <- 25
    max_haps_to_project <- 250
    
    set.seed(9910)
    h_all <- array(as.integer(runif(nSNPs * N_haps) > 0.5), c(nSNPs, N_haps)) ## 100 SNPs, 600 samples
    keep_samples_build <- sort(sample(1:N_haps, max_haps_to_build)) 
    h <- h_all[, keep_samples_build]
    
    c <- rowMeans(h) ## rows are SNPs
    h2 <- t(h - c) %*% (h - c)
    out <- eigen(h2, symmetric = TRUE)
    v <- out$vectors
    l <- out$values
    d <- diag(l)
    d_inv <- diag(1 / l)
    h3 <- v %*% d %*% t(v) ## ok, this

    ## can recreate
    expect_lt(max(abs(h3 - h2)), 1e-8)
    expect_lt(max(abs((t(h2) %*% v %*% d_inv)[, 1] - v[, 1])), 1e-8)

    ## can re-calculate first few eigenvalues
    expect_lt(max(abs((t(h2[, 1:3]) %*% v[, 1:4] %*% d_inv[1:4, 1:4])[1, ] - v[1, 1:4])), 1e-8)

    ## versus first four values for a sample
    expect_lt(max(abs((t(h2) %*% v[, 1:4] %*% d_inv[1:4, 1:4])[1, ] - v[1, 1:4])), 1e-8)

    ## check on larger group
    keep_samples_project <- sort(c(keep_samples_build, sample(setdiff(1:N_haps, keep_samples_build), max_haps_to_project - max_haps_to_build)))
    h_project <- h_all[, keep_samples_project]
    h2b <- t(h - c) %*% (h_project - c)
    ## check this is OK
    ib <- 1 ## i-build
    ip <- match(keep_samples_build, keep_samples_project)[ib] ## i-project
    expect_equal(h2[, ib], h2b[, ip])
    ## and can build!
    ## good
    expect_lt(max(abs((t(h2b[, ip]) %*% v[, 1:4] %*% d_inv[1:4, 1:4])[ib, ] - v[ib, 1:4])), 1e-8)

    eigen_cols_to_keep <- 10
    vs <- out$vectors[, 1:eigen_cols_to_keep]
    ls <- out$values[1:eigen_cols_to_keep]
    ds <- diag(ls)
    ds_inv <- diag(1 / ls)
    b <- (t(h2b) %*% vs %*% ds_inv)
    for(i in 1:eigen_cols_to_keep) {
        b[, i] <- b[, i] * ls[i]
    }
    ## local_K <- min(K, nrow(unique(b)))    
    ## out2 <- suppressWarnings(kmeans(b, centers = local_K, iter.max = 100, nstart = 10))

    ## can time using code like the below
    ## wrap whole thing in
    ## f <- function(N_haps) {
    
    ## local_K <- 100
    ## t <- system.time({
    ##     out2 <- suppressWarnings(kmeans(b, centers = local_K, iter.max = 100, nstart = 10))
    ## })
    ## return(t["elapsed"])
    ## }

    ## N_haps_list <- c(1000, 2000, 5000, 7500, 10000, 15000, 20000)
    ## out <- array(0, length(N_haps_list))
    ## names(out) <- as.character(N_haps_list)
    ## for(N_haps in N_haps_list) {
    ##     print(N_haps)
    ##     out[as.character(N_haps)] <- f(N_haps)
    ## }

    ## plot(N_haps_list, out)
    
    ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ## plot(v[, 1] + rnorm(n = 60) / 100, v[, 2] + rnorm(n = 60) / 100, col = cbPalette[out2$cluster])
    ## ## but proportional to values (multiply by that value)
    ## ##
    ## ## pdf("~/temp.pdf", height = 8, width = 20)
    ## par(mfrow = c(2, 5))
    ## for(i in 1:10) {
    ##     plot(v[, 1], v[, i], col = col)
    ## }
    ## dev.off()

    
})

test_that("can sample haps to use for eHaps start efficiently with binary format", {

    ## 1 - normal - less than all values, snps less than max_snps, haps less than both
    ## 2 - more SNPs than max_snps
    ## 3 - more haps than max_haps_to_build, but less than max_haps_to_project
    ## 4 - more haps to max_haps_to_project
    ## 5 - test case - massive - not run for tests
    verbose <- FALSE

    for(i_test in 1:4) {

        N_haps <- 200
        nSNPs <- 100
        max_snps <- 500
        max_haps_to_build <- 250
        max_haps_to_project <- 500
        
        if (i_test == 2) {
            max_snps <- round(nSNPs / 2)
        } else if (i_test == 3) {
            max_haps_to_build <- round(N_haps / 2)
        } else if (i_test == 4) {
            max_haps_to_build <- round(N_haps / 3)
            max_haps_to_project <- round(N_haps / 2)            
        } else if (i_test == 5) {
            N_haps <- 50000
            nSNPs <- 2000
            K <- 100
            max_snps <- 500
            max_haps_to_build <- 250
            max_haps_to_project <- 2500
        }

        if (verbose) {
            print(paste0("N_haps = ", N_haps))
            print(paste0("nSNPs = ", nSNPs))
            print(paste0("max_snps = ", max_snps))
            print(paste0("max_haps_to_build = ", max_haps_to_build))
            print(paste0("max_haps_to_project = ", max_haps_to_project))
        }

        K <- round(N_haps / 2)
        nRefSNPs <- nSNPs        
    
        set.seed(9910)
        reference_haps <- array(as.integer(runif(nSNPs * N_haps) > 0.5), c(nSNPs, N_haps))
        ref_alleleCount <- array(NA, c(nSNPs, 3))
        rh_in_L <- 1:nSNPs
        ref_alleleCount[rh_in_L, 1] <- rowSums(reference_haps)
        ref_alleleCount[rh_in_L, 2] <- ncol(reference_haps)
        ref_alleleCount[, 3] <- ref_alleleCount[, 1] / ref_alleleCount[, 2]
        
        rhi <- reference_haps
        rhi_t <- t(rhi)
        rhb_t <- make_rhb_t_from_rhi_t(rhi_t)
        rhb <- t(rhb_t)
        rm(reference_haps, rhi_t, rhi); gc(reset = TRUE); gc(reset = TRUE)
        
        cols_to_use <- sample_haps_to_use(
            rhb = rhb,
            ref_alleleCount = ref_alleleCount,
            N_haps = N_haps,
            nRefSNPs = nRefSNPs,
            K = K,
            max_snps = max_snps,
            max_haps_to_build = max_haps_to_build,
            max_haps_to_project = max_haps_to_project
        )
        expect_true(length(cols_to_use) > 0) ## placeholder
        
    }
    ## future concern - make this more efficient / faster
    ## not sure matters - could be kmeans is the bad thing etc
    ## h_all <- inflate_fhb(rhb, haps_to_get = keep_samples - 1, nSNPs = nSNPs)
    ## if (!is.na(keep_snps[1])) {
    ##    h_all <- h_all[keep_snps, ] ## meh
    ## }
    ##h2b <- t(h) %*% (h_all - c) 
    ## b <- (t(h2b) %*% vs %*% ds_inv) ## this is the potentially slow one!

    
})
