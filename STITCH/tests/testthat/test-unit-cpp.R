test_that("can do pseudoHaploid updates in C++", {

    ## note - with random numbers, below is kind of weird
    ## but just need to confirm equivalency
    srp <- c(0, 7, 9) # 0-based
    nReads <- length(srp)
    K <- 4
    nSNPs <- 10
    S <- 2
    set.seed(1)
    pRgivenH1_m <- array(runif(nReads * S), c(nReads, S))
    pRgivenH2_m <- array(runif(nReads * S), c(nReads, S))
    
    fbsoL <- lapply(1:2, function(iNor) {
        list_of_gamma_t <- lapply(1:S, function(s) { array(runif(K * nSNPs), c(K, nSNPs)) })
        list_of_eMatRead_t <- lapply(1:S, function(s) { array(runif(K * nReads), c(K, nReads)) })
        return(
            list(
                list_of_gamma_t = list_of_gamma_t,
                list_of_eMatRead_t = list_of_eMatRead_t
            )
        )
    })

    ## R version
    outR <- R_pseudoHaploid_update_9(
        pRgivenH1_m = pRgivenH1_m,
        pRgivenH2_m = pRgivenH2_m,
        list_of_eMatRead_t1 = fbsoL[[1]]$list_of_eMatRead_t,
        list_of_eMatRead_t2 = fbsoL[[2]]$list_of_eMatRead_t,
        list_of_gamma_t1 = fbsoL[[1]]$list_of_gamma_t,
        list_of_gamma_t2 = fbsoL[[2]]$list_of_gamma_t,
        K = K,
        srp = srp
    )

    ## in cpp
    outCPP <- pseudoHaploid_update_model_9(
        pRgivenH1_m = pRgivenH1_m,
        pRgivenH2_m = pRgivenH2_m,
        list_of_eMatRead_t1 = fbsoL[[1]]$list_of_eMatRead_t,
        list_of_eMatRead_t2 = fbsoL[[2]]$list_of_eMatRead_t,
        list_of_gamma_t1 = fbsoL[[1]]$list_of_gamma_t,
        list_of_gamma_t2 = fbsoL[[2]]$list_of_gamma_t,
        K = K,
        srp = srp
    )

    expect_equal(sum(abs(outR$pRgivenH1_m - outCPP$pRgivenH1_m)), 0)
    expect_equal(sum(abs(outR$pRgivenH2_m - outCPP$pRgivenH2_m)), 0)    

})


test_that("forwardBackwardDiploid and forwardBackwardHaploid work", {

    for(gridWindowSize in c(NA, 3)) {
        test_package <- make_fb_test_package(
            K = 4,
            nReads = 8,
            nSNPs = 10,
            gridWindowSize = gridWindowSize,
            S = 2
        )
        S <- test_package$S
        sampleReads <- test_package$sampleReads
        nSNPs <- test_package$nSNPs
        nGrids <- test_package$nGrids
        K <- test_package$K
        transMatRate_tc_H <- test_package$transMatRate_tc_H
        transMatRate_tc_D <- test_package$transMatRate_tc_D    
        alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
        eMatGrid_t <- test_package$list_of_eMatGrid_t[[1]]
        priorCurrent_m <- test_package$priorCurrent_m
        eHapsCurrent_tc <- test_package$eHapsCurrent_tc
        grid <- test_package$grid
        
        ## run through R function
        fbsoL <- run_forward_backwards(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_D = transMatRate_tc_D,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            method = "diploid",
            return_gamma = TRUE,
            grid = grid,
            suppressOutput = 1
        )
        gamma_t <- fbsoL[[1]]$list_of_gamma_t[[1]]
        expect_equal(ncol(gamma_t), nGrids)
        expect_equal(min(gamma_t) >= 0, TRUE)
        expect_equal(max(gamma_t) <= 1, TRUE)

        out <- get_default_hapProbs(pseudoHaploidModel = 9, sampleReads = sampleReads, S = S)
        pRgivenH1L_m <- out$pRgivenH1_m
        pRgivenH2L_m <- out$pRgivenH2_m
        
        ## run through R function
        fbsoL <- run_forward_backwards(
            sampleReads = sampleReads,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_H = transMatRate_tc_H,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            pRgivenH1_m = pRgivenH1L_m,
            pRgivenH2_m = pRgivenH2L_m,
            method = "pseudoHaploid",
            suppressOutput = 1,
            return_gamma = TRUE,
            grid = grid
        )
        print(names(fbsoL))
        ## basic checks
        gamma_t <- fbsoL[[1]]$list_of_gamma_t[[1]]
        expect_equal(ncol(gamma_t), nGrids)
        expect_equal(min(gamma_t) >= 0, TRUE)
        expect_equal(max(gamma_t) <= 1, TRUE)

    }
    
})


test_that("can sample one path from forwardBackwardDiploid", {

    set.seed(10)

    speed_test <- FALSE
    if (speed_test) {
        n_snps <- 25000 ## can set higher
        tmpdir <- "./"
        dir.create(tmpdir, showWarnings = FALSE)
        suppressOutput <- as.integer(0)
        n_reads <- round(n_snps / 10)  ## set to vary coverage
        gridWindowSize <- NA  ## 10 SNPs / grid
    } else {
        gridWindowSize <- NA
        n_snps <- 10
        n_reads <- n_snps * 2
        tmpdir <- tempdir()
        suppressOutput <- as.integer(1)
    }

    K <- 20
    phasemaster <- matrix(
        c(rep(0, n_snps), rep(1, n_snps)),
        ncol = K
    )
    file <- file.path(tmpdir, "package.RData")
    if (speed_test & file.exists(file)) {
        print("loading file!")
        load(file)
    } else {
        test_package <- make_fb_test_package(
            K = 4,
            nReads = 8,
            nSNPs = 10,
            gridWindowSize = 3,
            S = 2
        )
        if (speed_test) {
            save(sampleReads, test_package, file = file)
        }
    }

    S <- test_package$S
    sampleReads <- test_package$sampleReads
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    transMatRate_tc_D <- test_package$transMatRate_tc_D    
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    eMatGrid_t <- test_package$list_of_eMatGrid_t[[1]]
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    grid <- test_package$grid

    set.seed(40)

    ##out2 <- forwardBackwardDiploid_old(
    ##    sampleReads = sampleReads,
    ##    nReads = as.integer(length(sampleReads)),
    ##    pi = pi,
    ##    transMatRate = t(transMatRate),
    ##    alphaMat = t(alphaMat),
    ##    eHaps = t(eHaps),
    ##    maxDifferenceBetweenReads = as.double(1000),
    ##    maxEmissionMatrixDifference = as.double(1000),
    ##    Jmax = as.integer(10),
    ##    suppressOutput = suppressOutput,
    ##    return_a_sampled_path = TRUE,
    ##    blocks_for_output = array(NA, c(1, 1)),
    ##    whatToReturn = as.integer(0)
    ## )
    set.seed(50)
    out <- run_forward_backwards(
        sampleReads = sampleReads,
        priorCurrent_m = priorCurrent_m,
        transMatRate_tc_D = transMatRate_tc_D,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        eHapsCurrent_tc = eHapsCurrent_tc,
        method = "diploid",
        return_a_sampled_path = TRUE
    )

    ##expect_equal(out1$alphaHat_t, out2$alphaHat_t)
    ##expect_equal(out1$betaHat_t, out2$betaHat_t)
    ##expect_equal(out1$gammaK_t, out2$gammaK_t)
    ##expect_equal(out1$jUpdate_t, out2$jUpdate_t)
    ##expect_equal(out1$gammaUpdate_t, out2$gammaUpdate_t)
    ##print(mean(out1$jUpdate_t - out2$jUpdate_t))

    ## basically, these should be the same
    ## given the seed and the ridiculous good fit
    ## marginally_most_likely_path <- apply(out$gamma_t, 2, which.max) ## 1-based
    ## sampled_path <- out$sampled_path_diploid_t[3, ]
    ## 0-based,
    ## should be 1, 0 -> 1
    ## or 0, 1 -> 20
    ##     expect_equal(sum(sampled_path != 1 & sampled_path != 20), 0)

    if (speed_test) {
        print("test haploid")
    }
    out1 <- run_forward_backwards(
        sampleReads = sampleReads,
        priorCurrent_m = priorCurrent_m,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        eHapsCurrent_tc = eHapsCurrent_tc,
        method = "diploid-inbred",
        grid = grid
    )
    ## print("OLD")
    ## priorSum2 <- array(0, K)
    ## jUpdate_t2 <- array(0, c(K, nGrids - 1))
    ## gammaUpdate_t2 <- array(0, c(K, nSNPs, 2))
    ## hapSum_t2 <- array(0, c(K, nGrids))
    ## out2 <- forwardBackwardHaploid_TEMP(
    ##     sampleReads = sampleReads,
    ##     nReads = as.integer(length(sampleReads)),
    ##     pi = pi,
    ##     transMatRate = transMatRate_t_H,
    ##     alphaMat = t(alphaMat),
    ##     eHaps = t(eHaps),
    ##     maxDifferenceBetweenReads = as.double(1000),
    ##     maxEmissionMatrixDifference = as.double(1000),
    ##     Jmax = as.integer(10),
    ##     suppressOutput = 0,
    ##     blocks_for_output = array(NA, c(1, 1)),
    ##     return_extra = FALSE,
    ##     update_in_place = TRUE,
    ##     gammaUpdate_t = gammaUpdate_t2,
    ##     jUpdate_t = jUpdate_t2,
    ##     hapSum_t = hapSum_t2,
    ##     priorSum = priorSum2,
    ##     pass_in_alphaBeta = TRUE,
    ##     alphaHat_t = alphaHat_t,
    ##     betaHat_t = betaHat_t,
    ##     model = 1000,
    ##     run_pseudo_haploid = FALSE,
    ##     pRgivenH1 = array(0, 1),
    ##     pRgivenH2 = array(0, 1)
    ## )
    ## expect_equal(out1, out2)
    ## expect_equal(priorSum, priorSum2)
    ## expect_equal(jUpdate_t2, jUpdate_t)
    ## expect_equal(gammaUpdate_t2, gammaUpdate_t)
    ## expect_equal(hapSum_t, hapSum_t2)



})


test_that("can calculate eMatHapSNP and sample a haploid path", {

    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 10,
        gridWindowSize = 3,
        S = 2
    )
    S <- test_package$S
    sampleReads <- test_package$sampleReads
    nSNPs <- test_package$nSNPs
    nGrids <- test_package$nGrids
    K <- test_package$K
    transMatRate_tc_H <- test_package$transMatRate_tc_H
    transMatRate_tc_D <- test_package$transMatRate_tc_D    
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    eMatRead_t <- test_package$list_of_eMatRead_t[[1]]
    priorCurrent_m <- test_package$priorCurrent_m
    eHapsCurrent_tc <- test_package$eHapsCurrent_tc
    grid <- test_package$grid
    
    read_labels <- as.integer(runif(length(sampleReads)) < 0.5)

    ## I think this is not used
    path <- rcpp_sample_path(
        read_labels = read_labels,
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        maxDifferenceBetweenReads = 1000,
        Jmax = 10,
        priorCurrent_m = priorCurrent_m,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        s = 1
    )

    expect_equal(nrow(path), nGrids)


})
