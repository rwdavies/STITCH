test_that("phase estimation can work with no errors", {
    test <- cbind(
        c(0, 1, 0, 0),
        c(1, 0, 1, 1)
    )
    truth <- test
    pse <- calculate_pse(test, truth)
    expect_equal(
        pse,
        0
    )
})

test_that("phase estimation can work with a het wrongly assigned", {
    test <- cbind(
        c(0, 1, 0, 0),
        c(1, 0, 1, 1)
    )
    truth <- cbind(
        c(0, 0, 0, 0),
        c(1, 1, 1, 1)
    )
    pse <- calculate_pse(test, truth)
    expect_equal(
        pse,
        2 / 3
    )
})

test_that("phase estimation can work with a het wrongly assigned", {
    test <- cbind(
        c(0, 0, 0, 0, 0),
        c(1, 1, 1, 1, 1)
    )
    truth <- cbind(
        c(0, 0, 1, 1, 1),
        c(1, 1, 0, 0, 0)
    )
    pse <- calculate_pse(test, truth)
    expect_equal(
        pse,
        1 / 4
    )
})


test_that("phase estimation skips truth hom sites", {
    test <- cbind(
        c(0, 0, 0, 0, 0),
        c(1, 1, 0, 1, 1)
    )
    truth <- cbind(
        c(0, 0, 0, 1, 1),
        c(1, 1, 1, 1, 1)
    )
    pse <- calculate_pse(test, truth, seed = 1) ## is random
    expect_equal(
        pse,
        1 / 2
    )
})


test_that("phase estimation deals with non-integer test data", {
    test <- cbind(
        c(0.1, 0.2, 0.1, 0.1, 0.1),
        c(0.9, 0.8, 0.7, 0.6, 0.9)
    )
    truth <- cbind(
        c(0, 0, 0, 1, 1),
        c(1, 1, 1, 1, 1)
    )
    pse <- calculate_pse(test, truth, seed = 1) ## is random
    expect_equal(
        pse,
        0 / 3
    )
})


test_that("can proportion estimated haplotype counts ", {

    sampleReads <- list(
        list(
            0, 0, matrix(-10, ncol = 1), matrix(0, ncol = 1)
        ),
        list(
            1, 0,
            matrix(c(-20, -20), ncol = 1),
            matrix(c(0, 1), ncol = 1)            
        ),
        list(
            0, 1, matrix(-10, ncol = 1), matrix(1, ncol = 1)
        )
    )
    T <- 2    
    pRgivenH1 <- c(0.1, 0.1, 0.8)
    pRgivenH2 <- c(0.1, 0.9, 0.2)
    test <- estimate_read_proportions(
        sampleReads,
        pRgivenH1,
        pRgivenH2,
        T
    )
    
    truth <- array(0, c(2, 4))
    colnames(truth) <- c("ER1", "EA1", "ER2", "EA2")
    truth[1, ] <- c(
        9 / 10 * 0.5 + 99 / 100 * 0.1,
        1 / 30 * 0.5 + 1  / 300 * 0.1,
        9 / 10 * 0.5 + 99 / 100 * 0.9,
        1 / 30 * 0.5 + 1  / 300 * 0.9
    )
    truth[2, ] <- c(
        99 / 100 * 0.1 + 9 / 10 * 0.8,
        1  / 300 * 0.1 + 1 / 30 * 0.8,
        99 / 100 * 0.9 + 9 / 10 * 0.2,
        1  / 300 * 0.9 + 1 / 30 * 0.2
    )

    expect_equal(test, truth)

})


test_that("can write out VCF column with expected read counts", {

    T <- 2
    gp <- array(0, c(T, 3))
    gp[1, ] <- c(0, 0.1, 0.9)
    gp[2, ] <- c(0, 0.2, 0.8)
    read_proportions <- array(0, c(T, 4))
    colnames(read_proportions) <- c("ER1", "EA1", "ER2", "EA2")
    read_proportions[1, ] <- c(0.1, 0.2, 0.3, 0.4)
    read_proportions[2, ] <- c(10, 30, 40, 20)
    truth <- c(
        "1/1:0.000,0.100,0.900:1.900:0.100,0.200,0.300,0.400",
        "./.:0.000,0.200,0.800:1.800:10.000,30.000,40.000,20.000"
    )
    gp_t <- t(gp)
    rm(gp)
    
    test1 <- make_column_of_vcf(gp_t, read_proportions)
    test2 <- rcpp_make_column_of_vcf(
        gp_t = gp_t,
        use_read_proportions = TRUE,
        use_state_probabilities = FALSE,
        read_proportions = read_proportions,
        q_t = matrix(),
        add_x_2_cols = FALSE,
        x_t = matrix()
    )    
    expect_equal(test1, truth)
    expect_equal(test2, truth)    

})


test_that("can do method 3 phasing in R and Rcpp", {

    ## library("testthat"); library("STITCH"); library("rrbgen")
    ## ##    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    ## dir <- "~/proj/STITCH/"
    ## setwd(paste0(dir, "/STITCH/R"))
    ## a <- dir(pattern = "*R")
    ## b <- grep("~", a)
    ## if (length(b) > 0) {
    ##     a <- a[-b]
    ## }
    ## o <- sapply(a, source)
    ## setwd(dir)
    ## Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))
    
    ##
    test_package <- make_fb_test_package(
        K = 4,
        nReads = 8,
        nSNPs = 10,
        gridWindowSize = 3,
        S = 1,
        run_fb = TRUE
    )
    
    alphaMatCurrent_tc <- test_package$alphaMatCurrent_tc
    transMatRate_tc_H <- test_package$transMatRate_tc_H    
    alphaHat_t <- test_package$alphaHat_t
    betaHat_t <- test_package$betaHat_t
    gamma_t <- test_package$list_of_gamma_t$gamma_t
    eMatGrid_t <- test_package$eMatGrid_t
    c <- test_package$c
    phasing_n_votes <- 100
    nGrids <- ncol(eMatGrid_t)
    unifs_tc <- array(runif(phasing_n_votes * nGrids * 2), c(phasing_n_votes, nGrids, 2))
    
    out1 <- phase_sample_paths_method_3(    
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaHat_t = alphaHat_t,
        betaHat_t = betaHat_t,
        gamma_t = gamma_t,
        eMatGrid_t = eMatGrid_t,
        c = c,
        phasing_n_votes = phasing_n_votes,
        unifs_tc = unifs_tc
    )

    out2 <- rcpp_phase_sample_paths_method_3(    
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        alphaHat_t = alphaHat_t,
        betaHat_t = betaHat_t,
        gamma_t = gamma_t,
        eMatGrid_t = eMatGrid_t,
        c = c,
        phasing_n_votes = phasing_n_votes,
        unifs_tc = unifs_tc
    )

    expect_equal(out1, out2)

})
