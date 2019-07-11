test_that("boundind pseudoHaploid makes sense", {

    S <- 3
    nReads <- 4
    pRgivenH1_m <- array(0.5, c(4, 3))
    pRgivenH2_m <- array(0.5, c(4, 3))
    e <- 0.001
    ## each individually
    pRgivenH1_m[1, 1] <- e * 0.5
    pRgivenH2_m[1, 2] <- e * 0.5
    ## both together
    pRgivenH1_m[2, 2] <- e * 0.5
    pRgivenH2_m[2, 2] <- e * 0.5

    ## other way
    pRgivenH1_m[3, 1] <- 1 - e * 0.5
    pRgivenH2_m[3, 2] <- 1 - e * 0.5
    ## both together
    pRgivenH1_m[4, 2] <- 1 - e * 0.5
    pRgivenH2_m[4, 2] <- 1 - e * 0.5
    
    ##     
    out <- bound_pRgivenH(
        pRgivenH1_m = pRgivenH1_m,
        pRgivenH2_m = pRgivenH2_m,
        e = e
    )

    ##
    a <- out$pRgivenH1_m
    b <- out$pRgivenH2_m
    ## 
    expect_equal(a[1, 1], e);       expect_equal(b[1, 1], 0.5);
    expect_equal(a[1, 2], 0.5);     expect_equal(b[1, 2], e);
    ## together
    expect_equal(a[2, 2], e / 2);   expect_equal(b[2, 2], e / 2)
    
    ## other way
    expect_equal(a[3, 1], 1 - e);   expect_equal(b[3, 1], 0.5);
    expect_equal(a[3, 2], 0.5);     expect_equal(b[3, 2], 1 - e);
    ## together
    expect_equal(a[4, 2], 1 - e / 2);   expect_equal(b[4, 2], 1 - e / 2)
    

})

