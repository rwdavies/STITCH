test_that("boundind pseudoHaploid makes sense", {

    S <- 3
    nReads <- 6
    pRgivenH1_m <- array(0.5, c(nReads, 3))
    pRgivenH2_m <- array(0.5, c(nReads, 3))
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

    ## a value at 0
    pRgivenH1_m[5, 1] <- e * 0.5
    pRgivenH2_m[5, 1] <- 0
    pRgivenH1_m[5, 2] <- 1 - e * 0.5
    pRgivenH2_m[5, 2] <- 1

    ## 
    pRgivenH1_m[6, 1] <- 0.00000000000000000003060408
    pRgivenH2_m[6, 1] <- e
    pRgivenH1_m[6, 2] <- 0
    pRgivenH2_m[6, 2] <- e
    pRgivenH1_m[6, 3] <- 0
    pRgivenH2_m[6, 3] <- e
    
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

    ## now - even if both below, expect new lowest one to also be bounded
    expect_equal(b[5, 1] > 0, TRUE)
    expect_equal(a[5, 2] < 1, TRUE)

    ## no matter what, bound by e * e
    expect_equal(sum(a < (e * e)) + sum(a > (1 - e * e)), 0)
    expect_equal(sum(b < (e * e)) + sum(b > (1 - e * e)), 0)    

})

