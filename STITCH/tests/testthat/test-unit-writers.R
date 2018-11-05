test_that("can validate B_bit_prob", {
    for(B_bit_prob in c(8, 16, 24, 32)) {
        expect_null(validate_B_bit_prob(B_bit_prob, "bgen"))
    }
    expect_null(validate_B_bit_prob(8, "bgvcf"))
    expect_null(validate_B_bit_prob(9, "bgvcf"))
    expect_error(validate_B_bit_prob(9, "bgen"))
    expect_error(validate_B_bit_prob(-1, "bgen"))
    expect_error(validate_B_bit_prob("-1", "bgen"))
    expect_error(validate_B_bit_prob(FALSE, "bgen"))            
})

test_that("can write column of VCF output", {

    gp <- array(1 / 3, c(4, 3))
    gp[2, ] <- c(0.95, 0, 0.05)
    gp[3, ] <- c(0, 0.95, 0.05)    
    gp_t <- t(round(gp, 3))
    rm(gp)
    out <- make_column_of_vcf(gp_t)
    
    expect_equal(out[1], "./.:0.333,0.333,0.333:0.999")
    expect_equal(out[2], "0/0:0.950,0.000,0.050:0.100")

})

test_that("can write column of VCF output in C++", {

    set.seed(9449)
    nSNPs <- 100
    gp <- array(rbeta(nSNPs * 3, 0.2, 0.2), c(nSNPs, 3))
    gp_t <- t(gp / rowSums(gp))
    rm(gp)    
    out1 <- make_column_of_vcf(gp_t)
    out2 <- rcpp_make_column_of_vcf(
        gp_t = gp_t,
        use_read_proportions = FALSE,
        use_state_probabilities = FALSE,
        read_proportions = matrix(),
        q_t = matrix()
    )    
    expect_equal(out1, out2)

})

test_that("can write column of VCF output in C++ with posterior state probabilities", {

    set.seed(9449)
    nSNPs <- 100
    K <- 5
    gp <- array(rbeta(nSNPs * 3, 0.2, 0.2), c(nSNPs, 3))
    gp_t <- t(gp / rowSums(gp))
    q <- array(runif(nSNPs * K), c(nSNPs, K))
    q_t <- t(q / rowSums(q))
    gp_t[, 1] <- c(0.1, 0.2, 0.7) ## sums to 1
    q_t[, 1] <- c(0, 0.2, 0.1, 0.5, 1.2)
    out1 <- make_column_of_vcf(gp_t = gp_t, q_t = q_t)
    expect_equal(out1[1], "./.:0.100,0.200,0.700:1.600:0.000,0.200,0.100,0.500,1.200")
    
    out2 <- rcpp_make_column_of_vcf(
        gp_t = gp_t,
        use_read_proportions = FALSE,
        use_state_probabilities = TRUE,
        read_proportions = matrix(),
        q_t = q_t
    )    
    expect_equal(out1, out2)

})

