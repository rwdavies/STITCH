test_that("can validate B_bit_prob", {
    for(B_bit_prob in c(8, 16, 24, 32)) {
        expect_null(validate_B_bit_prob(B_bit_prob, "bgen"))
    }
    expect_null(validate_B_bit_prob(8, "gvcf"))
    expect_null(validate_B_bit_prob(9, "gvcf"))
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
    out2 <- rcpp_make_column_of_vcf(gp_t, 0, matrix())    
    expect_equal(out1, out2)

})

