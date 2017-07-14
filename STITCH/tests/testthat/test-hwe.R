## am here
## write tests that
## test for some arbitrary numbers c(0, 50, 100), c(100, 50, 0), c(0, 5, 10), c(10, 5, 0)
## as well as edge cases, c(0, 0, 0), c(0, 0, 1), etc

test_that("HWE in C++ is the same as in R", {

    ## have confirmed this fails
    mat <- rbind(
        c(0, 0, 0),
        c(0, 0, 1)
    )
    apply(mat, 1, function(x) {
        expect_equal(
            calculate_hwe_p(x),
            rcpp_calculate_hwe_p(x)
        )
    })

})
