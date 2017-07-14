test_that("HWE in C++ is the same as in R", {

    mat <- rbind(
        c(0, 0, 0),
        c(0, 0, 1),
        c(0, 1, 0),
        c(1, 0, 0),
        c(0, 1, 1),
        c(1, 0, 1),
        c(1, 1, 0),
        c(1, 1, 1),
        c(5, 10, 20),
        c(20, 10, 5),
        c(10000, 100, 10),
        c(10000, 101, 10),        
        c(10000, 102, 10),
        c(10000, 103, 10),        
        c(10000, 10000, 10)        
    )
    apply(mat, 1, function(x) {
        expect_equal(
            calculate_hwe_p(x),
            rcpp_calculate_hwe_p(x)
        )
    })

})
