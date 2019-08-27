test_that("profile reference haplotype binary options", {

    ## here, benchmark putative options
    ## below are numbers on my 2018 mac laptop
    ##                   expr   min     lq    mean median     uq    max neval
    ##                 f1_t() 274.0  293.0  299.40  302.0  303.0  325.0     5
    ##                   f1() 232.0  240.0  263.40  244.0  252.0  349.0     5
    ##                 f2_t() 367.0  388.0  412.60  399.0  402.0  507.0     5
    ##       f3_t(rhi_t, hap) 242.0  261.0  274.60  272.0  282.0  316.0     5
    ##           f3(rhi, hap) 147.0  156.0  184.60  183.0  195.0  242.0     5
    ##          f3b(rhb, hap) 980.0 1000.0 1018.00 1030.0 1040.0 1040.0     5
    ##         f3b2(rhb, hap) 111.0  119.0  148.00  153.0  174.0  183.0     5
    ##  rcpp_f3_t(rhi_t, hap)  93.0  107.0  111.80  113.0  120.0  126.0     5
    ##      rcpp_f3(rhi, hap)  32.8   37.2   47.94   46.3   54.4   69.0     5
    ##     rcpp_f3b(rhb, hap)  22.3   22.3   22.62   22.8   22.8   22.9     5
    ##  rcpp_f4_t(rhi_t, hap)  35.8   50.1   62.06   55.5   59.9  109.0     5
    ## rcpp_f4b_t(rhb_t, hap)  17.6   18.0   18.08   18.0   18.2   18.6     5
    ##      rcpp_f4(rhi, hap) 150.0  151.0  157.60  153.0  153.0  181.0     5
    ## 
    ## here "b" means the binary format and "t" means the transposed format
    ## these numbers suggest "rcpp_4b_t" is best
    ## i.e. store in the "transposed" (row = K/hap, col = SNP) binary format
    ## then, operate column wise
    ##     i.e. per binary SNP (akin to real SNP), go down the column
    ##     adding sum of absolute distance to each of the K outputs
    ## 
    
    library("testthat")
    library("microbenchmark")    ## want

    K <- 10000
    nSNPs <- 32 * 50
    bSNPs <- nSNPs / 32

    ## this is _t, with K in the rows
    rhi_t <- array(as.integer(round(runif(K * nSNPs))), c(K, nSNPs))
    ## want ~0.5 correlation - this works
    hap <- (rhi_t[1, ] + 1.5 * runif(nSNPs)) / 2.5
    cor(hap, rhi_t[1, ]) ** 2
    cor(hap, rhi_t[2, ]) ** 2
    rhi <- t(rhi_t)

    ## store bit-wise using integers
    rhb_t <- array(as.integer(0), c(K, nSNPs / 32))
    for(k in 1:K) {
        for(bs in 0:(nSNPs / 32 - 1)) {
            ##
            rhb_t[k, bs + 1] <- int_contract(rhi_t[k, 32 * bs + 1:32])
        }
    }
    rhb <- t(rhb_t)
    ## access!

    ## ## check
    ## expect_equal(rhi[, 1], int_expand(int_contract(rhi[, 1])))
    ## ## rcpp version
    ## o <- integer(nSNPs)
    ## rcpp_int_expand(o, int_contract(rhi[, 1]), bSNPs)
    ## expect_equal(rhi[, 1], o) 
    ## ## check edge cases
    ## x <- integer(nSNPs); x[] <- 0L
    ## expect_equal(x, int_expand(int_contract(x)))
    ## rcpp_int_expand(o, int_contract(x), bSNPs); expect_equal(x, o)
    ## ## 
    ## x <- integer(nSNPs); x[] <- 1L
    ## expect_equal(x, int_expand(int_contract(x)))
    ## rcpp_int_expand(o, int_contract(x), bSNPs); expect_equal(x, o)
    
    ## build 8bit raw version
    ## rhB_t <- array(

    ## apply
    f1_t <- function() {
        ## here, row is k
        apply(rhi_t, 1, function(y) sum(abs(hap - y)))
    }
    f1 <- function() {
        apply(rhi, 2, function(y) sum(abs(hap - y)))
    }
    f2_t <- function() {
        rowSums(abs(sweep(rhi_t, 2, hap)))
    }
    f3_t <- function(rhi_t, hap) {
        out <- integer(K)
        for(k in 1:K) {
            out[k] <- sum(abs(hap - rhi_t[k, ]))
        }
        return(out)
    }
    f3 <- function(rhi, hap) {
        out <- integer(K)
        for(k in 1:K) {
            out[k] <- sum(abs(hap - rhi[, k]))
        }
        return(out)
    }
    f3b <- function(rhb, hap) {
        out <- integer(K)
        for(k in 1:K) {
            out[k] <- sum(abs(hap - int_expand(rhb[, k])))
        }
        return(out)
    }
    f3b2 <- function(rhb, hap) {
        out <- integer(K)
        o <- integer(nSNPs)
        ## blargh - slightly faster
        for(k in 1:K) {
            rcpp_int_expand2(o, rhb, k - 1, bSNPs)
            out[k] <- sum(abs(hap - o))
        }
        return(out)
    }
    ## OK, am here
    ## what if I push sum / abs super close to this?
    ## binary version seems OK-ish. efficient for memory, although ~2X slower to perform comparison
    
    print("microbenchmark")
    a <- microbenchmark(
        f1_t(),
        f1(),
        f2_t(),
        f3_t(rhi_t, hap),
        f3(rhi, hap),
        f3b(rhb, hap),
        f3b2(rhb, hap), 
        rcpp_f3_t(rhi_t, hap),
        rcpp_f3(rhi, hap),
        rcpp_f3b(rhb, hap),        
        rcpp_f4_t(rhi_t, hap),
        rcpp_f4b_t(rhb_t, hap),        
        rcpp_f4(rhi, hap),
        times = 5
    )
    a[["time"]] <- signif(a[["time"]], 3)
    print(a)

    
    print("do checks")
    a <- f1_t()
    expect_equal(a, f1_t())
    expect_equal(a, f1())
    expect_equal(a, f2_t())
    expect_equal(a, f3_t(rhi_t, hap))
    expect_equal(a, f3(rhi, hap))
    expect_equal(a, f3b(rhb, hap))
    expect_equal(a, f3b2(rhb, hap))    
    
    print("do rcpp checks")
    expect_equal(a, rcpp_f3_t(rhi_t, hap))
    expect_equal(a, rcpp_f3(rhi, hap))
    expect_equal(a, rcpp_f3b(rhb, hap))    
    expect_equal(a, as.numeric(rcpp_f4_t(rhi_t, hap)))
    expect_equal(a, as.numeric(rcpp_f4b_t(rhb_t, hap)))    
    expect_equal(a, as.numeric(rcpp_f4(rhi, hap)))
    

    
})


## test_that("can work with efficient reference_haps", {

##     ## here, specifically want to test
##     int_contract <- function(x) {
##         return(packBits(x, type = "integer"))
##     }
##     int_expand <- function(x) {
##         o <- integer(nSNPs)
##         for(bs in 0:(nSNPs / 32 - 1)) {
##             o[32 * bs + 1:32] <- as.integer(intToBits(x[bs + 1]))
##         }
##         return(o)
##     }

##     set.seed(91)
##     nSNPs <- 32 * 3
##     bSNPs <- nSNPs / 32
##     original_o <- as.integer(round(runif(nSNPs)))
##     o <- integer(nSNPs)
##     x <- int_contract(original_o)
##     ## make x interesting
##     print("original_o")
##     f <- function(original_o) {
##         print(original_o[1:32 + 32 * 0])
##         print(original_o[1:32 + 32 * 1])
##         print(original_o[1:32 + 32 * 2])
##     }
##     f(original_o)
##     print("before")
##     f(o)
##     rcpp_int_expand(o, x, bSNPs)
##     print("after")
##     f(o)

## })


