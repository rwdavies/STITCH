test_that("can do pseudoHaploid updates in C++", {

    ## note - with random numbers, below is kind of weird
    ## but just need to confirm equivalency
    srp <- c(0, 7, 9) # 0-based
    n_reads <- length(srp)
    K <- 4
    T <- 10
    set.seed(1)
    pRgivenH1 <- runif(n_reads)
    pRgivenH2 <- runif(n_reads)
    fbsoL <- list(
        list(eMatHap = matrix(runif(K * n_reads), nrow = n_reads)),
        list(eMatHap = matrix(runif(K * n_reads), nrow = n_reads))
    )
    for(iNor in 1:2) {
        fbsoL[[iNor]]$gamma <- matrix(runif(K * T), nrow = T)
        fbsoL[[iNor]]$eMatHap_t <- t(fbsoL[[iNor]]$eMatHap)
        fbsoL[[iNor]]$gamma_t <- t(fbsoL[[iNor]]$gamma)
    }

    ## model 9 - original configuration
    x <- pRgivenH1/(pRgivenH1+pRgivenH2)
    y1 <- (fbsoL[[1]]$eMatHap - (1-x) * pRgivenH2) / x
    y2 <- (fbsoL[[2]]$eMatHap - x * pRgivenH1) / (1 - x)
    pRgivenH1_new <- rowSums(fbsoL[[1]]$gamma[srp + 1,] * y1)
    pRgivenH2_new <- rowSums(fbsoL[[2]]$gamma[srp + 1,] * y2)

    out1 <- pseudoHaploid_update_9(
        pRgivenH1 = pRgivenH1,
        pRgivenH2 = pRgivenH2,
        eMatHap_t1 = fbsoL[[1]]$eMatHap_t,
        eMatHap_t2 = fbsoL[[2]]$eMatHap_t,
        gamma_t1 = fbsoL[[1]]$gamma_t,
        gamma_t2 = fbsoL[[2]]$gamma_t,
        K = K,
        srp = srp
    )

    out2 <- pseudoHaploid_update_model_9(
        pRgivenH1 = pRgivenH1,
        pRgivenH2 = pRgivenH2,
        eMatHap_t1 = fbsoL[[1]]$eMatHap_t,
        eMatHap_t2 = fbsoL[[2]]$eMatHap_t,
        gamma_t1 = fbsoL[[1]]$gamma_t,
        gamma_t2 = fbsoL[[2]]$gamma_t,
        K = K,
        srp = srp
    )

    expect_equal(sum(abs(out1$pRgivenH1 - pRgivenH1_new)), 0)
    expect_equal(sum(abs(out1$pRgivenH2 - pRgivenH2_new)), 0)
    expect_equal(sum(abs(out2$pRgivenH1 - pRgivenH1_new)), 0)
    expect_equal(sum(abs(out2$pRgivenH2 - pRgivenH2_new)), 0)


})


test_that("forwardBackwardDiploid and forwardBackwardHaploid work", {

    n_snps <- 10 ## set to 10000 to check times better
    K <- 20

    phasemaster <- matrix(
        c(rep(0, n_snps), rep(1, n_snps)),
        ncol = K
    )
    data_package <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_snps * 2,
        seed = 2,
        chr = 10,
        K = K,
        phasemaster = phasemaster,
        reads_span_n_snps = 3,
        n_cores = 1
    )
    ##tmpdir = "/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/"
    ##save(data_package, file = "/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/test.RData")
    ##load("/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/test.RData")

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = data_package$L,
        pos = data_package$pos,
        nSNPs = data_package$nSNPs,
        bam_files = data_package$bam_files,
        N = 1,
        sampleNames = data_package$sample_names,
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = data_package$chr,
        chrStart = 1,
        chrEnd = max(data_package$pos[, 2]) + 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    eHaps <- array(runif(n_snps * K), c(n_snps, K))
    sigma <- runif(n_snps - 1)
    alphaMat <- array(runif((n_snps - 1) * K), c(n_snps - 1, K))
    x <- sigma
    transMatRate <- cbind(x ** 2, x * (1 - x), (1 - x) ** 2)
    pi <- runif(K) / K


    out <- forwardBackwardDiploid(
        sampleReads = sampleReads,
        nReads = as.integer(length(sampleReads)),
        pi = pi,
        transMatRate = t(transMatRate),
        alphaMat = t(alphaMat),
        eHaps = t(eHaps),
        maxDifferenceBetweenReads = as.double(1000),
        maxEmissionMatrixDifference = as.double(1e4),
        whatToReturn = as.integer(0),
        Jmax = as.integer(10),
        suppressOutput = as.integer(1)
    )



    pRgivenH1L <- runif(length(sampleReads))
    pRgivenH2L <- runif(length(sampleReads))

    out <- forwardBackwardHaploid(
        sampleReads = sampleReads,
        nReads = as.integer(length(sampleReads)),
        Jmax = as.integer(10),
        pi = pi,
        pRgivenH1 = pRgivenH1L,
        pRgivenH2 = pRgivenH2L,
        pState = eHaps,
        eHaps = t(eHaps),
        alphaMat = t(alphaMat),
        transMatRate = t(transMatRate),
        maxDifferenceBetweenReads = as.double(1000),
        maxEmissionMatrixDifference = as.double(1e4),        
        whatToReturn = as.integer(0),
        suppressOutput=as.integer(1),
        model=as.integer(9)
    )


})


test_that("can sample one path from forwardBackwardDiploid", {

    set.seed(10)
    n_snps <- 10 ## set to 10000 to check times better
    K <- 20

    phasemaster <- matrix(
        c(rep(0, n_snps), rep(1, n_snps)),
        ncol = K
    )
    data_package <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_snps * 2,
        seed = 2,
        chr = 10,
        K = K,
        phasemaster = phasemaster,
        reads_span_n_snps = 3,
        n_cores = 1
    )
    ##tmpdir = "/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/"
    ##save(data_package, file = "/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/test.RData")
    ##load("/data/smew1/rdavies/stitch_development/STITCH_v1.2.7_development/cppdir/test.RData")

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = data_package$L,
        pos = data_package$pos,
        nSNPs = data_package$nSNPs,
        bam_files = data_package$bam_files,
        N = 1,
        sampleNames = data_package$sample_names,
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = data_package$chr,
        chrStart = 1,
        chrEnd = max(data_package$pos[, 2]) + 100
    )

    set.seed(40)
    load(file_sampleReads(tempdir(), 1, regionName))
    eHaps <- array(runif(n_snps * K), c(n_snps, K))
    sigma <- rep(0.999, n_snps - 1)
    alphaMat <- array(1 / K / K, c(n_snps - 1, K))
    x <- sigma
    transMatRate <- cbind(x ** 2, x * (1 - x), (1 - x) ** 2)
    pi <- runif(K) / K
    eHaps[, 1] <- 0.01
    eHaps[, 2] <- 0.99

    out <- forwardBackwardDiploid(
        sampleReads = sampleReads,
        nReads = as.integer(length(sampleReads)),
        pi = pi,
        transMatRate = t(transMatRate),
        alphaMat = t(alphaMat),
        eHaps = t(eHaps),
        maxDifferenceBetweenReads = as.double(1000),
        maxEmissionMatrixDifference = as.double(1000),
        whatToReturn = as.integer(0),
        Jmax = as.integer(10),
        suppressOutput = as.integer(1),
        return_a_sampled_path = as.integer(1)
    )

    ## basically, these should be the same
    ## given the seed and the ridiculous good fit
    marginally_most_likely_path <- apply(out$gamma_t, 2, which.max) ## 1-based


    sampled_path <- out$sampled_path_diploid[, 3]
    ## 0-based,
    ## should be 1, 0 -> 1
    ## or 0, 1 -> 20
    expect_equal(sum(sampled_path != 1 & sampled_path != 20), 0)

})


test_that("can calculate eMatHapSNP and sample a path", {

    n_snps <- 10 ## set to 10000 to check times better
    K <- 4
    phasemaster <- matrix(
        c(rep(0, n_snps), rep(1, n_snps)),
        ncol = K
    )
    data_package <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_snps * 2,
        seed = 2,
        chr = 10,
        K = K,
        phasemaster = phasemaster,
        reads_span_n_snps = 3,
        n_cores = 1
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = data_package$L,
        pos = data_package$pos,
        nSNPs = data_package$nSNPs,
        bam_files = data_package$bam_files,
        N = 1,
        sampleNames = data_package$sample_names,
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = data_package$chr,
        chrStart = 1,
        chrEnd = max(data_package$pos[, 2]) + 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    eHaps <- array(runif(n_snps * K), c(n_snps, K))

    eMatHap_t <- rcpp_make_eMatHap_t(
        sampleReads = sampleReads,
        nReads = length(sampleReads),
        eHaps_t = t(eHaps),
        maxDifferenceBetweenReads = 1000,
        Jmax = 10
    )

    sigma <- runif(n_snps - 1)
    alphaMat <- array(runif((n_snps - 1) * K), c(n_snps - 1, K))
    x <- sigma
    transMatRate <- cbind(x, 1 - x)
    pi <- runif(K) / K

    read_labels <- as.integer(runif(length(sampleReads)) < 0.5)

    path <- rcpp_sample_path(
        read_labels = read_labels,
        eMatHap_t = eMatHap_t,
        sampleReads = sampleReads,
        nReads = length(sampleReads),
        eHaps_t = t(eHaps),
        maxDifferenceBetweenReads = 1000,
        Jmax = 10,
        pi = pi,
        transMatRate_t = t(transMatRate),
        alphaMat_t = t(alphaMat)
    )

})


