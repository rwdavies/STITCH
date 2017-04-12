test_that("downsample can remove a read to ensure downsampleToCov is respected", {

    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    sampleReads <- list(
        list(1, 0, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1)),
        list(1, 0, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1))        
    )
    sampleReadsInfo <- data.frame(
        qname = c("read1", "read2"),
        strand = c("+", "-")
    )

    set.seed(1)
    out <- downsample(
        sampleReads,
        iBam = 1,
        downsampleToCov = 1,
        sampleNames = "sample",
        sampleReadsInfo = sampleReadsInfo,
        verbose = FALSE
    )

    expect_equal(
        out$sampleReads,
        sampleReads[2]
    )
    expect_equal(out$sampleReadsInfo, sampleReadsInfo[2,] )

})


test_that("downsample can remove a read to ensure downsampleToCov is respected at all SNPs, not just lead SNPs", {

    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    sampleReads <- list(
        list(1, 0, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1)),
        list(1, 1, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1))        
    )
    sampleReadsInfo <- data.frame(
        qname = c("read1", "read2"),
        strand = c("+", "-")
    )

    set.seed(1)
    out <- downsample(
        sampleReads,
        iBam = 1,
        downsampleToCov = 1,
        sampleNames = "sample",
        sampleReadsInfo = sampleReadsInfo,
        verbose = FALSE
    )

    expect_equal(
        out$sampleReads,
        sampleReads[2]
    )
    expect_equal(out$sampleReadsInfo, sampleReadsInfo[2,] )

})



test_that("downsampling twice has no effect", {

    sampleReads <- list(
        list(1, 0, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1)),
        list(1, 0, matrix(c(-25, -25), ncol = 1), matrix(c(0, 1), ncol = 1))        
    )
    sampleReadsInfo <- data.frame(
        qname = c("read1", "read2"),
        strand = c("+", "-")
    )

    set.seed(1)
    out1 <- downsample(
        sampleReads,
        iBam = 1,
        downsampleToCov = 1,
        sampleNames = "sample",
        sampleReadsInfo = sampleReadsInfo,
        verbose = FALSE
    )
    out2 <- downsample(
        out1$sampleReads,
        iBam = 1,
        downsampleToCov = 1,
        sampleNames = "sample",
        sampleReadsInfo = out1$sampleReadsInfo,
        verbose = FALSE
    )

    expect_equal(
        out1$sampleReads,
        out2$sampleReads        
    )
    expect_equal(
        out1$sampleReadsInfo,
        out2$sampleReadsInfo    
    )


})



test_that("downsample can remove reads in a more complicated scenario to ensure that downsampleToCov is respected at all SNPs, not just lead SNPs", {

    ## make lots of reads
    ## make sure nothing weird happens
    set.seed(1)
    T <- 20
    nReads <- 100
    sampleReads <- lapply(1:nReads, function(i) {
        x <- sort(sample(1:T, 5))
        y <- sample(10:20, 5) * sample(c(-1, 1), 5, replace = TRUE)
        return(
            list(
                length(x) - 1,
                round(median(x)),
                matrix(y, ncol = 1),
                matrix(x, ncol = 1)
            )
        )
    })
    sampleReadsInfo <- data.frame(
        qname = paste0("read", 1:nReads),
        strand = sample(c("-", "+"), nReads, replace = TRUE)
    )

    set.seed(2)
    out <- downsample(
        sampleReads,
        iBam = 1,
        downsampleToCov = 10,
        sampleNames = "sample",
        sampleReadsInfo = sampleReadsInfo,
        verbose = FALSE
    )
    old_depth_per_SNP <- get_depth_per_SNP_for_sampleReads(sampleReads = sampleReads)
    new_depth_per_SNP <- get_depth_per_SNP_for_sampleReads(sampleReads = out$sampleReads)
    expect_equal(sum(new_depth_per_SNP > 10), 0)

})

test_that("downsample works with fake sampleReads ", {
    
    sampleReads <- make_fake_sampleReads()
    nReads <- length(sampleReads)
    sampleReadsInfo <- data.frame(
        qname = paste0("read", 1:nReads),
        strand = sample(c("-", "+"), nReads, replace = TRUE)
    )
    
    out <- downsample(
        sampleReads,
        iBam = 1,
        downsampleToCov = 10,
        sampleNames = "sample",
        sampleReadsInfo = sampleReadsInfo,
        verbose = FALSE
    )
    expect_equal(out$sampleReads, sampleReads)
    expect_equal(out$sampleReadsInfo, sampleReadsInfo)

})

