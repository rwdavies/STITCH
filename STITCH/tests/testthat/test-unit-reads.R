## base qualities in character form, from 25-35 inclusive
## bqs <- paste0(sapply(25 + 33 + 0:10, function(x) rawToChar(as.raw(x))), collapse = "")


test_that("C++ APIs give same sample name for BAM files", {
    sample_name <- "jimmy"
    bam_name <- make_simple_bam(
        file_stem = file.path(tempdir(), sample_name),
        sam = make_simple_sam_text(sample_name = sample_name)
    )
    samtools_name <- get_sample_name_from_bam_file_using_external_samtools(bam_name)
    seqLib_name <- get_sample_name_from_bam_file_using_SeqLib(bam_name)
    expect_equal(samtools_name, seqLib_name)
})


test_that("C++ APIs give same sample name for CRAM files", {
    sample_name <- "jimmy"
    cram_name <- make_simple_cram(
        file_stem = file.path(tempdir(), sample_name),
        sam = make_simple_sam_text(sample_name = sample_name)
    )$cram_file
    samtools_name <- get_sample_name_from_bam_file_using_external_samtools(cram_name)
    seqLib_name <- get_sample_name_from_bam_file_using_SeqLib(cram_name)
    expect_equal(samtools_name, seqLib_name)
})



test_that("an error is thrown when supplied BAM file does not exist", {

    file_that_does_not_exist <- tempfile()

    expect_error(
        get_sample_names_from_bam_or_cram_files(
            file_that_does_not_exist,
            nCores = 1,
            file_type = "BAM",
            verbose = FALSE
        ),
        paste0("Cannot find BAM file:", file_that_does_not_exist)
    )

})



test_that("an error is thrown when supplied BAM file does not have @RG tag", {

    sample_name <- "jimmy"
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), sample_name),
        sam = make_simple_sam_text(sample_name = sample_name, include_rg_tag = FALSE)
    )

    expect_error(
        get_sample_names_from_bam_or_cram_files(
            files = bam_file,
            nCores = 1,
            file_type = "BAM",
            verbose = FALSE
        ),
        paste0("There is no RG tag with sample name in file:", bam_file)
    )

})

test_that("an error is thrown when supplied BAM file has an RG tag but no SM tag", {

    sample_name <- "jimmy"
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), sample_name),
        sam = make_simple_sam_text(
            sample_name = sample_name,
            include_rg_tag = FALSE,
            include_rg_tag_with_no_sm = TRUE
        )
    )

    expect_error(
        get_sample_names_from_bam_or_cram_files(
            files = bam_file,
            nCores = 1,
            file_type = "BAM",
            verbose = FALSE
        ),
        paste0("The RG tags do not contain SM entries for file:", bam_file)
    )

})


test_that("sample names can be properly retrieved from BAM files", {

    sample_names <- c("file1", "file2")
    bam_names <- sapply(sample_names, function(sample_name) {
        bam_name <- make_simple_bam(
            file_stem = file.path(tempdir(), sample_name),
            sam = make_simple_sam_text(sample_name = sample_name)
        )
        return(bam_name)
    })

    bamlist <- tempfile()
    write_names_to_disk(bam_names, bamlist)

    acquired_sample_names <- get_sample_names(
        bamlist = bamlist,
        save = FALSE,
        verbose = FALSE
    )$sampleNames

    expect_equal(
        acquired_sample_names,
        sample_names
    )

})

test_that("sample names can be properly retrieved from CRAM files", {

    sample_names <- c("file1", "file2")
    cram_names <- sapply(sample_names, function(sample_name) {
        cram_name <- make_simple_cram(
            file_stem = file.path(tempdir(), sample_name),
            sam = make_simple_sam_text(sample_name = sample_name)
        )$cram_file
        return(cram_name)
    })

    cramlist <- tempfile()
    write_names_to_disk(cram_names, cramlist)

    acquired_sample_names <- get_sample_names(
        cramlist = cramlist,
        save = FALSE,
        verbose = FALSE
    )$sampleNames

    expect_equal(
        acquired_sample_names,
        sample_names
    )

})



test_that("sample names can be properly retrieved, even if @RG is found somewhere else in the header, like in a program statement", {

    skip("something in newer samtools broke this test")
    
    sample_names <- c("file1", "file2")
    bam_names <- sapply(sample_names, function(sample_name) {
        bam_name <- make_simple_bam(
            file_stem = file.path(tempdir(), sample_name),
            sam = make_simple_sam_text(
                sample_name = sample_name,
                entries = list(
                    c("@PG", "ID:bwa", "PN:bwa", "VN:0.7.12-r1039", "CL:bwa samse ref.fa -r @RG\\tID:7369_8x15\\tSM:", sample_name)
                )
            )
        )
        return(bam_name)
    })

    bamlist <- tempfile()
    write_names_to_disk(bam_names, bamlist)

    acquired_sample_names <- get_sample_names(
        bamlist = bamlist,
        save = FALSE,
        verbose = FALSE
    )$sampleNames

    expect_equal(
        acquired_sample_names,
        sample_names
    )

})






test_that("BAM with one read can be properly interpreted", {

    ## recall - whole ref is A otherwise
    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("G", "G", "G", "G", "G", "G", "G"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf

    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "60",
                  "6M", "*", "0", "0",
                  "GACCGA", ":;<=>?")
            ),
            chr
        )
    )

    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            2, 5,
            matrix(c(-25, 27, -29), ncol = 1),
            matrix(c(4, 5, 6), ncol = 1)
        )
    )


    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-one-read",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))

    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})






test_that("BAM with two reads can be properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "9", "60",
                  "6M", "*", "0", "0",
                  "AACCTT", ":;<=>?"), # 25-30
                c("r002", "0", chr, "11", "60",
                  "3M", "*", "0", "0",
                  "CCT", "@AB") # 31, 32, 33
            ),
            chr
        )
    )
    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            2, 1,
            matrix(c(-25, 27, 29), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        ),
        list(
            1, NA,
            matrix(c(31, 33), ncol = 1),
            matrix(c(1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-two-reads",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    sampleReads[[2]][[2]] <- NA ## blank out: random
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with reads but only mapping quality 0 reads can be properly interpreted", {

    ## recall - whole ref is A otherwise
    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("G", "G", "G", "G", "G", "G", "G"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf

    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "0",
                  "6M", "*", "0", "0",
                  "GACCGA", "::::::")
            ),
            chr
        )
    )

    ## only mq0, so expect empty
    expected_sample_reads <- make_fake_sampleReads()

    sink("/dev/null")
    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-mq-0",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )
    sink()

    load(file_sampleReads(tempdir(), 1, regionName))

    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with two part split read is properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "7", "60",
                  "3M", "*", "0", "0",
                  "AAA", ":;<"), #25, 26, 27
                c("r001", "0", chr, "11", "60",
                  "3M", "*", "0", "0",
                  "CCT", "=>?") # 28, 29, 30
            ),
            chr
        )
    )
    expected_sample_reads <- list(
        list(
            2, 1,
            matrix(c(-27, 28, 30), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-two-part-split",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with two part overlapping split read is properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "9", "60",
                  "6M", "*", "0", "0",
                  "AAAAAA", ":;<=>?"), #25-30
                c("r001", "0", chr, "10", "60",
                  "5M", "*", "0", "0",
                  "ACATA", "@ABCD") # 31-35
            ),
            chr
        )
    )
    expected_sample_reads <- list(
        list(
            4, 1,
            matrix(c(-25, -27, -29, 32, 34), ncol = 1),
            matrix(c(0, 1, 2, 1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-two-part-split-overlapping",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with three part split read is properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "9", "60",
                  "2M", "*", "0", "0",
                  "AA", ":;"), # 25, 26
                c("r001", "0", chr, "10", "60",
                  "2M", "*", "0", "0",
                  "CC", "<="), # 27, 28
                c("r001", "0", chr, "13", "60",
                  "2M", "*", "0", "0",
                  "TT", ">?") # 29, 30
            ),
            chr
        )
    )
    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            2, 1,
            matrix(c(-25, 28, 29), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-three-part",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})



test_that("BAM with three part split read that maps very far apart is removed", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 101),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r000", "0", chr, "9", "60",
                  "2M", "*", "0", "0",
                  "AA", ":;"),
                c("r001", "0", chr, "9", "60",
                  "2M", "*", "0", "0",
                  "AA", ":;"), # 25, 26
                c("r001", "0", chr, "10", "60",
                  "2M", "*", "0", "0",
                  "CC", "<="), # 27, 28
                c("r001", "0", chr, "100", "60",
                  "2M", "*", "0", "0",
                  "TT", ">?"), # 29, 30,
                c("r001", "0", chr, "100", "60",
                  "2M", "*", "0", "0",
                  "AA", "::")
            ),
            chr
        )
    )
    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            0, 0,
            matrix(c(-25), ncol = 1),
            matrix(c(0), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-three-part",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
        iSizeUpperLimit = 20
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})



test_that("BAM with several informative and uninformative reads is properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("A", "A", "A", "A", "A", "A", "A"),
        ALT = c("G", "C", "T", "C", "C", "C", "C")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r000", "0", chr, "5", "60",
                  "2M", "*", "0", "0",
                  "AA", ":;"),
                c("r001", "0", chr, "9", "60",
                  "2M", "*", "0", "0",
                  "AA", ";:"),
                c("r002", "0", chr, "11", "60",
                  "2M", "*", "0", "0",
                  "CC", "AB"),
                c("r001", "0", chr, "11", "60",
                  "2M", "*", "0", "0",
                  "CC", "<="),
                c("r002", "0", chr, "13", "60",
                  "2M", "*", "0", "0",
                  "TT", "BC"),
                c("r003", "0", chr, "21", "60",
                  "2M", "*", "0", "0",
                  "AA", "::")
            ),
            chr
        )
    )
    expected_sample_reads <- list(
        list( ##r001
            1, NA, ## NA as can be random
            matrix(c(-26, 27), ncol = 1),
            matrix(c(0, 1), ncol = 1)
        ),
        list( ## r002
            1, NA, ## NA as is random between 1, 2
            matrix(c(32, 33), ncol = 1),
            matrix(c(1, 2), ncol = 1)
        ),
        list(
            0, 6,
            matrix(-25, ncol = 1),
            matrix(6, ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-several",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
        save_sampleReadsInfo = TRUE
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    sampleReads[[1]][[2]] <- NA ## manually blank out
    sampleReads[[2]][[2]] <- NA ## manually blank out

    expect_equal(
        sampleReads,
        expected_sample_reads
    )

    ## check read length is OK
    load(file_sampleReadsInfo(tempdir(), 1, regionName))

})




test_that("BAM with no reads is properly handled", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c()
            ),
            chr
        )
    )

    expected_sample_reads <- make_fake_sampleReads()

    sink("/dev/null")

    regionName <- "region-name"
    output <- loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-no-reads",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    sink()

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with no informative reads is properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13),
        REF = c("A", "A", "A"),
        ALT = c("G", "C", "T")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "20", "60",
                  "2M", "*", "0", "0",
                  "AA", "::"),
                c("r001", "0", chr, "22", "60",
                  "2M", "*", "0", "0",
                  "CC", "::"),
                c("r001", "0", chr, "24", "60",
                  "2M", "*", "0", "0",
                  "TT", "::")
            ),
            chr
        )
    )

    expected_sample_reads <- make_fake_sampleReads()

    regionName <- "region-name"
    output <- capture.output(loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-no-informative-reads",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    ))

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})


test_that("BAM with one read can use bases in soft clipping", {

    ## recall - whole ref is A otherwise
    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("G", "G", "G", "G", "G", "G", "G"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )

    ## https://samtools.github.io/hts-specs/SAMv1.pdf
    ## recall that soft clipped reads give first non-soft clipped position
    ## so below, 17 is the start of the non-soft clipped part
    ## so the read including soft clipping starts at 15
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "60",
                  "2S2M2S", "*", "0", "0",
                  "ccCCgg", ":;<=>?") # bqs 25-30
            ),
            chr
        )
    )

    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            2, 4,
            matrix(c(25, 27, -29), ncol = 1),
            matrix(c(3, 4, 5), ncol = 1)
        )
    )


    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-1-read-soft-clipped",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
        useSoftClippedBases = TRUE
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(sampleReads, expected_sample_reads)

})



test_that("BAM with one read can remove bases in soft clipping", {

    ## recall - whole ref is A otherwise
    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("G", "G", "G", "G", "G", "G", "G"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )

    ## https://samtools.github.io/hts-specs/SAMv1.pdf
    ## recall that both soft and hard clipped reads have the same number of bases as qualities scores
    ## recall that soft clipped reads give first non-soft clipped position
    ## so below, 14 is the start of the non-soft clipped part
    ## so the read including soft clipping starts at (1-based) 10
    ## include different qualities as well to confirm exact position
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "14", "60",
                  "4S4M3S", "*", "0", "0",
                  "ggccGGGGccc", ":;<=>?@ABCD")
            ),
            chr
        )
    )

    expected_sample_reads <- list(
        list(
            1, NA,
            matrix(c(-30, -32), ncol = 1),
            matrix(c(3, 4), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-BAM-soft",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
        useSoftClippedBases = FALSE
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    sampleReads[[1]][2] <- NA
    expect_equal(sampleReads, expected_sample_reads)

})


test_that("BAM with hard clipped bases are not used", {

    ## recall - whole ref is A otherwise
    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("G", "G", "G", "G", "G", "G", "G"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )

    ## https://samtools.github.io/hts-specs/SAMv1.pdf
    ## so hard clipped bases are not included as sequence or in the pos
    ## so 17 means 17, and C is the first base there
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "60",
                  "5H4M5H", "*", "0", "0",
                  "CCGG", ":;<=")
            ),
            chr
        )
    )

    expected_sample_reads <- list(
        list(
            1, NA,
            matrix(c(25, -27), ncol = 1),
            matrix(c(4, 5), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-BAM-hard",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    sampleReads[[1]][2] <- NA
    expect_equal(sampleReads, expected_sample_reads)

})



test_that("CRAM with one read can be properly interpreted", {

    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("A", "A", "A", "A", "A", "A", "A"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )
    ## https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf

    out <- make_simple_cram(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "60",
                  "6M", "*", "0", "0",
                  "AACCAA", ":;<=>?")
            ),
            chr
        ),
        pos = pos
    )
    cram_file <- out$cram_file
    ref <- out$ref

    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    expected_sample_reads <- list(
        list(
            2, 5,
            matrix(c(-25, 27, -29), ncol = 1),
            matrix(c(4, 5, 6), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        cram_files = cram_file,
        reference = ref,
        N = 1,
        sampleNames = "test-name-0",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))

    expect_equal(
        sampleReads,
        expected_sample_reads
    )

})



test_that("can handle a cigar string of *", {

    chr <- 10
    pos <- cbind(
        CHR = rep(chr, 7),
        POS = c(9, 11, 13, 15, 17, 19, 21),
        REF = c("A", "A", "A", "A", "A", "A", "A"),
        ALT = c("C", "C", "C", "C", "C", "C", "C")
    )
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "17", "60",
                  "6M", "*", "0", "0",
                  "CCCCCC", "::::::"),
                c("r001", "0", chr, "17", "60",
                  "*", "=", "0", "0",
                  "NNNNNN", "######")
            ),
            chr
        )
    )
    expected_sample_reads <- list(
        list(
            2, 5,
            matrix(c(25, 25, 25), ncol = 1),
            matrix(c(4, 5, 6), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        nSNPs = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-0",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    expect_equal(
        sampleReads,
        expected_sample_reads
    )


})


test_that("can properly save qname in light of iSizeUpperLimit", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13, 15, 17, 19),
        REF = rep("A", 6),
        ALT = rep("C", 6)
    )
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r00A", "0", chr, "10", "60",
                  "2M", "*", "0", "0",
                  "AA", "::"),
                c("r001", "0", chr, "10", "60",
                  "2M", "*", "0", "0",
                  "AA", "::"),
                c("r002", "0", chr, "10", "60",
                  "2M", "*", "0", "0",
                  "CC", "::"),
                c("r002", "0", chr, "14", "60",
                  "2M", "*", "0", "0",
                  "AA", "::"),
                c("r001", "0", chr, "18", "60",
                  "2M", "*", "0", "0",
                  "CC", "::")
            ),
            chr
        )
    )

    regionName <- "region-name"
    for(iSizeUpperLimit in c(100000, 5)) {

        output <- loadBamAndConvert(
            iBam = 1,
            L = as.integer(pos[, 2]),
            pos = pos,
            nSNPs = as.integer(nrow(pos)),
            bam_files = bam_file,
            N = 1,
            sampleNames = "test-name-no-informative-reads",
            inputdir = tempdir(),
            regionName = regionName,
            tempdir = tempdir(),
            chr = chr,
            chrStart = 1,
            chrEnd = 100,
            save_sampleReadsInfo = TRUE,
            iSizeUpperLimit = iSizeUpperLimit
        )

        load(file_sampleReads(tempdir(), 1, regionName))
        load(file_sampleReadsInfo(tempdir(), 1, regionName))
        
        expect_equal(
            sampleReads[sampleReadsInfo[, "qname"] == "r002"][[1]][[1]],
            1 ## 0-based, means 2 reads
        )
        if (iSizeUpperLimit == 100000) {
            expect_equal(
                sampleReads[sampleReadsInfo[, "qname"] == "r001"][[1]][[1]],
                1 ## 0-based, means 2 reads
            )
        }

    }
})
