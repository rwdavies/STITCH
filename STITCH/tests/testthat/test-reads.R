test_that("loading windows can be properly calculated for NA chrStart and chrEnd", {
    ## if NA, load entire chromosome?
    out <- determine_loading_windows(
        chrStart = NA,
        chrEnd = NA,
        chrLength = 10,
        width_of_loading_window = 1e6
    )
    expect_equal(out$start, c(1))
    out <- expect_equal(out$end, c(10))
})

test_that("loading windows can be properly calculated for given chrStart and chrEnd", {
    out <- determine_loading_windows(
        chrStart = 1,
        chrEnd = 10,
        chrLength = NA,
        width_of_loading_window = 6
    )
    expect_equal(out$start, c(1, 7))
    expect_equal(out$end, c(6, 10))
})

test_that("loading windows can be properly calculated for given chrStart and chrEnd with exact matching", {
    out <- determine_loading_windows(
        chrStart = 1,
        chrEnd = 7,
        chrLength = NA,
        width_of_loading_window = 7
    )
    expect_equal(out$start, c(1))
    expect_equal(out$end, c(7))
})

test_that("loading windows can be properly calculated for given chrStart and chrEnd with window of 1 base", {
    out <- determine_loading_windows(
        chrStart = 1,
        chrEnd = 8,
        chrLength = NA,
        width_of_loading_window = 7
    )
    expect_equal(out$start, c(1, 8))
    expect_equal(out$end, c(7, 8))
})

test_that("loading windows can be properly calculated for given chrStart and chrEnd with really big window size", {
    out <- determine_loading_windows(
        chrStart = 1,
        chrEnd = 9,
        chrLength = NA,
        width_of_loading_window = 100
    )
    expect_equal(out$start, c(1))
    expect_equal(out$end, c(9))
})


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





test_that("results from different windows are properly merged", {

    sampleReadsAcrossRegions <- list(
        list(
            sampleReadsRaw = list(
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(0, ncol = 1),
                    0
                ),
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(1, ncol = 1),
                    1
                )
            ),
            qname = c("read1", "read2"),
            strand = c("-", "+")
        ),
        list(
            sampleReadsRaw = list(
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(2, ncol = 1),
                    0 ## now changed
                ),
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(5, ncol = 1),
                    1
                )
            ),
            qname = c("read1", "read1"),
            strand = c("-", "-")
        )
    )
    expected_sampleReadsRaw <- list(
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(0, ncol = 1),
            0
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1),
            1
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(2, ncol = 1),
            2
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(5, ncol = 1),
            3
        )
    )

    out <- combineReadsAcrossRegions(sampleReadsAcrossRegions)
    expect_equal(out$sampleReadsRaw, expected_sampleReadsRaw)
    expect_equal(out$strand, c("-", "+", "-", "-"))    
    expect_equal(out$qname, c("read1", "read2", "read1", "read1"))

})




test_that("results from different windows are properly merged with unnecessary qname", {

    sampleReadsAcrossRegions <- list(
        list(
            sampleReadsRaw = list(
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(0, ncol = 1),
                    0
                ),
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(1, ncol = 1),
                    1
                )
            ),
            qname = c("read1", "read2", "read3"),
            strand = c("-", "+", "-")
        ),
        list(
            sampleReadsRaw = list(
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(2, ncol = 1),
                    0 ## now changed
                ),
                list(
                    0, 0,
                    matrix(-10, ncol = 1),
                    matrix(5, ncol = 1),
                    1
                )
            ),
            qname = c("read1", "read1", "read4"),
            strand = c("-", "-", "+")
        )
    )
    expected_sampleReadsRaw <- list(
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(0, ncol = 1),
            0
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1),
            1
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(2, ncol = 1),
            3
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(5, ncol = 1),
            4
        )
    )

    out <- combineReadsAcrossRegions(sampleReadsAcrossRegions)
    expect_equal(out$sampleReadsRaw, expected_sampleReadsRaw)
    expect_equal(out$strand, c("-", "+", "-", "-", "-", "+"))    
    expect_equal(out$qname, c("read1", "read2", "read3", "read1", "read1", "read4"))

})



test_that("reads are appropriately merged from sampleReadsRaw to sampleReads", {

    ## note - with long reads, can have reads split into
    ## multiple parts
    qname <- c("read1", "read2", "read1", "read1")
    strand <- c("-", "+", "-", "-")
    sampleReadsRaw <- list(
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(0, ncol = 1),
            0
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1),
            1
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(2, ncol = 1),
            2
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(5, ncol = 1),
            3
        )
    )
    expected_sampleReads <- list(
        list(
            2, 0,
            matrix(c(-10, -10, -10), ncol = 1),
            matrix(c(0, 2, 5), ncol = 1)            
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1)            
        )
    )

    
    out <- merge_reads_from_sampleReadsRaw(
        sampleReadsRaw,
        qname,
        strand
    )
    
    sampleReads <- out$sampleReads
    sampleReadsInfo <- out$sampleReadsInfo
    expect_equal(
        sampleReads,
        expected_sampleReads
    )
    expect_equal(
        sampleReadsInfo,
        data.frame(
            qname = c("read1", "read2"),
            strand = c("-", "+")
        )
    )
    
})



test_that("sampleReadsRaw can be converted to sampleReads and also strand information obtained", {

    ## note - with long reads, can have reads split into
    ## multiple parts
    qname <- c("read3", "read1", "read2", "read3")
    strand <- c("-", "+", "*", "-")
    sampleReadsRaw <- list(
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(0, ncol = 1),
            0
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1),
            1
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(2, ncol = 1),
            2
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(5, ncol = 1),
            3
        )
    )
    expected_sampleReads <- list(
        list(
            1, 0,
            matrix(c(-10, -10), ncol = 1),
            matrix(c(0, 5), ncol = 1)            
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(1, ncol = 1)            
        ),
        list(
            0, 0,
            matrix(-10, ncol = 1),
            matrix(2, ncol = 1)            
        )
    )
    
    out <- merge_reads_from_sampleReadsRaw(
        sampleReadsRaw,
        qname,
        strand
    )
    
    sampleReads <- out$sampleReads
    sampleReadsInfo <- out$sampleReadsInfo
    expect_equal(
        sampleReads,
        expected_sampleReads
    )
    expect_equal(
        sampleReadsInfo,
        data.frame(
            qname = c("read3", "read1", "read2"),
            strand = c("-", "+", "*")
        )
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
                  "GACCGA", "::::::")
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
            matrix(c(-25, 25, -25), ncol = 1),
            matrix(c(4, 5, 6), ncol = 1)
        )
    )

    
    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-8",
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
                  "AACCTT", "::::::"),
                c("r002", "0", chr, "11", "60",
                  "3M", "*", "0", "0",
                  "CCT", ":::")
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
            matrix(c(-25, 25, 25), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        ),
        list(
            1, NA,
            matrix(c(25, 25), ncol = 1),
            matrix(c(1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-7",
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
        T = as.integer(nrow(pos)),
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
                  "AAA", ":::"),
                c("r001", "0", chr, "11", "60",
                  "3M", "*", "0", "0",
                  "CCT", ":::")
            ), 
            chr
        )
    )
    expected_sample_reads <- list(
        list(
            2, 1,
            matrix(c(-25, 25, 25), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-6",
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
                  "AA", "::"),
                c("r001", "0", chr, "11", "60",
                  "2M", "*", "0", "0",
                  "CC", "::"),
                c("r001", "0", chr, "13", "60",
                  "2M", "*", "0", "0",
                  "TT", "::")
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
            matrix(c(-25, 25, 25), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-5",
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
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-4",
        inputdir = tempdir(),
        regionName = regionName,
        tempdir = tempdir(),
        chr = chr,
        chrStart = 1,
        chrEnd = 100,
        width_of_loading_window = 10 ## make difficult!
    )
    
    load(file_sampleReads(tempdir(), 1, regionName))
    sampleReads[[1]][[2]] <- NA ## manually blank out
    sampleReads[[2]][[2]] <- NA ## manually blank out    

    expect_equal(
        sampleReads,
        expected_sample_reads
    )

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
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-3",
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
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-2",
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
                  "ccCCgg", "::::::")
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
            matrix(c(25, 25, -25), ncol = 1),
            matrix(c(3, 4, 5), ncol = 1)
        )
    )

    
    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
        bam_files = bam_file,
        N = 1,
        sampleNames = "test-name-1",
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
    bqs <- paste0(sapply(58 + 0:10, function(x) rawToChar(as.raw(x))), collapse = "")
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                c("r001", "0", chr, "14", "60",
                  "4S4M3S", "*", "0", "0",
                  "ggccGGGGccc", bqs)
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
        T = as.integer(nrow(pos)),
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
                  "CCGG", "::::")
            ), 
            chr
        )
    )
    
    expected_sample_reads <- list(
        list(
            1, NA,
            matrix(c(25, -25), ncol = 1),
            matrix(c(4, 5), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
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



test_that("sample CRAM can be properly interpreted", {
    
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
                  "AACCAA", "::::::")
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
            matrix(c(-25, 25, -25), ncol = 1),
            matrix(c(4, 5, 6), ncol = 1)
        )
    )

    regionName <- "region-name"
    loadBamAndConvert(
        iBam = 1,
        L = as.integer(pos[, 2]),
        pos = pos,
        T = as.integer(nrow(pos)),
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
