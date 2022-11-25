## 1 ST-J00101:180:HGVV3BBXY:5:1113:18030:28516
## 2 99
## 3 chr1
## 4 10058
## 5 12
## 6 42M1I9M2I54M42S
## 7 =
## 8 10027
## 9 74
## 10 CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCTAACCCCTAACCCCAACCCCTACCCCTAACCCAACCCCAACCCCCACCCCTACCCCTACCC
## 11 AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<JJJJJJ7F-FJJ-FAJFJ-<FFJJ-A-FJF--AFJJ--FJJJ--FJJJJ-<-AJJ7---77A---7-A---7-F-7-7<7---7AA---7AJ---<<<---7<

## rest follow TAG:TYPE:VALUE
## https://samtools.github.io/hts-specs/SAMv1.pdf

## XA:Z:chr1,+10349,4M1I103M42S,3;chr4,-190122800,35M7I12M1I42M1I52M,15;chr15,-101981001,29S18M1D50M1I6M1I15M1I29M,5;chr1,+10189,49M1I10M8D36M1I53M,18;chr3,+10643,74M1D19M5D57M,16;chr17,+136325,6M1D97M47S,4;chr13,-114354114,55S11M1D16M1D68M,2;chrX,-156029975,17M5D12M1D24M1I30M1D54M1I11M,15;chr22,-50808374,56S41M1I6M1I45M,2;chr1,-248946295,55S45M1D4M1D46M,3;chr1,-248946021,67S30M1I52M,1;chr1,-248946123,29S67M1I6M1I46M,8;chrX,-156030395,55S42M1I6M1I45M,3;chr5,+11739,66M1D12M1D25M47S,5;chr18,-80262892,69S35M1I45M,1;chr1,-248946155,54M1I45M17D4M1D46M,27;chr9,+10195,16M1I19M1D45M5I15M49S,7;chr20,-64287251,51S34M1I12M5D52M,8;chrX,+222348,42M1I6M1I40M60S,3;chr11,+191800,25S59M8D8M1I28M29S,10;chr4,-190122659,47S41M2D9M5D6M1I46M,10;chr17_GL383563v3_alt,+76325,6M1D97M47S,4;
## MC:Z:45S105M
## MD:Z:80T18T5
## PG:Z:MarkDuplicates
## RG:Z:Run180_Coriell-P2
## NM:i:5
## AS:i:80
## XS:i:93
## BX:Z:A08C56B87D20
## QX:Z:AAFFFJJJJJJJJ+AAFFFFJJJJJF
## RX:Z:TGTGCTATTGCCT+TCCATGCGGAAT
## BX: Chromium barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences. Use this for analysis.
## QX: Raw Chromium barcode read quality. Phred scores as reported by sequencer.
## RX: Raw Chromium barcode sequence. This read is subject to sequencing errors. Do not use for analysis.



test_that("can rescue bxtag", {

    ## so options
    ## 1: 2 reads, fine
    ## 2: 3 reads, fine
    ## 3: 2 reads, one is 00, so rescued
    ## 4: 2 reads, both 00, so kept 00
    ## 5: 2 reads, they disagree, so set to bad and missing

    qnameInteger_ord <- c(
        c(1, 1),
        c(2, 2, 2),
        c(3, 3),
        c(4, 4),
        c(5, 5)
    )
    bxtag_ord <- c(
        c("A01", "A01"),
        c("A02", "A02", "A02"),
        c("A01", "A00"),
        c("A00", "A00"),
        c("A05", "A06")
    )
    bxtag_ord_init <- bxtag_ord
    ## force make copy
    a <- bxtag_ord_init[1]
    bxtag_ord_init[1] <- "TEMP"
    bxtag_ord_init[1] <- a
    
    expected_bxtag_ord <- c(
        c("A01", "A01"),
        c("A02", "A02", "A02"),
        c("A01", "A01"),
        c("A00", "A00"),
        c("00", "00")
    )
    expected_bxtag_bad_ord <- array(FALSE, length(bxtag_ord))
    expected_bxtag_bad_ord[qnameInteger_ord == 4] <- TRUE    
    expected_bxtag_bad_ord[qnameInteger_ord == 5] <- TRUE

    bxtag_bad_ord <- rcpp_evaluate_bxtag(
        qnameInteger_ord = qnameInteger_ord,
        bxtag_ord = bxtag_ord
    )
    
    expect_equal(as.logical(bxtag_bad_ord), as.logical(expected_bxtag_bad_ord))
    expect_equal(expected_bxtag_ord, bxtag_ord)

})


test_that("can get bx tag into R", {

    chr <- 10
    pos <- cbind(
        CHR = c(chr, chr, chr),
        POS = c(9, 11, 13, 15, 17, 19),
        REF = rep("A", 6),
        ALT = rep("C", 6)
    )
    chrStart <- 1
    chrEnd <- 20
    f <- function(read_name, start, tag) {
        return(
            c(read_name, "0", chr, start, "60",
              "2M", "*", "0", "0",
              "AA", "::", tag)
        )
    }
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            list(
                f("r1", "10", "BX:Z:A1"),
                f("r1", "14", "BX:Z:A1"),
                f("r2", "16", "BX:Z:A1"),
                f("r3", "18", "BX:Z:A2"),
                f("r4", "18", NULL),
                f("r5", "18", "BX:Z:A3")
            ),
            chr
        )
    )
    
    useSoftClippedBases <- FALSE
    bqFilter <- 17
    iSizeUpperLimit <- 1000
    ref <- as.character(pos[, "REF"])
    alt <- as.character(pos[, "ALT"])
    nSNPs <- nrow(pos)
    L <- as.integer(pos[, "POS"])
    region <- paste0(chr, ":", chrStart, "-", chrEnd)
    file_name <- bam_file
    reference <- ""
    save_sampleReadsInfo <- FALSE
    
    out <- get_sampleReadsRaw_from_SeqLib(
        useSoftClippedBases = useSoftClippedBases,
        bqFilter = bqFilter,
        iSizeUpperLimit = iSizeUpperLimit,
        ref = ref,
        alt = alt,
        nSNPs = nSNPs,
        L = L,
        region = region,
        file_name = file_name,
        reference = reference,
        save_sampleReadsInfo = save_sampleReadsInfo,
        use_bx_tag = TRUE
    )

    expect_equal(out$bxtag, c("A1", "A1", "A1", "A2", "", "A3"))
    expect_equal(out$qname, c("r1", "r1", "r2", "r3", "r4", "r5"))
    
})





test_that("can make sampleReadsRaw incorporate bx tag properly for complicated example", {
    
    chr <- 10
    chrStart <- 1
    chrEnd <- 130
    L <- seq(chrStart, chrEnd, 2)
    pos <- cbind(
        CHR = rep(chr, length(L)),
        POS = L,
        REF = rep("A", length(L)),
        ALT = rep("C", length(L))
    )
    ## 
    ## read1, 3 read pairs, first one has 3 reads, middle pair fails iSizeUpperLimit difference, otherwise bx tag fine
    ## read2, 2 read pair, one pair has one with 00 tag but rescued
    ## read3, 1 read pairs, both have 00 tag, so together but no other
    ## read4, 2 read pairs, bx too large between them so not done

    ## 
    ## make expected output here
    ## 1: 0-based number of SNPs
    ## 2: 0-based central SNP
    ## 3: bq, - = ref, + = alt
    ## 4: 0-based position in pos
    ## note, don't bother checking central SNP
    make_expect_read <- function(bqs, poss) {
        stopifnot(length(bqs) == length(poss))
        return(list(
            length(bqs) - 1, NA,
            matrix(bqs, ncol = 1),
            matrix(match(poss, L) - 1, ncol = 1)
        ))
    }
    expected_sample_reads <- list(
        make_expect_read(c(-9, -10, -11, -14, -15), c(1, 5, 9, 11, 15)),
        make_expect_read(c(-16, -17, -18, -19), c(31, 39, 41, 43)),
        make_expect_read(c(-20, -21), c(61, 63)),
        make_expect_read(c(-22, -23), c(71, 77)),
        make_expect_read(c(-24, -25), c(121, 125))
    )
    
    bxTagUpperLimit <- 20
    iSizeUpperLimit <- 15
    f <- function(read_name, start, tag, bq) {
        return(
            c(read_name, "0", chr, start, "60",
              "2M", "*", "0", "0",
              "AA", paste0(rep(rawToChar(as.raw(bq + 33)), 2), collapse = ""), tag)
        )
    }

    to_sam <- list(
        f("r1_A01", "1", "BX:Z:A01", 9), ## first read here has three parts which is OK
        f("r1_A01", "5", "BX:Z:A01", 10), 
        f("r1_A01", "9", "BX:Z:A01", 11),
        f("r2_A01", "5", "BX:Z:A01", 12), ## these two filtered out
        f("r2_A01", "50", "BX:Z:A01", 13),## these two filtered out
        f("r3_A01", "11", "BX:Z:A01", 14),
        f("r3_A01", "15", "BX:Z:A01", 15),
        f("r1_A02", "31", "BX:Z:A02B01", 16),
        f("r1_A02", "39", "BX:Z:A00B00", 17), ## rescued
        f("r2_A02", "41", "BX:Z:A02B01", 18),
        f("r2_A02", "43", "BX:Z:A02B01", 19),
        f("r1_A03", "61", "BX:Z:A03B00", 20), ## both 00-tagged
        f("r1_A03", "63", "BX:Z:A03B00", 21),
        f("r1_A04", "71", "BX:Z:A04", 22),
        f("r1_A04", "77", "BX:Z:A04", 23),
        f("r2_A04", "121", "BX:Z:A04", 24), ## this is far beyond the 10 (11) of above
        f("r2_A04", "125", "BX:Z:A04", 25)
    )
    expected_bxtags <- c("A01", "A02B01", "A03B00", "A04", "A04") ## 2 of last one as broken up
    to_sam <- to_sam[order(as.numeric(sapply(to_sam, function(x) x[[4]])))]
    
    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            to_sam,
            chr
        )
    )

    set.seed(919)
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
        chrStart = chrStart,
        chrEnd = chrEnd,
         bqFilter = 5,
        use_bx_tag = TRUE,
        iSizeUpperLimit = iSizeUpperLimit,
        bxTagUpperLimit = bxTagUpperLimit,
        save_sampleReadsInfo = TRUE
    )

    load(file_sampleReads(tempdir(), 1, regionName))
    ## nuke central read for this check - randomness
    for(iRead in 1:length(sampleReads)) {
        sampleReads[[iRead]][[2]] <- NA
    }
    expect_equal(
         sampleReads,
         expected_sample_reads
     )

    load(file_sampleReadsInfo(tempdir(), 1, regionName))
    expect_equal(sum(sampleReadsInfo[, "bxtag"] != expected_bxtags), 0)
    
    ## save(sampleReads,          expected_sample_reads, file = "~/temp2.RData")

    ## load("~/temp2.RData")
    ## print("---------------------sampleReads-------------")
    ## iRead <- 2
    ## print(sampleReads[[iRead]])
    ## print("---------------------expected sampleReads-------------")    
    ## print(expected_sample_reads[[iRead]])
    
})















test_that("can resize appropriately", {

    chr <- 10
    chrStart <- 1
    chrEnd <- 130
    L <- seq(chrStart, chrEnd, 2)
    pos <- cbind(
        CHR = rep(chr, length(L)),
        POS = L,
        REF = rep("A", length(L)),
        ALT = rep("C", length(L))
    )
    ## 

    ## make three reads, all same bx tag
    make_expect_read <- function(bqs, poss) {
        stopifnot(length(bqs) == length(poss))
        return(list(
            length(bqs) - 1, 7,
            matrix(bqs, ncol = 1),
            matrix(match(poss, L) - 1, ncol = 1)
        ))
    }
    expected_sample_reads <- list(
        make_expect_read(c(-9, -9, -10, -10, -11, -11), 1 + c(10, 12, 14, 16, 18, 20))
    )
    
    bxTagUpperLimit <- 20
    iSizeUpperLimit <- 15
    f <- function(read_name, start, tag, bq) {
        return(
            c(read_name, "0", chr, start, "60",
              "4M", "*", "0", "0",
              "AAAA", paste0(rep(rawToChar(as.raw(bq + 33)), 4), collapse = ""), tag)
        )
    }

    to_sam <- list(
        f("r1_A01", "10", "BX:Z:A01", 9),
        f("r1_A01", "14", "BX:Z:A01", 10), 
        f("r1_A01", "18", "BX:Z:A01", 11)
    )
    expected_bxtags <- c("A01", "A01", "A01")
    to_sam <- to_sam[order(as.numeric(sapply(to_sam, function(x) x[[4]])))]

    bam_file <- make_simple_bam(
        file_stem = file.path(tempdir(), "simple"),
        sam = make_simple_sam_text(
            to_sam,
            chr
        )
    )

    set.seed(919)
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
        chrStart = chrStart,
        chrEnd = chrEnd,
         bqFilter = 5,
        use_bx_tag = TRUE,
        iSizeUpperLimit = iSizeUpperLimit,
        bxTagUpperLimit = bxTagUpperLimit,
        save_sampleReadsInfo = TRUE,
        maxnSNPInRead = 3
    )

    load(file_sampleReads(tempdir(), 1, regionName))

    ## check here
    expect_equal(
         sampleReads,
         expected_sample_reads
    )

})
