### test pos
test_that("pos is properly loaded and interpreted", {
    posfile <- tempfile()
    pos <- make_posfile(posfile)
    loaded_pos <- get_and_validate_pos(posfile, chr = 1)
    expect_equal(loaded_pos, pos)
})

test_that("an error is thrown when posfile is supplied but does not exist", {
    posfile <- file.path(tempfile(), "not_a_file")
    expect_error(
        validate_posfile(posfile),
        paste0("Cannot find supplied posfile:", posfile)
    )
})

test_that("pos works for chromosome with chr in it", {
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = 3, chr = rep("chr1", 3))
    loaded_pos <- get_and_validate_pos(posfile, chr = "chr1")
    expect_equal(loaded_pos, pos)
})


test_that("pos with more than one chromosome throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        chr = c(1, 1, 2),
        n_snps = 3
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file needs to be unique to chromosome and you supplied chr=1 but chromosome 2 was observed in row 3"
    )
})


test_that("pos chr column disagreeing with chr variable throws an erro", {
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = 2, chr = c(10, 10))
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file needs to be unique to chromosome and you supplied chr=1 but chromosome 10 was observed in row 1"
    )
})


test_that("pos column 2 not being an integer throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        L = c(10, "A"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        paste0("pos file column 2 needs to be integer valued between 1 and ", .Machine$integer.max, " but in row 2 A was observed")
    )
})


test_that("pos out of order throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        n_snps = 4,
        L = c(9, 10, 5, 6),
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file column 2 needs to be sorted on position with increasing positions between rows but row number 2 has position 10 and row number 3 has position 5"
    )
})



test_that("pos file with two repeat positions throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        L = c(5, 10, 10, 11),
        n_snps = 4
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file column 2 needs to be sorted on position with increasing positions between rows but row number 2 has position 10 and row number 3 has position 10"
    )
})


test_that("pos file with indel in ref throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        refs = c("A", "AA"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file ref column entry AA in row 2 contains is not one or A, C, G or T. STITCH is only supported for bi-allelic SNPs"
    )
})

test_that("pos file with other character in ref throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        refs = c("A", "X"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file ref column entry X in row 2 contains is not one or A, C, G or T. STITCH is only supported for bi-allelic SNPs"
    )
})


test_that("pos file with indel in alt throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        alts = c("A", "AA"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file alt column entry AA in row 2 contains is not one or A, C, G or T. STITCH is only supported for bi-allelic SNPs"
    )
})

test_that("pos file with other character in alt throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        alts = c("A", "X"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file alt column entry X in row 2 contains is not one or A, C, G or T. STITCH is only supported for bi-allelic SNPs"
    )
})

test_that("pos file with same character twice an error", {
    posfile <- tempfile()
    pos <- make_posfile(
        posfile,
        refs = c("A", "C"),
        alts = c("A", "G"),
        n_snps = 2
    )
    expect_error(
        get_and_validate_pos(posfile, chr = 1),
        "pos file row 1 has reference base A which is the same as alternate base A, which is not a bi-allelic SNP"
    )
})



### test gen
test_that("gen returns null if no file is given", {
    expect_equal(get_and_validate_gen(""), NULL)
})


test_that("an error is thrown if genfile does not exist", {
    genfile <- file.path(tempfile(), "not_a_file")
    expect_error(validate_genfile(genfile), paste0("Cannot find supplied genfile:", genfile))
})


test_that("gen is properly loaded and interpreted", {
    genfile <- tempfile()
    gen <- make_genfile(genfile)
    loaded_gen <- get_and_validate_gen(genfile)
    expect_identical(loaded_gen, gen)
})

test_that("gen works if there is only 1 sample", {
    genfile <- tempfile()
    gen <- make_genfile(genfile, n_samples = 1)
    loaded_gen <- get_and_validate_gen(genfile)
    expect_identical(loaded_gen, gen)
})

test_that("gen file with no header throws an error", {
    genfile <- tempfile()
    gen <- make_genfile(genfile, include_header = FALSE)
    expect_error(
        get_and_validate_gen(genfile),
        "The header for the genfile is either invalid or missing. The first invalid entry is in position 1 and is entry 0"
    )
})

test_that("gen throws an error if entries other than 0, 1, 2 or NA are used", {
    genfile <- tempfile()
    gen <- make_genfile(genfile, include_header = FALSE, vals = c(5, 5))
    expect_error(
        get_and_validate_gen(genfile),
        "Ineligible entry for genfile column 1 at position 1 with entry 5. Acceptable entries are 0, 1, 2 or NA"
    )
})




### test phase
test_that("phase returns null if default input", {
    expect_equal(
        get_and_validate_phase(""),
        NULL
    )
})



test_that("phase is properly loaded and interpreted", {
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile)
    loaded_phase <- get_and_validate_phase(phasefile)
    expect_identical(loaded_phase, phase)
})

test_that("phase works if there is only 1 sample", {
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, n_samples = 1)
    loaded_phase <- get_and_validate_phase(phasefile)
    expect_identical(loaded_phase, phase)
})


test_that("phase throws an error if no header", {
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, include_header = FALSE)
    expect_error(
        get_and_validate_phase(phasefile),
        paste0("The header for the phasefile is either invalid or missing. The first invalid entry is in position 1 and is entry ", paste0(phase[1, 1, 1], "|", phase[1, 1, 2])),
        fixed = TRUE
    )
})

test_that("phase throws an error if bad split character", {
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, split_character = ":")
    expect_error(
        get_and_validate_phase(phasefile),
        paste0(
            "Unable to split column 1 of phasefile at position 1 with entry ",
            paste0(phase[1, 1, 1], ":", phase[1, 1, 2]),
            " due to lack of field separator |"
        ),
        fixed = TRUE
    )
})

test_that("phase throws an error if values other than 0 or 1 are used", {
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, vals = c(0, 1, 2), seed = 1)
    if (length(RNGkind()) > 2 && RNGkind()[3] == "Rejection") {
        ## 3.6.0 and beyond
        expected_error <- "The phasefile contains entries other than 0, 1 or NA. One such entry is in column 1 and row 2  with value 2|0"
    } else {
        expected_error <- "The phasefile contains entries other than 0, 1 or NA1. One such entry is in column 1 and row 4  with value 2|0"
    }
    expect_error(
        get_and_validate_phase(phasefile),
        expected_error,
        fixed = TRUE
    )
})



## all 3
test_that("pos, gen and phase can be loaded", {
    n_snps <- 3
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = n_snps)
    genfile <- tempfile()
    gen <- make_genfile(genfile, n_snps = n_snps)
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, n_snps = n_snps)
    out <- get_and_validate_pos_gen_and_phase(
        posfile,
        genfile,
        phasefile,
        chr = 1
    )
    expect_equal(out$nSNPs, n_snps)
    expect_equal(out$gen, gen)
    expect_equal(out$pos, pos)
    expect_equal(out$phase, phase)
})


test_that("gen and phase can be omitted", {
    n_snps <- 3
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = n_snps)
    out <- get_and_validate_pos_gen_and_phase(
        posfile,
        chr = 1
    )
    expect_equal(out$nSNPs, n_snps)
    expect_equal(out$pos, pos)
})

test_that("gen can be omitted", {
    n_snps <- 3
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = n_snps)
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, n_snps = n_snps)
    out <- get_and_validate_pos_gen_and_phase(
        posfile = posfile,
        phasefile = phasefile,
        chr = 1
    )
    expect_equal(out$nSNPs, n_snps)
    expect_equal(out$pos, pos)
    expect_equal(out$phase, phase)
})

test_that("phase can be omitted", {
    n_snps <- 3
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = n_snps)
    genfile <- tempfile()
    gen <- make_genfile(genfile, n_snps = n_snps)
    out <- get_and_validate_pos_gen_and_phase(
        posfile = posfile,
        genfile = genfile,
        chr = 1
    )
    expect_equal(out$nSNPs, n_snps)
    expect_equal(out$pos, pos)
    expect_equal(out$gen, gen)
})


test_that("pos different sized than gen throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = 4)
    genfile <- tempfile()
    gen <- make_genfile(genfile, n_snps = 3)
    expect_error(
        get_and_validate_pos_gen_and_phase(
            posfile = posfile,
            genfile = genfile,
            chr = 1
        ),
        "posfile and genfile must have the same number of SNPs but posfile has 4 SNPs and genfile has 3 SNPs (and 1 header row)",
        fixed = TRUE
    )
})

test_that("pos different sized than phase throws an error", {
    posfile <- tempfile()
    pos <- make_posfile(posfile, n_snps = 4)
    phasefile <- tempfile()
    phase <- make_phasefile(phasefile, n_snps = 5)
    expect_error(
        get_and_validate_pos_gen_and_phase(
            posfile = posfile,
            phasefile = phasefile,
            chr = 1
        ),
        "posfile and phasefile must have the same number of SNPs but posfile has 4 SNPs and phasefile has 5 SNPs (and 1 header row)",
        fixed = TRUE
    )
})


## test sample name matching
test_that("gen name is properly matched", {
    sampleNames <- paste0("samp", 1:10)
    genfile <- tempfile()
    gen <- make_genfile(
        genfile, n_samples = 3,
        header = c("samp1", "samp8", "samp5")
    )
    out <- match_gen_and_phase_to_samples(
        sampleNames,
        gen
    )
    expect_equal(
        out$highCovInLow,
        c(1, 8, 5)
    )
})

## test sample name matching
test_that("gen name is properly matched with one sample", {
    sampleNames <- paste0("samp", 1:10)
    genfile <- tempfile()
    gen <- make_genfile(
        genfile, n_samples = 1,
        header = c("samp7")
    )
    out <- match_gen_and_phase_to_samples(
        sampleNames,
        gen
    )
    expect_equal(
        out$highCovInLow,
        7
    )
})

test_that("phase name is properly matched", {
    sampleNames <- paste0("samp", 1:10)
    phasefile <- tempfile()
    phase <- make_phasefile(
        phasefile, n_samples = 3,
        header = c("samp2", "samp7", "samp3")
    )
    out <- match_gen_and_phase_to_samples(
        sampleNames,
        phase = phase
    )
    expect_equal(
        out$samples_with_phase,
        c(2, 7, 3)
    )
})

test_that("phase name is properly matched with one sample", {
    sampleNames <- paste0("samp", 1:10)
    phasefile <- tempfile()
    phase <- make_phasefile(
        phasefile, n_samples = 1,
        header = c("samp9")
    )
    out <- match_gen_and_phase_to_samples(
        sampleNames,
        phase = phase
    )
    expect_equal(
        out$samples_with_phase,
        9
    )
})


test_that("can test-drive NA phase", {
    
    sampleNames <- paste0("samp", 1:10)
    phasefile <- tempfile()
    n_snps <- 10
    write_row_as_NA <- array(FALSE, n_snps)
    write_row_as_NA[5] <- TRUE
    phaseIN <- make_phasefile(
        phasefile,
        n_samples = 1,
        header = c("samp9"),
        write_row_as_NA = write_row_as_NA
    )
    phaseIN[write_row_as_NA, , ] <- NA
    ## 
    phaseOUT <- get_and_validate_phase(phasefile)
    expect_equal(phaseIN, phaseOUT)

})
