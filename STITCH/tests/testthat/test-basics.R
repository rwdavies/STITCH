test_that("dependency checker works", {

    expect_error(
        check_program_dependency("not_a_program"),
        paste0(
            "The program not_a_program is not available in the PATH. ",
            "STITCH requires not_a_program to function. ",
            "Please make not_a_program available from the PATH"
        )
    )

})

test_that("vcf_output_name throws an error when too short", {
    expect_error(validate_vcf_output_name("A"), "vcf_output_name must have at least 8 characters and end with .vcf.gz, and you have supplied vcf_output_name:A")
})

test_that("vcf_output_name throws an error if no .vcf.gz", {
    expect_error(validate_vcf_output_name("A.vcf"), "vcf_output_name must have at least 8 characters and end with .vcf.gz, and you have supplied vcf_output_name:A.vcf")
})

test_that("vcf_output_name can be OK", {
    expect_equal(validate_vcf_output_name("A.vcf.gz"), NULL)
})

test_that("if vcf_output_name does not have a folder, it is placed in outputdir", {
    vcf_output_name <- "jimmy.vcf.gz"
    outputdir <- "/path/to/folder"
    regionName <- "20.10.20"
    output_vcf <- get_output_vcf(vcf_output_name, outputdir, regionName)
    expect_equal(
        output_vcf,
        file.path(outputdir, "jimmy.vcf.gz")
    )
})

test_that("if vcf_output_name has a folder, this is the path to the file", {
    vcf_output_name <- file.path("path", "to", "jimmy.vcf.gz")
    outputdir <- "/path/to/folder"
    regionName <- "20.10.20"
    output_vcf <- get_output_vcf(vcf_output_name, outputdir, regionName)
    expect_equal(
        output_vcf,
        vcf_output_name
    )
})

test_that("if vcf_output_name is NULL, vcf_output_name is constructed in the default way", {
    vcf_output_name <- NULL
    outputdir <- "/path/to/folder"
    regionName <- "20.10.20"
    output_vcf <- get_output_vcf(vcf_output_name, outputdir, regionName)
    expect_equal(
        output_vcf,
        file.path(
            outputdir,
            paste0("stitch.", regionName, ".vcf.gz")
        )
    )
})


test_that("an error is thrown when writing to an unwrittable directory", {

    unwritable_outputdir <- tempfile("test")
    dir.create(unwritable_outputdir)
    system(paste0("chmod 000 ", unwritable_outputdir))
    unwritable_outputdir <- file.path(unwritable_outputdir, "outputdir")

    expect_error(
        initialize_directories(
            tempdir = tempdir(),
            keepTempDir = FALSE,
            outputdir = unwritable_outputdir
        ),
        paste0(
            "Unable to make the required directory ",
            unwritable_outputdir,
            " while running. You can try re-starting STITCH, but if the problem consists, please contact your system administrator."
        )
    )

})


test_that("an error is thown when bamlist does not exist", {

    bamlist_that_does_not_exist <- tempfile()
    expect_error(
        validate_bamlist_and_cramlist_for_input_generation(
            regenerateInput = TRUE,
            bamlist = bamlist_that_does_not_exist
        ),
        paste0("Cannot find bamlist:", bamlist_that_does_not_exist)
    )

})


test_that("an error is thown when cramlist does not exist", {

    cramlist_that_does_not_exist <- tempfile()
    ref <- tempfile()
    system(paste0("touch ", ref))

    expect_error(
        validate_bamlist_and_cramlist_for_input_generation(
            regenerateInput = TRUE,
            cramlist = cramlist_that_does_not_exist,
            reference = ref
        ),
        paste0("Cannot find cramlist:", cramlist_that_does_not_exist)
    )

})


test_that("can print allele count properly for empty results", {
    alleleCount <- array(0, c(10, 3))
    print_allele_count(alleleCount, N = 10)
})

test_that("can print allele count properly for normal results", {
    set.seed(1)
    alleleCount <- array(runif(30), c(10, 3))
    print_allele_count(alleleCount, N = 10)
})

test_that("can print allele count for variable results", {
    alleleCount <- array(0, c(100, 3))
    alleleCount[, 2] <- seq(1, 1e6, length.out = 100)
    print_allele_count(alleleCount, N = 10)
})

test_that("allele count throws an error if there is an internal problem", {
    alleleCount <- array(NA, c(5, 3))
    alleleCount[, 2:3] <- c(2, 3)
    expect_error(
        check_allele_count_OK(alleleCount),
        "An internal problem with STITCH occured. There are NAs in the allele count. This may indicate a problem converting the input BAM/CRAM data into the internal STITCH format. The first error occured for SNP number 1 with entry NA, 2, 3"
    )
})
