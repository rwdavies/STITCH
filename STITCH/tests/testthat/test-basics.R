test_that("reference legend file with duplicate entries throws an error", {
    reference_legend <- data.frame(
        position = c(10, 10, 20),
        a0 = c("A", "A", "T"),
        a1 = c("A", "A", "G")
    )
    expect_error(
        validate_reference_legend(reference_legend),
        "There are 1 duplicate row ids. One such example is 10 A A"
    )
})


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
