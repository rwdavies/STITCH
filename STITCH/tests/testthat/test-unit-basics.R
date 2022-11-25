test_that("reference legend file with duplicate entries throws an error", {
    reference_legend <- data.frame(
        position = c(10, 10, 20),
        a0 = c("A", "A", "T"),
        a1 = c("A", "A", "G")
    )
    expect_error(
        validate_reference_legend(reference_legend, "dummy.txt"),
        "The reference legend file dummy.txt  column 2 needs to be sorted on position with increasing positions between rows but row number 1 has position 10 and row number 2 has position 10"
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
    alleleCount[, c(1, 3)] <- NA
    print_allele_count(alleleCount, N = 10) 
})

test_that("can print allele count properly for normal results", {
    set.seed(1)
    alleleCount <- array(runif(30), c(10, 3))
    alleleCount[, c(1, 3)] <- NA    
    print_allele_count(alleleCount, N = 10) 
})

test_that("can print allele count for variable results", {
    alleleCount <- array(NA, c(100, 3))
    alleleCount[, 2] <- seq(1, 1e6, length.out = 100)
    print_allele_count(alleleCount, N = 10) 
})


test_that("can validate method", {
    expect_equal(validate_method("diploid"), NULL)
    expect_equal(validate_method("pseudoHaploid"), NULL)
    expect_equal(validate_method("diploid-inbred"), NULL)
    expect_error(validate_method("haploid"))
})

test_that("can validate K", {
    
    expect_null(validate_K(1))
    expect_null(validate_K(10))
    expect_null(validate_K(100))
    
    expect_error(validate_K(0))
    expect_error(validate_K(0.5))
    expect_error(validate_K(4.2))
    expect_error(validate_K(-10))
    expect_error(validate_K("10"))
    expect_error(validate_K("word"))
    

})
