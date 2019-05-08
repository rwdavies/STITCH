## suppressWarnings(RNGversion("3.5.0"))

n_snps <- 10
reads_span_n_snps <- 6
chr <- 1
n_reads <- round(5 / (reads_span_n_snps / n_snps)) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")

phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
phasemaster[2, ] <- c(1, 0)
phasemaster[3, ] <- c(0, 0)
phasemaster[4, ] <- c(1, 1)
phasemaster[7, ] <- c(1, 0)    
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 5,
    chr = "chr5",
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster
)


test_that("STITCH pseudoHaploid works under default parameters", {
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()
        set.seed(357)

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,                
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            method = "pseudoHaploid",
            output_format = output_format
        )
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.5,
            min_info = 0.6
        )

    }
    
})


test_that("STITCH pseudoHaploid works with outputSNPBlockSize", {

    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,                
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            method = "pseudoHaploid",
            output_format = output_format,
            gridWindowSize = 4,
            outputSNPBlockSize = 3            
        )
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.5,
            min_info = 0.60
        )

    }
    
})


test_that("STITCH pseudoHaploid works with switchModelIteration", {

    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- file.path(make_unique_tempdir(), "wer wer2")    
        set.seed(10)

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,        
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            method = "pseudoHaploid",
            switchModelIteration = 39,
            niterations = 40,
            output_format = output_format
        )

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})


test_that("STITCH pseudoHaploid works with a single sample and two cores", {

    phasemaster <- matrix(c(c(0, 0, 0), c(1, 1, 1)), ncol = 2)
    data_package3 <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = 3,
        n_reads = 4,
        seed = 1,
        chr = "chrWER",
        K = 2,
        phasemaster = phasemaster
    )
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- file.path(make_unique_tempdir(), "wer wer2")    
        set.seed(103)

        STITCH(
            chr = data_package3$chr,
            bamlist = data_package3$bamlist,
            posfile = data_package3$posfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 2,
            method = "pseudoHaploid",
            output_format = output_format
        )

        check_output_against_phase(
            file = file.path(outputdir, paste0("stitch.", data_package3$chr, extension[output_format])),
            data_package = data_package3,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})

