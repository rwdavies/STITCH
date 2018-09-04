n_snps <- 10
reads_span_n_snps <- 6
chr <- 1
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
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
    seed = 3,
    chr = "chr5",
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster
)




test_that("STITCH diploid works with snap to grid", {

    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

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
            outputBlockSize = 3,
            gridWindowSize = 2,
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




test_that("STITCH diploid works with snap to grid with a buffer", {

    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        set.seed(1257)
        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,        
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            outputBlockSize = 3,
            gridWindowSize = 2,
            regionStart = 3,
            regionEnd = 8,
            buffer = 1,
            output_format = output_format
        )
        
        regionName <- paste0(data_package$chr, ".", 3, ".", 8)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = 3:8,
            tol = 0.2
        )

    }

})


test_that("STITCH pseudoHaploid works with grid", {

    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        set.seed(1301)

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,        
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            gridWindowSize = 3,
            method = "pseudoHaploid",
            output_format = output_format
        )

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.5, ## not good!
            min_info = 0.6
        )

    }

})


test_that("STITCH diploid works with snap to grid with downsampleToCov", {

    ## here need something that gets re-downsampled
    ## need more snps to make this happen easily
    outputdir <- make_unique_tempdir()

    n_snps <- 100
    phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    data_package <- make_acceptance_test_data_package(
        n_samples = 4,
        n_snps = n_snps,
        n_reads = n_snps * 4,
        seed = 1355,
        chr = "chr5",
        K = 2,
        reads_span_n_snps = 2,
        phasemaster = phasemaster
    )
    
    for(output_format in c("bgvcf", "bgen")) {

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
            outputBlockSize = 3,
            gridWindowSize = 5,
            output_format = output_format,
            downsampleToCov = 2 ## insanely low!
        )

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 2, ## basically, disable. only care about whether it worked with this toy example
            min_info = 0
        )

        ## here additional testing is performed
        if (output_format == "bgvcf") {
            vcf <- read.table(file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])))
            ## want eaf and info-score
            results <- t(sapply(1:nrow(vcf), function(i_snp) {
                x <- strsplit(as.character(vcf[i_snp, 8]), ";")[[1]]
                sapply(c("EAF", "INFO_SCORE"), function(z) {
                    y <- grep(z, x)
                    expect_true(length(y) == 1)
                    z2 <- as.numeric(strsplit(x[y], paste0(z, "="))[[1]][2])
                    return(z2)
                })
            }))
            no_freq <- (results[, "EAF"] <= 0.001) | (results[, "EAF"] >= (1-0.001))
            ## check that only those with no frequency get an INFO score of 1
            expect_equal(sum(results[no_freq, "INFO_SCORE"] != 1), 0)
            expect_equal(sum(results[!no_freq, "INFO_SCORE"] == 1), 0)            
        }
        
    }

})
