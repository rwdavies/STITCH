if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    ## dir <- "~/Google Drive/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    o <- sapply(dir(pattern = "*R"), source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))
    

}

run_only_one_acceptance_test <- TRUE
run_acceptance_tests <- TRUE

make_unique_tempdir <- function() {
    ## make every folder have a space!
    x <- tempfile(pattern = "folder", fileext = "wer wer2")
    dir.create(x)
    return(x)
}

n_snps <- 10
reads_span_n_snps <- 6
chr <- 1
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")

if (run_acceptance_tests | run_only_one_acceptance_test) {
    phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemaster[2, ] <- c(1, 0)
    phasemaster[3, ] <- c(0, 0)    
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
}

if (run_acceptance_tests) {
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        chr = chr
    )
    phasemasterC <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemasterC[1, ] <- c(1, 0)
    phasemasterC[4, ] <- c(0, 0)        
    phasemasterC[5, ] <- c(1, 0)    
    phasemasterC[6, ] <- c(1, 0)    
    data_package_crams <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 3,
        chr = chr,
        K = 2,
        phasemaster = phasemasterC,
        reads_span_n_snps = reads_span_n_snps,        
        use_crams = TRUE
    )
    chr <- "X"
    phasemasterX <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemasterX[3, ] <- c(1, 0)
    phasemasterX[5, ] <- c(0, 0)
    phasemasterX[6, ] <- c(1, 0)    
    data_packageX <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 3,
        chr = chr,
        K = 2,
        reads_span_n_snps = reads_span_n_snps,    
        phasemaster = phasemasterX
    )
    refpackX <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 4,
        reference_populations = c("CEU", "GBR", "CHB"),
        chr = data_packageX$chr
    )
}


test_that("STITCH diploid works under default parameters", {

    skip_test_if_TRUE(run_acceptance_tests | run_only_one_acceptance_test)
    outputdir <- make_unique_tempdir()
    
    sink("/dev/null")

    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        outputBlockSize = 3
    )

    sink()

    vcf <- read.table(
        file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH diploid can write to bgen", {

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()
    
    sink("/dev/null")

    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,        
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        output_format = "bgen"
    )

    sink()

    out_gp <- rrbgen::rrbgen_load(bgen_file = file.path(outputdir, paste0("stitch.", data_package$chr, ".bgen")))

    expect_equal(as.character(out_gp$var_info[, "chr"]), as.character(data_package$pos[, "CHR"]))
    expect_equal(as.character(out_gp$var_info[, "position"]), as.character(data_package$pos[, "POS"]))
    expect_equal(as.character(out_gp$var_info[, "ref"]), as.character(data_package$pos[, "REF"]))
    expect_equal(as.character(out_gp$var_info[, "alt"]), as.character(data_package$pos[, "ALT"]))        
    check_bgen_gp_against_phase(
        gp = out_gp$gp,
        phase = data_package$phase,
        tol = 0.2
    )

    var_info <- read.table(
        file.path(outputdir, paste0("stitch.", data_package$chr, ".bgen.per_snp_stats.txt.gz")),
        header = TRUE
    )
    expect_equal(nrow(var_info), 10)
    expect_equal((var_info[, "INFO_SCORE"] > 0.98), rep(TRUE, nrow(var_info)))    
    
})



## afterwards, everything should test both bgen and bgvcf

test_that("STITCH diploid works with multiple cores", {

    skip_test_if_TRUE(run_acceptance_tests | run_only_one_acceptance_test)
    
    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")

        outputdir <- make_unique_tempdir()        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 4,
            output_format = output_format
        )

        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )
        
    }

})

test_that("STITCH diploid works under default parameters when outputdir has a space in it", {

    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        sink("/dev/null")

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
            outputBlockSize = 3,
            output_format = output_format            
        )
    
        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )
    }

})




test_that("STITCH diploid works under default parameters with nCores = 40 and N = 25", {
    
    skip_test_if_TRUE(run_acceptance_tests)

    phasemaster25 <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
    phasemaster25[3, ] <- c(1, 0)
    phasemaster25[7, ] <- c(0, 1)    
    phasemaster25[10, ] <- c(1, 1)    
    data_package25 <- make_acceptance_test_data_package(
        n_samples = 25,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = 2,
        phasemaster = phasemaster25
    )
    
    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()
    
        sink("/dev/null")
        set.seed(10)

        STITCH(
            chr = data_package25$chr,
            bamlist = data_package25$bamlist,
            posfile = data_package25$posfile,
            genfile = data_package25$genfile,
            outputdir = outputdir,
            K = 4,
            nGen = 100,
            nCores = 40,
            output_format = output_format
        )

        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package25$chr, extension[output_format])),
            data_package25,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})


test_that("STITCH diploid works under default parameters with N = 25 and nCores = 40", {

    skip_test_if_TRUE(run_acceptance_tests)
    data_package25 <- make_acceptance_test_data_package(
        n_samples = 25,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = "2",
        K = 2,
        phasemaster = phasemaster
    )

    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(319)
        STITCH(
            chr = data_package25$chr,
            bamlist = data_package25$bamlist,
            posfile = data_package25$posfile,
            genfile = data_package25$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 40,
            output_format = output_format
        )

        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package25$chr, extension[output_format])),
            data_package25,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }


})



test_that("STITCH diploid works with regionStart, regionEnd and buffer", {

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")
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
            regionStart = 3,
            regionEnd = 7,
            buffer = 1,
            output_format = output_format
        )

        sink()

        regionName <- paste0(data_package$chr, ".", 3, ".", 7)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = 3:7,
            tol = 0.2
        )

    }

})


test_that("STITCH pseudoHaploid works under default parameters", {
    
    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        sink("/dev/null")

        outputdir <- file.path(make_unique_tempdir(), "wer wer2")    
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
        
        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.5
        )

    }
    
})


test_that("STITCH pseudoHaploid works with switchModelIteration", {
    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        sink("/dev/null")

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

        sink()

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
    
    skip_test_if_TRUE(run_acceptance_tests)

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
        
        sink("/dev/null")

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

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package3$chr, extension[output_format])),
            data_package3,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})


test_that("STITCH can initialize with reference data", {

    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        sink("/dev/null")

        outputdir <- make_unique_tempdir()
        set.seed(478)
        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,                
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            K = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )
        
        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }


})


test_that("STITCH can initialize with reference data with defined regionStart and regionEnd", {

    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        sink("/dev/null")

        outputdir <- make_unique_tempdir()
        set.seed(519)
        
        regionStart <- 3
        regionEnd <- 4
        buffer <- 1
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,        
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            K = 2,
            nGen = 100,
            nCores = 1,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            output_format = output_format
        )

        sink()

        which_snps <- regionStart:regionEnd
        regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = which_snps,
            tol = 0.2
        )

        load(file.path(outputdir, "RData", paste0("hwe.",regionName,".RData")))

        ## check hweCount - looks OK!
        ## probably just need this the once?
        g <- data_package$phase[3:4, , 1] + data_package$phase[3:4, , 2]
        expect_equal(rowSums(g == 0), hweCount[, 1])
        expect_equal(rowSums(g == 1), hweCount[, 2])
        expect_equal(rowSums(g == 2), hweCount[, 3])
    }

})

test_that("STITCH can initialize with reference data for certain populations", {
    
    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(578)
    
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_sample_file = refpack$reference_sample_file,
            reference_populations = refpack$reference_populations[1],
            K = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})


test_that("STITCH can initialize with reference data for certain populations for defined regionStart and regionEnd", {

    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(621)
        
        regionStart <- 5
        regionEnd <- 7
        buffer <- 1
        regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)

        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            reference_sample_file = refpack$reference_sample_file,
            reference_populations = refpack$reference_populations[1],
            K = 2,
            nGen = 100,
            nCores = 1,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            output_format = output_format
        )

        sink()
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = regionStart:regionEnd,
            tol = 0.2
        )

    }
        

})



test_that("STITCH can initialize on chromosome X with reference data", {
    
    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(621)
        
        STITCH(
            chr = data_packageX$chr,
            bamlist = data_packageX$bamlist,
            posfile = data_packageX$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackX$reference_haplotype_file,
            reference_legend_file = refpackX$reference_legend_file,
            K = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_packageX$chr, extension[output_format])),
            data_packageX,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})


test_that("STITCH can initialize on chromosome X with reference data with defined regionStart and regionEnd", {
    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(621)

        regionStart <- 2
        regionEnd <- 6
        buffer <- 1
        
        STITCH(
            chr = data_packageX$chr,
            bamlist = data_packageX$bamlist,
            posfile = data_packageX$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackX$reference_haplotype_file,
            reference_legend_file = refpackX$reference_legend_file,
            K = 2,
            nGen = 100,
            nCores = 1,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            output_format = output_format
        )

        sink()

        regionName <- paste0(data_packageX$chr, ".", regionStart, ".", regionEnd)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_packageX,
            output_format,
            which_snps = regionStart:regionEnd,
            tol = 0.2
        )

    }


})


test_that("STITCH can initialize on chromosome X with reference data for certain populations", {
    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    
        sink("/dev/null")

        set.seed(753)
        
        STITCH(
            chr = data_packageX$chr,
            bamlist = data_packageX$bamlist,
            posfile = data_packageX$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackX$reference_haplotype_file,
            reference_legend_file = refpackX$reference_legend_file,
            reference_sample_file = refpackX$reference_sample_file,
            reference_populations = refpackX$reference_populations[1],
            K = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_packageX$chr, extension[output_format])),
            data_packageX,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )
    }

})




test_that("STITCH can initialize on chromosome X with reference data for certain populations for defined regionStart and regionEnd", {

    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    

        set.seed(795)
        regionStart <- 3
        regionEnd <- 4
        buffer <- 1
        
        sink("/dev/null")

        STITCH(
            chr = data_packageX$chr,
            bamlist = data_packageX$bamlist,
            posfile = data_packageX$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackX$reference_haplotype_file,
            reference_legend_file = refpackX$reference_legend_file,
            reference_sample_file = refpackX$reference_sample_file,
            reference_populations = refpackX$reference_populations[1],
            K = 2,
            nGen = 100,
            nCores = 1,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            output_format = output_format
        )

        sink()

        regionName <- paste0(data_packageX$chr, ".", regionStart, ".", regionEnd)        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_packageX,
            output_format,
            which_snps = regionStart:regionEnd,
            tol = 0.2
        )

    }

})


test_that("STITCH diploid works under default parameters using CRAM files", {

    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()    

        set.seed(843)
        
        sink("/dev/null")

        STITCH(
            chr = data_package_crams$chr,
            cramlist = data_package_crams$cramlist,
            reference = data_package_crams$ref,
            posfile = data_package_crams$posfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            output_format = output_format
        )

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package_crams$chr, extension[output_format])),
            data_package_crams,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }

})





test_that("STITCH with generateInputOnly actually only generates input", {

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()
    sink("/dev/null")

    set.seed(10)
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        generateInputOnly = TRUE
    )

    expect_equal(
        FALSE,
        file.exists(file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")))
    )
    inputdir_contents <- dir(file.path(outputdir, "input"))
    expect_equal(
        sort(inputdir_contents),
        sort(paste0("sample.", 1:10, ".input.", data_package$chr, ".RData"))
    )


})



test_that("STITCH can generate input in VCF format", {

    skip_test_if_TRUE(run_acceptance_tests)    
    outputdir <- make_unique_tempdir()    

    sink("/dev/null")

    set.seed(10)
    n_snps <- 5
    chr <- 10
    phasemaster <- matrix(c(0, 1, 0, 0, 1, 0, 1, 1, 0, 0), ncol = 2)
    n_reads <- n_snps * 10
    n_samples <- 10
    data_package_local <- make_acceptance_test_data_package(
        n_samples = n_samples,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 1,
        chr = chr,
        K = 2,
        phasemaster = phasemaster
    )

    STITCH(
        chr = data_package_local$chr,
        bamlist = data_package_local$bamlist,
        posfile = data_package_local$posfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        outputInputInVCFFormat = TRUE
    )

    ## expect not to find normal VCF output, but new one
    expect_equal(
        FALSE,
        file.exists(file.path(outputdir, paste0("stitch.", data_package_local$chr, ".vcf.gz")))
    )
    expect_equal(
        TRUE,
        file.exists(file.path(outputdir, paste0("stitch.input.", data_package_local$chr, ".vcf.gz")))
    )

    ## since this is comparatively high coverage
    ## results should be OK
    vcf <- read.table(
        file.path(outputdir, paste0("stitch.input.", data_package_local$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    ## check read counts are OK
    ## in this, > 0 read counts as appropriate
    for (iSample in 1:n_samples) {
        rc <- sapply(strsplit(vcf[, iSample + 9], ":"), I)[2, ]
        rc2 <- sapply(strsplit(rc, ","), I) ## ref, alt
        genotype <- data_package_local$phase[, iSample, 1] + data_package_local$phase[, iSample, 2]
        ## ref counts - expect no alternate reads
        expect_equal(sum(as.integer(rc2[1, genotype == 0]) == 0), 0)
        expect_equal(sum(as.integer(rc2[2, genotype == 0]) > 0), 0)
        ## alt counts - expect both to never be 0
        expect_equal(sum(as.integer(rc2[2, genotype == 1]) == 0), 0)
        expect_equal(sum(as.integer(rc2[2, genotype == 1]) == 0), 0)
        ## hom alt counts - expect no reference reads
        expect_equal(sum(as.integer(rc2[1, genotype == 2]) > 0), 0)
        expect_equal(sum(as.integer(rc2[2, genotype == 2]) == 0), 0)
    }

})


test_that("STITCH can impute with reference panels with only 1 iteration if the initialize with reference data", {

    skip_test_if_TRUE(run_acceptance_tests)

    n_snps <- 5
    chr <- 10
    set.seed(10)
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 1,
        reference_populations = c("CEU", "GBR"),
        chr = chr
    )
    phasemaster <- refpack$reference_haplotypes
    data_package <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = 4,
        phasemaster = phasemaster
    )

    for(output_format in c("bgvcf", "bgen")) {
        
        outputdir <- make_unique_tempdir()
        sink("/dev/null")
        set.seed(1008)
    
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            refillIterations = NA,
            shuffleHaplotypeIterations = NA,
            K = 4,
            nGen = 100,
            nCores = 1,
            niterations = 1,
            output_format = output_format
        )
        
        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.3
        )

    }

})

test_that("STITCH errors if niterations=1 with reference panel and posfile is not the same as reference legend", {
    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()
    sink("/dev/null")

    n_snps <- 5
    chr <- 10
    set.seed(10)
    refpack <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = 1,
        reference_populations = c("CEU", "GBR"),
        chr = chr
    )
    K <- 4
    phasemaster <- matrix(c(rep(0, n_snps + 1), rep(1, n_snps + 1)), ncol = K)
    data_package <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps + 1,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = K,
        phasemaster = phasemaster
    )

    expect_error(
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            refillIterations = NA,
            shuffleHaplotypeIterations = NA,
            K = 4,
            nGen = 100,
            nCores = 1,
            niterations = 1
        ),
        "You have selected to use reference haplotypes with niterations=1, which requires exact matching of reference legend SNPs and posfile SNPs. However, posfile SNP with pos-ref-alt 6-A-G was not found in reference legend"
    )
    sink()

})


test_that("STITCH can initialize with reference data with snap to grid", {

    skip_test_if_TRUE(run_acceptance_tests)

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()
    
        sink("/dev/null")
        set.seed(1098)
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            reference_haplotype_file = refpack$reference_haplotype_file,
            reference_legend_file = refpack$reference_legend_file,
            K = 2,
            nGen = 100,
            nCores = 1,
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



test_that("STITCH diploid works with snap to grid", {

    skip_test_if_TRUE(run_acceptance_tests)
    
    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()
    
        sink("/dev/null")
        set.seed(10)
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
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




test_that("STITCH diploid works with snap to grid", {

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")
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

        sink()

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

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")
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
        
        sink()
        
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

    skip_test_if_TRUE(run_acceptance_tests)
    outputdir <- make_unique_tempdir()

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")

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

        sink()
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.5 ## not good!
        )

    }

})


test_that("STITCH diploid works with snap to grid with downsampleToCov", {

    ## here need something that gets re-downsampled
    ## need more snps to make this happen easily
    skip_test_if_TRUE(run_acceptance_tests)
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

        sink("/dev/null")
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

        sink()

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 2 ## basically, disable. only care about whether it worked with this toy example
        )
        
    }

})
