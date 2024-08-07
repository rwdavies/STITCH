if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/proj/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

}


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



test_that("STITCH diploid works under default parameters", {

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
        outputBlockSize = 3
    )

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




## afterwards, everything should test both bgen and bgvcf

test_that("STITCH diploid works with multiple cores", {

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
            nCores = 4,
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


test_that("STITCH diploid works with completely split reads", {

    for(output_format in c("bgvcf", "bgen")) {

        set.seed(99)        
        outputdir <- make_unique_tempdir()
        output_filename <- paste0("jjj", extension[output_format])        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            output_format = output_format,
            nCores = 3,
            readAware = FALSE,
            output_filename = output_filename
        )

        check_output_against_phase(
            file = file.path(outputdir, output_filename),            
            data_package,
            output_format,
            tol = 0.2,
            min_info = 0.84 ## fluke so much poorer performance 3.6.0
        )
        
    }

})

test_that("STITCH diploid works with genfile with entirely NA columns", {

    outputdir <- make_unique_tempdir()

    local_genfile <- tempfile()
    gen <- read.table(data_package$genfile, header = TRUE)
    gen[, 1] <- NA
    write.table(gen, file = local_genfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = local_genfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        outputBlockSize = 3
    )

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


test_that("STITCH diploid works under default parameters when outputdir has a space or a tilda in it", {

    ## not 100% sure I can write to this in all systems, so check it first
    outputdirs <- c(
        file.path(make_unique_tempdir(), "wer wer2"),
        tempfile(tmpdir = "~/")
    )
    ## check can write to tilda
    dir.create(outputdirs[2])
    testfile <- file.path(outputdirs[2], "temp.txt")
    cat("", file = testfile)
    expect_equal(file.exists(testfile), TRUE)
    unlink(outputdirs[2])
        
    for(outputdir in outputdirs) {

        for(output_format in c("bgvcf", "bgen")) {
        
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
            
            check_output_against_phase(
                file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
                data_package,
                output_format,
                which_snps = NULL,
                tol = 0.2
            )

        }
        unlink(outputdir, recursive = TRUE)
    }

})


test_that("STITCH diploid works under default parameters when tempdir is relative", {

    outputdir <- make_unique_tempdir()
    setwd(tempdir())    
    temp_folder <- "temp"
    system(paste0("mkdir temp"))
        
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,        
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        tempdir = temp_folder
    )
    
    check_output_against_phase(
        file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        data_package,
        output_format,
        which_snps = NULL,
        tol = 0.2
    )
    
    unlink(outputdir, recursive = TRUE)

})






test_that("STITCH diploid works under default parameters with nCores = 40 and N = 25", {
    
    n_snps <- 10
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
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package25$chr, extension[output_format])),
            data_package25,
            output_format,
            which_snps = NULL,
            tol = 0.2,
            min_info = 0.85
        )

    }

})


test_that("STITCH diploid works under default parameters with N = 25 and nCores = 40", {

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

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package25$chr, extension[output_format])),
            data_package25,
            output_format,
            which_snps = NULL,
            tol = 0.2,
            min_info = 0.85
        )

    }


})



test_that("STITCH diploid works with regionStart, regionEnd and buffer", {

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
            regionStart = 3,
            regionEnd = 7,
            buffer = 1,
            output_format = output_format
        )

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



test_that("STITCH diploid works under default parameters using CRAM files", {

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










test_that("STITCH can get sample names from a file", {

    ## here need something that gets re-downsampled
    ## need more snps to make this happen easily
    outputdir <- make_unique_tempdir()

    sampleNames_file <- tempfile()
    N <- length(data_package$sample_names)
    newNames <- sapply(1:N, function(i) return(paste0(sample(letters, 10), collapse = "")))
    write.table(
        matrix(newNames, ncol = 1),
        file = sampleNames_file,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    for(output_format in c("bgvcf", "bgen")) {

        sink("/dev/null")
        set.seed(10)
        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            output_format = output_format,
            sampleNames_file = sampleNames_file
        )
        
        file <- file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format]) )
        if (output_format == "bgvcf") {
            x <- system(paste0("gunzip -c ", shQuote(file), " | grep CHROM"), intern = TRUE)
            output_sample_names <- strsplit(x, "\t")[[1]][-c(1:9)]
        } else {
            output_sample_names <- rrbgen::rrbgen_load_samples(file)
        }

        expect_equivalent(output_sample_names, newNames)
    }
    
})


test_that("STITCH works in a situation with grid, buffer, outputBlockSize, etc", {

    ## lots of samples, K, varying coverages, has grid, etc
    set.seed(9952)    
    n_snps <- 50
    reads_span_n_snps <- 4
    chr <- 1
    K <- 2
    extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")
    phasemaster <- matrix(
        c(
            sample(c(0, 1), n_snps, replace = TRUE),
            sample(c(0, 1), n_snps, replace = TRUE)
        ),
        ncol = K
    )
    phasemaster[10, 1] <- 1
    phasemaster[10, -1] <- 0
    phasemaster[15, ] <- 0
    regionStart <- 4
    regionEnd <- 36
    L <- 1:50
    buffer <- 2
    n_samples <- 30
    n_reads <- round(seq(5 * n_samples, 0, length.out = n_samples))
    data_package <- make_acceptance_test_data_package(
        n_samples = n_samples,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 6,
        chr = "chr5",
        K = K,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster,
        L = L
    )
    new_genfile <- tempfile()
    write.table(
        read.table(data_package$genfile, header = TRUE)[, c(5, 8, 10, 12)],
        file = new_genfile,
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )

    for(output_format in c("bgvcf", "bgen")) {    

        outputdir <- make_unique_tempdir()        
        output_filename <- paste0("jimmy", extension[output_format])
        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = new_genfile,
            outputdir = outputdir,
            output_filename = output_filename,
            output_format = output_format,
            K = 2,
            nGen = 100,
            nCores = 2,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            gridWindowSize = 3,
            outputSNPBlockSize = 8,
            inputBundleBlockSize = 5
        )

        which_snps <- which((regionStart <= L) & (L <= regionEnd))
        who <- which(n_reads >= 10 )
        check_output_against_phase(
            file = file.path(outputdir, output_filename),
            data_package = data_package,
            output_format = output_format,
            which_snps = which_snps,
            tol = 0.2,
            who = who,
            min_info = 0.945
        )
        
    }

})

test_that("STITCH can generate interim plots", {

    for(method in get_available_methods()) {

        outputdir <- make_unique_tempdir()
        chr <- data_package$chr
        niterations <- 5
        STITCH(
            chr = chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            method = method,
            nCores = 1,
            niterations = niterations,
            plotHapSumDuringIterations = TRUE,
            shuffleHaplotypeIterations = NA,
	    refillIterations = NA
        )

        ## only care about plumbing and making plots
        d <- dir(file.path(outputdir, "plots"))
        expect_equal(length(grep(paste0("hapSum.", chr, ".iteration."), d)), niterations - 1)
        expect_equal(length(grep(paste0("alphaMat.", chr, ".iteration."), d)), 2 * (niterations - 1))

    }

})

test_that("STITCH works when using a tempdir for writing useTempdirWhileWriting", {

    outputdir <- make_unique_tempdir()

    for(inputBundleBlockSize in c(NA, 5)) {
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,
            outputdir = outputdir,
            K = 2,
            nGen = 100,
            nCores = 1,
            useTempdirWhileWriting = TRUE,
            inputBundleBlockSize = inputBundleBlockSize
        )
        
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
    }

})

test_that("STITCH can rebundle inputs", {

    outputdir <- make_unique_tempdir()
    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 2,        
        inputBundleBlockSize = 2
    )
    unlink(file.path(outputdir, paste0("stitch.chr", chr, ".vcf.gz")))

    STITCH(
        chr = data_package$chr,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 3,
        inputBundleBlockSize = 4,
        regenerateInput = FALSE,
        originalRegionName = data_package$chr,
        regenerateInputWithDefaultValues = TRUE
    )
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


test_that("STITCH can avoid rebundling inputs when only one bundle", {

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
        inputBundleBlockSize = length(data_package$bam_files)
    )
    unlink(file.path(outputdir, paste0("stitch.chr", chr, ".vcf.gz")))

    STITCH(
        chr = data_package$chr,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        K = 2,
        nGen = 100,
        nCores = 3,
        inputBundleBlockSize = 4,
        regenerateInput = FALSE,
        originalRegionName = data_package$chr,
        regenerateInputWithDefaultValues = TRUE
    )
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
