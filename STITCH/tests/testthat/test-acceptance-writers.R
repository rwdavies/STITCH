if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/Google Drive/STITCH/"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*R")
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


test_that("STITCH can generate input in VCF format", {

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



test_that("STITCH diploid can write to bgen", {

    outputdir <- make_unique_tempdir()
    
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
        output_format = "bgen"
    )

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
    
})



test_that("STITCH can write out haplotype probabilities", {

    set.seed(10)
    method <- "diploid"
    
    for(K in c(2, 3)) {

        outputdir <- make_unique_tempdir()        
        STITCH(
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            genfile = data_package$genfile,        
            outputdir = outputdir,
            K = K,
            nGen = 100,
            nCores = 1,
            output_format = "bgvcf",
            output_haplotype_dosages = TRUE,
            method = method
        )

        ## includes checking that the sum is correct
        check_output_against_phase(
            file = file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2,
            min_info = 0.85
        )

        ## check format line exists
        format_line <- system(paste0("gunzip -c ", shQuote(file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz"))), " | grep 'ID=HD'"), intern = TRUE)
        expect_equal(length(format_line), 1)
        ## check format line has proper K
        expect_equal(as.numeric(strsplit(strsplit(format_line, "Number=")[[1]][2], ",")[[1]][1]), K)

        file <- file.path(outputdir, paste0("stitch.", data_package$chr, ".vcf.gz"))
        vcf <- read.table(file)
        load(file.path(outputdir, "RData", "EM.all.chr5.RData"))
        transMatRate_t_H <- get_transMatRate(method = "diploid-inbred", sigmaCurrent = sigmaCurrent)
        transMatRate_t_D <- get_transMatRate(method = "diploid", sigmaCurrent = sigmaCurrent)    
        for(i_sample in 1:10) {
            ## run fbd normally, check OK
            load(file.path(outputdir, "input", paste0("sample.", i_sample, ".input.chr5.RData")))            
            q <- sapply(strsplit(as.character(vcf[, 9 + i_sample]), ":"), I)[4, ]
            q_t <- sapply(strsplit(q, ","), as.numeric)
            ## now, re-do
            fbsoL <- run_forward_backwards(
                sampleReads = sampleReads,
                method = method,
                K = K,
                priorCurrent = priorCurrent,
                alphaMatCurrent_t = alphaMatCurrent_t,
                eHapsCurrent_t = eHapsCurrent_t,
                transMatRate_t_H = transMatRate_t_H,
                transMatRate_t_D = transMatRate_t_D, 
                Jmax = 10,
                maxDifferenceBetweenReads = 1000,
                maxEmissionMatrixDifference = 1e6,
                output_haplotype_dosages = TRUE,
                grid = grid
            )$fbsoL
            ## check result
            expect_equal(max(abs(q_t - fbsoL[[1]]$gammaEK_t * 2)) > 1e-3, FALSE)
        }

    }
    
})
