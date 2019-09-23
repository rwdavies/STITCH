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
original_chr <- "chr5"
chr <- original_chr
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")


phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
phasemaster[2, ] <- c(1, 0)
phasemaster[3, ] <- c(0, 0)
phasemaster[7, ] <- c(1, 0)
phasemaster_original <- phasemaster
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 3,
    chr = original_chr,
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster
)
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = chr,
    phasemaster = phasemaster
)



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
    chr = data_packageX$chr,
    phasemaster = phasemasterX
)


test_that("STITCH can initialize with reference data", {

    for(output_format in c("bgvcf", "bgen")) {

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

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
            data_package,
            output_format,
            which_snps = NULL,
            tol = 0.2
        )

    }


})


test_that("STITCH can initialize with reference data for S > 1", {

    for(output_format in c("bgvcf", "bgen")) {

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
            S = 3,
            nGen = 100,
            nCores = 1,
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


test_that("STITCH can initialize with reference data including genetic map", {

    output_format <- "bgvcf"
    outputdir <- make_unique_tempdir()
    set.seed(3393)

    STITCH(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        genfile = data_package$genfile,
        outputdir = outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        genetic_map_file = refpack$reference_genetic_map_file,
        K = 2,
        S = 3,
        gridWindowSize = 3,
        nGen = 100,
        nCores = 1,
        output_format = output_format
    )
    
    check_output_against_phase(
        file.path(outputdir, paste0("stitch.", data_package$chr, extension[output_format])),
        data_package,
        output_format,
        which_snps = NULL,
        tol = 0.2
    )


})


test_that("STITCH can initialize with reference data with three sizes of K vs number of haps", {

    ## these are diploid
    n_samples_per_pop <- 4

    output_format <- "bgvcf"
    L_modified <- 10 + seq(1, 3 * n_snps, 3)
    data_packageL <- make_acceptance_test_data_package(
        n_samples = 10,
        n_snps = n_snps,
        n_reads = n_reads,
        L = L_modified,
        seed = 3,
        chr = original_chr,
        K = 2,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster
    )
    
    for(i_ref in 1:2) {
        
        if (i_ref == 1) {
            ## perfect fit
            refpackL <- make_reference_package(
                L = L_modified,
                n_snps = n_snps,
                n_samples_per_pop = n_samples_per_pop,
                reference_populations = c("CEU"),
                chr = 1,
                phasemaster = phasemaster
            )
        } else if (i_ref == 2) {
            
            ## overlap some
            L_hole <- 10 + seq(1, 3 * (n_snps + 4), 2)
            n_snps_hole <- length(L_hole)
            phasemaster_hole <- matrix(c(rep(0, n_snps_hole), rep(1, n_snps_hole)), ncol = 2)
            phasemaster[c(1, 3, 6), ] <- c(1, 0)
            phasemaster[c(2, 4, 7), ] <- c(0, 0)
            phasemaster[c(3, 5, 8), ] <- c(1, 0)
            t1 <- match(L_hole, L_modified)
            phasemaster_hole[!is.na(t1), ] <- phasemaster[t1[!is.na(t1)], ]
            refpackL <- make_reference_package(
                n_snps = n_snps_hole,
                L = L_hole,
                n_samples_per_pop = n_samples_per_pop,
                reference_populations = c("CEU"),
                chr = original_chr,
                phasemaster = phasemaster_hole
            )
            
        }
        
        for(i_scenario in 1:3) {

            outputdir <- make_unique_tempdir()
            plotHapSumDuringIterations <- FALSE
            ## outputdir <- paste0("~/temp.scenario", i_scenario, "/"); plotHapSumDuringIterations <- TRUE
            
            ## scenarios are: too few, exactly right amount, too many
            K <- list(n_samples_per_pop * 4, n_samples_per_pop * 2, n_samples_per_pop)[[i_scenario]]
            STITCH(
                chr = data_packageL$chr,
                bamlist = data_packageL$bamlist,
                posfile = data_packageL$posfile,
                genfile = data_packageL$genfile,
                outputdir = outputdir,
                reference_haplotype_file = refpackL$reference_haplotype_file,
                reference_legend_file = refpackL$reference_legend_file,
                K = K,
                S = 2,
                nGen = 100,
                nCores = 1,
                output_format = output_format,
                plotHapSumDuringIterations = plotHapSumDuringIterations
            )
            
            check_output_against_phase(
                file.path(outputdir, paste0("stitch.", data_packageL$chr, extension[output_format])),
                data_packageL,
                output_format,
                which_snps = NULL,
                tol = 0.2
            )

        }
    }
    
})


test_that("STITCH can initialize with reference data with defined regionStart and regionEnd", {

    for(output_format in c("bgvcf", "bgen")) {

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

        which_snps <- regionStart:regionEnd
        regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = which_snps,
            tol = 0.2
        )

        load(file.path(outputdir, "RData", paste0("EM.all.",regionName,".RData")))

        ## check hweCount - looks OK!
        ## probably just need this the once?
        g <- data_package$phase[3:4, , 1] + data_package$phase[3:4, , 2]
        expect_equal(rowSums(g == 0), hweCount[, 1])
        expect_equal(rowSums(g == 1), hweCount[, 2])
        expect_equal(rowSums(g == 2), hweCount[, 3])

    }

})

test_that("STITCH can initialize with reference data for certain populations", {

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()

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

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()

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

        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_package,
            output_format,
            which_snps = regionStart:regionEnd,
            tol = 0.2
        )

    }


})



test_that("STITCH can initialize with reference data for certain populations for defined regionStart and regionEnd with non-perfect overlap of SNPs", {

    n_snps <- 200
    L <- 1:200
    regionStart <- 30
    regionEnd <- 170
    buffer <- 7
    
    set.seed(90)

    for(i_test in 1:2) {
        if (i_test == 1) {
            L_data <- sort(sample(L, 120))
            L_haps <- sort(sample(L, 60))
        } else if (i_test == 2) {
            L_data <- sort(sample(L, 60))
            L_haps <- sort(sample(L, 120))
        }
        
        phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
        w <- sample(L, 100)
        phasemaster[w, ] <- 1 - phasemaster[w, ]
        w <- sample(L, 10)
        phasemaster[w, ] <- 0
        w <- sample(L, 10)
        phasemaster[w, ] <- 1
    
        data_packageL <- make_acceptance_test_data_package(
            n_samples = 10,
            n_snps = length(L_data),
            n_reads = length(L_data) * 4,
            seed = 3,
            chr = "chr5",
            K = 2,
            L = L_data,
            reads_span_n_snps = 3,
            phasemaster = phasemaster[L_data, ]
        )
        refpackL <- make_reference_package(
            n_snps = length(L_haps),
            n_samples_per_pop = 4,
            reference_populations = c("CEU", "GBR", "CHB"),
            chr = "chr5",
            L = L_haps,
            phasemaster = phasemaster[L_haps, ]
        )
        
        output_format <- "bgvcf"
        outputdir <- make_unique_tempdir()
        set.seed(621 * i_test)
        
        regionName <- paste0(data_packageL$chr, ".", regionStart, ".", regionEnd)
        
        STITCH(
            chr = data_packageL$chr,
            bamlist = data_packageL$bamlist,
            posfile = data_packageL$posfile,
            genfile = data_packageL$genfile,
            outputdir = outputdir,
            reference_haplotype_file = refpackL$reference_haplotype_file,
            reference_legend_file = refpackL$reference_legend_file,
            reference_sample_file = refpackL$reference_sample_file,
            reference_populations = refpackL$reference_populations[1],
            K = 2,
            nGen = 100,
            nCores = 1,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            output_format = output_format
        )
        
        check_output_against_phase(
            file = file.path(outputdir, paste0("stitch.", regionName, ".vcf.gz")),
            data_package = data_packageL,
            output_format = output_format,
            which_snps = match(
                L_data[(regionStart) <= L_data & L_data <= (regionEnd)],
                L_data
            ),
            tol = 0.2
        )
    }

})


test_that("STITCH can initialize on chromosome X with reference data with defined regionStart and regionEnd", {

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()

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

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()

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

    for(output_format in c("bgvcf", "bgen")) {

        outputdir <- make_unique_tempdir()

        set.seed(795)
        regionStart <- 3
        regionEnd <- 4
        buffer <- 1
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
            output_format = output_format,
            genfile = data_packageX$genfile
        )

        regionName <- paste0(data_packageX$chr, ".", regionStart, ".", regionEnd)
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", regionName, extension[output_format])),
            data_packageX,
            output_format,
            which_snps = regionStart:regionEnd,
            tol = 0.2,
            min_info = 0.90
        )

    }

})


test_that("STITCH can impute with reference panels with only 1 iteration if the initialize with reference data", {

    n_snps <- 20
    set.seed(6)
    n_samples_per_pop <- 4
    Kbuild <- 4
    phasemasterL <- array(0, c(n_snps, Kbuild))
    phasemasterL[, 1] <- rep(c(0, 0, 1, 1), n_snps)[1:n_snps]
    phasemasterL[, 2] <- rep(c(0, 1), n_snps)[1:n_snps]
    phasemasterL[, 3] <- rep(c(1, 0), n_snps)[1:n_snps]
    phasemasterL[, 4] <- rep(c(1, 1, 0, 0), n_snps)[1:n_snps]    
    ## 
    refpackL <- make_reference_package(
        n_snps = n_snps,
        n_samples_per_pop = n_samples_per_pop,
        reference_populations = c("CEU"),
        chr = original_chr,
        phasemaster = phasemasterL
    )
    data_packageL <- make_acceptance_test_data_package(
        reads_span_n_snps = 4,
        n_samples = 10,
        n_snps = n_snps,
        seed = 5,
        chr = original_chr,
        n_reads = n_snps * 2,
        K = ncol(phasemasterL),
        phasemaster = phasemasterL
    )
    output_format <- "bgvcf"
    
    for(i_scenario in 1:3) {

        outputdir <- make_unique_tempdir()
        Krun <- list(n_samples_per_pop * 4, n_samples_per_pop * 2, n_samples_per_pop)[[i_scenario]]        

        STITCH(
            chr = data_packageL$chr,
            bamlist = data_packageL$bamlist,
            posfile = data_packageL$posfile,
            genfile = data_packageL$genfile,            
            outputdir = outputdir,
            reference_haplotype_file = refpackL$reference_haplotype_file,
            reference_legend_file = refpackL$reference_legend_file,
            refillIterations = NA,
            shuffleHaplotypeIterations = NA,
            K = Krun,
            nGen = 100,
            nCores = 1,
            niterations = 1,
            output_format = output_format
        )
        
        check_output_against_phase(
            file.path(outputdir, paste0("stitch.", data_packageL$chr, extension[output_format])),
            data_packageL,
            output_format,
            which_snps = NULL,
            tol = 0.25,
            min_info = 0.85
        )

    }


})

test_that("STITCH errors if niterations=1 with reference panel and posfile is not the same as reference legend", {

    outputdir <- make_unique_tempdir()

    n_snps <- 5
    chr <- 700
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

})


test_that("STITCH can initialize with reference data with snap to grid", {

    for(output_format in c("bgvcf", "bgen")) {

        set.seed(1098)
        outputdir <- make_unique_tempdir()
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
            output_format = output_format,
            keepInterimFiles = TRUE,
            gridWindowSize = 2
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
