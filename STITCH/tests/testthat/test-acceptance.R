run_acceptance_tests <- FALSE

n_snps <- 5
chr <- 10
phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
data_package <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = 4,
    seed = 1,
    chr = chr,
    K = 2,
    phasemaster = phasemaster
)
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = chr
)

data_package_crams <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = 4,
    seed = 1,
    chr = chr,
    K = 2,
    phasemaster = phasemaster,
    use_crams = TRUE
)


chr <- "X"
phasemasterX <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
data_packageX <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = 4,
    seed = 1,
    chr = chr,
    K = 2,
    phasemaster = phasemasterX
)
refpackX <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = data_packageX$chr
)



test_that("STITCH diploid works under default parameters", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        K = 2,
        nGen = 100,
        nCores = 1
    )

    sink()

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH diploid works under default parameters with nCores = 40 and N = 25", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    data_package25 <- make_acceptance_test_data_package(
        n_samples = 25,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = 2,
        phasemaster = phasemaster
    )

    STITCH(
        tempdir = tempdir(),
        chr = data_package25$chr,
        bamlist = data_package25$bamlist,
        posfile = data_package25$posfile,
        outputdir = data_package25$outputdir,
        K = 4,
        nGen = 100,
        nCores = 40
    )

    sink()

    vcf <- read.table(
        file.path(data_package25$outputdir, paste0("stitch.", data_package25$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package25$phase,
        tol = 0.2
    )

})


test_that("STITCH diploid works under default parameters with N = 25 and nCores = 40", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    data_package25 <- make_acceptance_test_data_package(
        n_samples = 25,
        n_snps = n_snps,
        n_reads = 4,
        seed = 1,
        chr = chr,
        K = 2,
        phasemaster = phasemaster
    )

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package25$chr,
        bamlist = data_package25$bamlist,
        posfile = data_package25$posfile,
        outputdir = data_package25$outputdir,
        K = 2,
        nGen = 100,
        nCores = 40
    )

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH pseudoHaploid works under default parameters", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    set.seed(11)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        method = "pseudoHaploid"
    )

    sink()

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH pseudoHaploid works with switchModelIteration", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    set.seed(12)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        method = "pseudoHaploid",
        switchModelIteration = 39,
        niterations = 40
    )

    sink()

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH pseudoHaploid works with a single sample and two cores", {
    skip_test_if_TRUE(run_acceptance_tests)

    phasemaster <- matrix(c(c(0, 0, 0), c(1, 1, 1)), ncol = 2)
    data_package <- make_acceptance_test_data_package(
        n_samples = 1,
        n_snps = 3,
        n_reads = 4,
        seed = 1,
        chr = 10,
        K = 2,
        phasemaster = phasemaster
    )

    sink("/dev/null")

    set.seed(11)

    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        K = 2,
        nGen = 100,
        nCores = 2,
        method = "pseudoHaploid"
    )

    sink()

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH can initialize with reference data", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        K = 2,
        nGen = 100,
        nCores = 1
    )
    sink()


    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH can initialize with reference data with defined regionStart and regionEnd", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    set.seed(10)

    regionStart <- 3
    regionEnd <- 4
    buffer <- 1
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        K = 2,
        nGen = 100,
        nCores = 1,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    sink()

    vcf <- read.table(
        file.path(
            data_package$outputdir,
            paste0("stitch.", data_package$chr, ".", regionStart, ".", regionEnd, ".vcf.gz")
        ),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    keep <- data_package$pos[, "POS"] >= regionStart & data_package$pos[, "POS"] <= regionEnd
    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase[keep, , ],
        tol = 0.2
    )

})

test_that("STITCH can initialize with reference data for certain populations", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations[1],
        K = 2,
        nGen = 100,
        nCores = 1
    )

    sink()

    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})


test_that("STITCH can initialize with reference data for certain populations for defined regionStart and regionEnd", {
    skip_test_if_TRUE(run_acceptance_tests)
    sink("/dev/null")

    regionStart <- 3
    regionEnd <- 4
    buffer <- 1
    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        reference_sample_file = refpack$reference_sample_file,
        reference_populations = refpack$reference_populations[1],
        K = 2,
        nGen = 100,
        nCores = 1,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    sink()

    vcf <- read.table(
        file.path(
            data_package$outputdir,
            paste0("stitch.", data_package$chr, ".", regionStart, ".", regionEnd, ".vcf.gz")
        ),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    keep <- data_package$pos[, "POS"] >= regionStart & data_package$pos[, "POS"] <= regionEnd
    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase[keep, ,],
        tol = 0.2
    )

})



test_that("STITCH can initialize on chromosome X with reference data", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_packageX$chr,
        bamlist = data_packageX$bamlist,
        posfile = data_packageX$posfile,
        outputdir = data_packageX$outputdir,
        reference_haplotype_file = refpackX$reference_haplotype_file,
        reference_legend_file = refpackX$reference_legend_file,
        K = 2,
        nGen = 100,
        nCores = 1
    )

    sink()

    vcf <- read.table(
        file.path(data_packageX$outputdir, paste0("stitch.", data_packageX$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_packageX$phase,
        tol = 0.2
    )

})


test_that("STITCH can initialize on chromosome X with reference data with defined regionStart and regionEnd", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    regionStart <- 3
    regionEnd <- 4
    buffer <- 1
    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_packageX$chr,
        bamlist = data_packageX$bamlist,
        posfile = data_packageX$posfile,
        outputdir = data_packageX$outputdir,
        reference_haplotype_file = refpackX$reference_haplotype_file,
        reference_legend_file = refpackX$reference_legend_file,
        K = 2,
        nGen = 100,
        nCores = 1,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    sink()

    vcf <- read.table(
        file.path(
            data_packageX$outputdir,
            paste0("stitch.", data_packageX$chr, ".", regionStart, ".", regionEnd, ".vcf.gz")
        ),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    keep <- data_packageX$pos[, "POS"] >= regionStart & data_packageX$pos[, "POS"] <= regionEnd
    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase[keep, , ],
        tol = 0.2
    )

})


test_that("STITCH can initialize on chromosome X with reference data for certain populations", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_packageX$chr,
        bamlist = data_packageX$bamlist,
        posfile = data_packageX$posfile,
        outputdir = data_packageX$outputdir,
        reference_haplotype_file = refpackX$reference_haplotype_file,
        reference_legend_file = refpackX$reference_legend_file,
        reference_sample_file = refpackX$reference_sample_file,
        reference_populations = refpackX$reference_populations[1],
        K = 2,
        nGen = 100,
        nCores = 1
    )

    sink()

    vcf <- read.table(
        file.path(data_packageX$outputdir, paste0("stitch.", data_packageX$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_packageX$phase,
        tol = 0.2
    )

})




test_that("STITCH can initialize on chromosome X with reference data for certain populations for defined regionStart and regionEnd", {
    skip_test_if_TRUE(run_acceptance_tests)

    regionStart <- 3
    regionEnd <- 4
    buffer <- 1

    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_packageX$chr,
        bamlist = data_packageX$bamlist,
        posfile = data_packageX$posfile,
        outputdir = data_packageX$outputdir,
        reference_haplotype_file = refpackX$reference_haplotype_file,
        reference_legend_file = refpackX$reference_legend_file,
        reference_sample_file = refpackX$reference_sample_file,
        reference_populations = refpackX$reference_populations[1],
        K = 2,
        nGen = 100,
        nCores = 1,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer
    )

    sink()

    vcf <- read.table(
        file.path(
            data_packageX$outputdir,
            paste0("stitch.", data_packageX$chr, ".", regionStart, ".", regionEnd, ".vcf.gz")
        ),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    keep <- data_packageX$pos[, "POS"] >= regionStart & data_packageX$pos[, "POS"] <= regionEnd
    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase[keep, , ],
        tol = 0.2
    )

})


test_that("STITCH diploid works under default parameters using CRAM files", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    STITCH(
        tempdir = tempdir(),
        chr = data_package_crams$chr,
        cramlist = data_package_crams$cramlist,
        reference = data_package_crams$ref,
        posfile = data_package_crams$posfile,
        outputdir = data_package_crams$outputdir,
        K = 2,
        nGen = 100,
        nCores = 1
    )

    vcf <- read.table(
        file.path(data_package_crams$outputdir, paste0("stitch.", data_package_crams$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.2
    )

})





test_that("STITCH with generateInputOnly actually only generates input", {
    skip_test_if_TRUE(run_acceptance_tests)

    sink("/dev/null")

    set.seed(10)
    local_outputdir <- tempfile()
    dir.create(local_outputdir)
    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = local_outputdir,
        K = 2,
        nGen = 100,
        nCores = 1,
        generateInputOnly = TRUE
    )

    expect_equal(
        FALSE,
        file.exists(file.path(local_outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")))
    )
    inputdir_contents <- dir(file.path(local_outputdir, "input"))
    expect_equal(
        sort(inputdir_contents),
        sort(paste0("sample.", 1:10, ".input.", data_package$chr, ".RData"))
    )


})


test_that("STITCH can impute with reference panels with only 1 iteration if the initialize with reference data", {
    skip_test_if_TRUE(run_acceptance_tests)
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

    STITCH(
        tempdir = tempdir(),
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = data_package$outputdir,
        reference_haplotype_file = refpack$reference_haplotype_file,
        reference_legend_file = refpack$reference_legend_file,
        refillIterations = NA,
        shuffleHaplotypeIterations = NA,
        K = 4,
        nGen = 100,
        nCores = 1,
        niterations = 1
    )
    sink()


    vcf <- read.table(
        file.path(data_package$outputdir, paste0("stitch.", data_package$chr, ".vcf.gz")),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    check_vcf_against_phase(
        vcf = vcf,
        phase = data_package$phase,
        tol = 0.3 ## meh, a bit looser
    )

})

test_that("STITCH errors if niterations=1 with reference panel and posfile is not the same as reference legend", {
    skip_test_if_TRUE(run_acceptance_tests)
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
            tempdir = tempdir(),
            chr = data_package$chr,
            bamlist = data_package$bamlist,
            posfile = data_package$posfile,
            outputdir = data_package$outputdir,
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
