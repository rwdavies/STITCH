n_snps <- 10
reads_span_n_snps <- 6
chr <- 1
n_reads <- 5 / (reads_span_n_snps / n_snps) ## want about 5X / sample
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")


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
refpack <- make_reference_package(
    n_snps = n_snps,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    chr = chr,
    phasemaster = phasemaster
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
    chr = data_packageX$chr,
    phasemaster = phasemasterX
)



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
