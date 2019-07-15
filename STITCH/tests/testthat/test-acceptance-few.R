## edge cases with very few SNPs

n_snps <- 5
n_reads <- 30
reads_span_n_snps <- 2
phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
phasemaster[2, ] <- c(1, 0)
phasemaster[4, ] <- c(1, 0)
L_few <- 6:10
data_packages_few <- lapply(c(FALSE, TRUE), function(samples_are_inbred) {
    make_acceptance_test_data_package(
        n_samples = 20,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 3,
        chr = "chr5",
        K = 2,
        reads_span_n_snps = reads_span_n_snps,
        phasemaster = phasemaster,
        L = L_few,
        samples_are_inbred = samples_are_inbred
    )
})
names(data_packages_few) <- c("outbred", "inbred")
extension <- c("bgvcf" = ".vcf.gz", "bgen" = ".bgen")
data_package_few <- data_packages_few[["outbred"]]


## fail if no SNPs to impute
test_that("STITCH throws an error when no SNPs in region to impute, or if fewer than 2 SNPs to impute", {

    for(i in 1:2) {
        outputdir <- make_unique_tempdir()        
        expect_error(
            STITCH(
                chr = data_package_few$chr,
                bamlist = data_package_few$bamlist,
                posfile = data_package_few$posfile,
                genfile = data_package_few$genfile,
                outputdir = outputdir,
                K = 2,
                nGen = 100,
                nCores = 1,
                regionStart = c(11, 10)[i],
                regionEnd = c(100, 100)[i],
                buffer = c(5, 0)[i]
            ),
            c(
                "There are no SNPs to impute", 
                "There are fewer than 2 SNPs to impute, i.e. there is 1 SNP to impute. In this case, imputation is really just genotyping. STITCH could support genotyping but does not, and note that this kind of defeats the point of imputation. Please use your favourite genotyper e.g. GATK to genotype these SNPs. If you strongly disagree please file a bug report and this can be re-examined"
            )[i]
        )
    }

})

## OK if no SNPs in left buffer, 1 in central, 1 in left buffer
test_that("STITCH works with very few SNPs in central region and buffer", {

    for(output_haplotype_dosages in c(FALSE, TRUE)) {
        for(output_format in c("bgvcf", "bgen")) {

            for(method in setdiff(get_available_methods(), "pseudoHaploid")) {

                if (output_haplotype_dosages) {
                    Ss <- 1
                } else {
                    Ss <- c(3, 1)
                }

                for(S in Ss) {
                    for(i in 1:3) {

                        if (i == 1) {
                            regionStart <- 10; regionEnd <- 100; buffer <- 1 ## 1 1 0
                        } else if (i == 2) {
                            regionStart <- 9; regionEnd <- 100; buffer <- 0 ## 0 2 0
                        } else if (i == 3) {
                            regionStart <- 1; regionEnd <- 6; buffer <- 1 ## 0 1 1
                        }

                        outputdir <- make_unique_tempdir()
                        if (method == "diploid-inbred") {
                            data_package_few <- data_packages_few[["inbred"]]
                        } else {
                            data_package_few <- data_packages_few[["outbred"]]
                        }
                        set.seed(9919)
                        output_filename <- paste0("jimmy", extension[output_format])
                        STITCH(
                            chr = data_package_few$chr,
                            bamlist = data_package_few$bamlist,
                            posfile = data_package_few$posfile,
                            genfile = data_package_few$genfile,        
                            outputdir = outputdir,
                            K = 2,
                            S = S,
                            nGen = 100,
                            nCores = 1,
                            output_format = output_format,
                            regionStart = regionStart,
                            regionEnd = regionEnd,
                            buffer = buffer,
                            output_filename = output_filename,
                            method = method,
                            niterations = 10,
                            shuffleHaplotypeIterations = 2,
                            refillIterations = 4
                        )
                        check_output_against_phase(
                            file = file.path(outputdir, output_filename),
                            data_package = data_package_few,
                            output_format,
                            which_snps = which((regionStart <= L_few) & (L_few<= regionEnd)),
                            tol = 0.2,
                            min_info = 0.80
                        )
                        
                    }

                }

            }
        }
        
    }


})
