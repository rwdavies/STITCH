if (1 == 0) {
    
    library("testthat"); library("STITCH"); library("rrbgen")
    ##    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    dir <- "~/proj/STITCH/"
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


n_snps <- 5
n_reads <- 20
reads_span_n_snps <- 2
phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
phasemaster[2, ] <- c(1, 0)
phasemaster[4, ] <- c(1, 0)
L_few <- 6:10
data_package_few <- make_acceptance_test_data_package(
    n_samples = 10,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 3,
    chr = "chr5",
    K = 2,
    reads_span_n_snps = reads_span_n_snps,
    phasemaster = phasemaster,
    L = L_few
)


## OK if no SNPs in left buffer, 1 in central, 1 in left buffer
test_that("STITCH works with very few SNPs in central region and buffer", {

    for(output_haplotype_dosages in c(FALSE, TRUE)) {
        for(output_format in c("bgvcf", "bgen")) {

            ## for(method in get_available_methods()) {
            for(method in "pseudoHaploid") {

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
                            file =  file.path(outputdir, output_filename),
                            data_package = data_package_few,
                            output_format,
                            which_snps = which((regionStart <= L_few) & (L_few<= regionEnd)),
                            tol = 0.2,
                            min_info = 0.95
                        )
                        
                    }

                }

            }
        }
        
    }


})
