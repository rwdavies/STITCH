#!/usr/bin/env Rscript


stitch_dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
        stitch_dir <- getwd()
    }
}


library("testthat"); library("STITCH"); library("rrbgen")
dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/"
#dir <- "~/proj/STITCH/"
setwd(paste0(dir, "/STITCH/R"))
a <- dir(pattern = "*.R")
b <- grep("~", a)
if (length(b) > 0) {
a <- a[-b]
}
o <- sapply(a, source)
setwd(dir)
Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

setwd(stitch_dir)


# General variables - modify as appropriate
if (Sys.info()["sysname"] == "Darwin") {
operating_system <- "mac"
} else if(Sys.info()["sysname"] == "Linux"){
operating_system <- "linux"
} else {
stop("Unable to determine OS - please find this error and manually choose OS to properly download data")
}
tempdir <- tempdir() # try /dev/shm/ or put on local fast disk if possible
n_cores <- detectCores() # change as appropriate
inputBundleBlockSize <- NA
# Set directories for testing
# IMPORTANT NOTE - make these new directories, as they are over-written before being written to
mouse_datadir <- file.path(getwd(), "test-data/mouse_data/")
dir.create(mouse_datadir)
mouse_resultsdir <- file.path(getwd(), "test-results/mouse_tests/")
human_datadir <- file.path(getwd(), "test-data/human_data/")
human_matched_to_reference_datadir <- file.path(getwd(), "test-data/human_data_matched_to_reference/")
human_resultsdir <- file.path(getwd(), "test-results/human_tests/")
dir.create(mouse_datadir)
dir.create(mouse_resultsdir)
dir.create(human_datadir)
dir.create(human_resultsdir)
dir.create(human_matched_to_reference_datadir)


# Mouse variables - does not need modifying
mouse_bamlist <- paste0(mouse_datadir, "bamlist.txt")
mouse_genfile <- paste0(mouse_datadir, "gen.txt")
mouse_posfile <- paste0(mouse_datadir, "pos.txt")
mouse_K <- 4
mouse_nGen <- 100
mouse_chr <- "chr19"

# Human variables - does not need modifying
human_genfile <- paste0(human_datadir, "gen_sequencing.txt")
human_posfile <- paste0(human_datadir, "pos.txt")
human_K <- 10
human_nGen <- 4 * 20000 / human_K
human_reference_sample_file <- paste0(human_datadir, "1000GP_Phase3.sample")
human_reference_legend_file <- paste0(human_datadir, "1000GP_Phase3_chr20.legend.gz")
human_reference_haplotype_file <- paste0(human_datadir, "1000GP_Phase3_chr20.hap.gz")


ancillary_http <- "https://mus.well.ox.ac.uk/rwdavies/STITCH/"
ancillary_http <- "https://www.well.ox.ac.uk/~rwdavies/ancillary/"

get_and_untar_if_tgz_file <- function(file) {
if (file.exists(file) == FALSE) {
if (operating_system == "mac") {
system(paste0('curl ', ancillary_http, file, ' -o "', file, '"'))
} else {
system(paste0("wget ", ancillary_http, file))
}
}
if (substr(file, nchar(file) - 3, nchar(file)) == ".tgz")
system(paste0("tar -xzvf ", file))
}

setwd(human_datadir)
##get_and_untar_if_tgz_file("STITCH_human_example_2016_10_18.tgz")

## Get example human reference panel data
## Note - downloaded and extracted from 1000GP_Phase3.tgz at https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html, October 14, 2016
files <- c("1000GP_Phase3.sample", "1000GP_Phase3_chr20.hap.gz", "1000GP_Phase3_chr20.legend.gz")
for(file in files) {
get_and_untar_if_tgz_file(file)
}

## Get example of human data with reference data at exactly the same sites
setwd(human_matched_to_reference_datadir)
##get_and_untar_if_tgz_file("STITCH_human_reference_example_2018_07_11.tgz")











##
##

library("proftools")
library("STITCH")

## change directory to one up from scripts, no matter how this was called
profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)

profile_start <- Sys.time()

###############################################
outputdir <- paste0(human_resultsdir, "with_ref_panel/")
unlink(outputdir, recursive = TRUE)
system(paste0("rsync -a ", shQuote(human_datadir), "/* ", shQuote(outputdir)))
outputdir <- paste0(human_resultsdir, "with_ref_panel/")
dir.create("/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/test-results/human_tests/with_ref_panel//plots/", recursive = TRUE)
##
## startup
##
chr = "20"
out <- get_and_validate_pos_gen_and_phase(
    posfile = human_posfile,
    chr = chr,
    verbose = TRUE
)
pos <- out$pos
nSNPs <- out$nSNPs
L <- pos[, 2]
regionName <- chr
alleleCount <- array(0, c(nSNPs, 3))
gridWindowSize <- NA
cM <- NULL
##
## set up grid
##
out <- assign_positions_to_grid(
    L = L,
    gridWindowSize = gridWindowSize,
    cM = cM
)
grid <- out$grid
grid_distances <- out$grid_distances
L_grid <- out$L_grid
nGrids <- out$nGrids
snps_in_grid_1_based <- out$snps_in_grid_1_based
cM_grid <- out$cM_grid

##
    reference_populations = c("CHB", "CHS", "CHD")
    reference_haplotype_file = human_reference_haplotype_file
    reference_sample_file = human_reference_sample_file
    reference_legend_file = human_reference_legend_file
    reference_phred = 20
    reference_iterations = 40
    nSNPs = nSNPs
    K = 100
    S = 2
    L  = L
    pos = pos
    inputBundleBlockSize = NA
    nCores = 1
    regionName = regionName
    alleleCount = alleleCount
    expRate = 0.5
    nGen = 4 * 20000 / 100
    tempdir = tempdir()
    outputdir = outputdir
    pseudoHaploidModel = 9
    alphaMatThreshold = 1e-4
    emissionThreshold = 1e-4
    maxRate = 100
    minRate = 0.1
    regionStart = 1000000
    regionEnd = 1100000
    buffer = 10000
    niterations = 40
    grid = grid
    grid_distances = grid_distances
    nGrids = nGrids
    reference_shuffleHaplotypeIterations = c(4, 8, 12, 16)
    L_grid = L_grid
    plot_shuffle_haplotype_attempts =TRUE
    shuffle_bin_radius = 5000
    snps_in_grid_1_based = snps_in_grid_1_based
    plotHapSumDuringIterations = FALSE
    cM_grid = cM_grid

##
## the thing I care about!
##
out <- initialize_parameters(
    reference_populations = c("CHB", "CHS", "CHD"),
    reference_haplotype_file = human_reference_haplotype_file,
    reference_sample_file = human_reference_sample_file,
    reference_legend_file = human_reference_legend_file,
    reference_phred = 20,
    reference_iterations = 40,
    nSNPs = nSNPs,
    K = 100,
    S = 2,
    L  = L,
    pos = pos,
    inputBundleBlockSize = NA,
    nCores = 1,
    regionName = regionName,
    alleleCount = alleleCount,
    expRate = 0.5,
    nGen = 4 * 20000 / 100,
    tempdir = tempdir(),
    outputdir = outputdir,
    pseudoHaploidModel = 9,
    alphaMatThreshold = 1e-4,
    emissionThreshold = 1e-4,
    maxRate = 100,
    minRate = 0.1,
    regionStart = 1000000,
    regionEnd = 1100000,
    buffer = 10000,
    niterations = 40,
    grid = grid,
    grid_distances = grid_distances,
    nGrids = nGrids,
    reference_shuffleHaplotypeIterations = c(4, 8, 12, 16),
    L_grid = L_grid,
    plot_shuffle_haplotype_attempts =TRUE,
    shuffle_bin_radius = 5000,
    snps_in_grid_1_based = snps_in_grid_1_based,
    plotHapSumDuringIterations = FALSE,
    cM_grid = cM_grid
)
##################
profile_end <- Sys.time()

Rprof(NULL)

pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- Sys.getenv("OUTPUT_PLOT")
if (output_plot == "") {
setwd(stitch_dir)
output_plot <- file.path("profile.pdf")
}
pdf(output_plot, height = 48, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
title(title, outer=TRUE)
dev.off()
system("rsync -av profile.pdf gru:~/")
print(profile_end - profile_start)



quit()




STITCH(
       outputdir = outputdir,
		 method = "diploid",
		 originalRegionName = "20.1000000.1100000",
		 regenerateInput = FALSE,
		 regionStart = 1000000,
		 regionEnd = 1100000,
		 buffer = 10000,
		 niterations = 1,
		 chr = "20",
		 inputBundleBlockSize = 100,
		 reference_populations = c("CHB", "CHS", "CHD"),
		 reference_haplotype_file = human_reference_haplotype_file,
		 reference_sample_file = human_reference_sample_file,
		 reference_legend_file = human_reference_legend_file,
		 shuffleHaplotypeIterations = NA,
		 refillIterations = NA,
		 keepInterimFiles = TRUE,
		 genfile = human_genfile,
		 posfile = human_posfile,
		 K = human_K,
		 tempdir = tempdir,
		 nCores = 1,
		 nGen = human_nGen
		 )
