#!/usr/bin/env Rscript

## How to use this file:
## This file contains fully working examples of how to use STITCH, both to serve as a guide on the functionality of STITCH, and to verify that supported functionality works on a target machine. You should be able to run this script through running ./scripts/all-tests.sh on your machine (although consider changing the nCores variable below as appropriate).
## Alternatively, to explore this script interactively, first, ensure you can run STITCH in R using the guide given in the README. Then, you can copy and paste Sections A-C into an R session, and then run specific examples from Sections D and E as appropriate.
## Note that only Linux and Mac OS are supported (not Windows).


# This file is broken into the following high level components
# Section A - Set common variables / defaults for mice, humans or both
# Section B - List the examples found in the rest of the script
# Section C - Download data
# Section D - Mouse examples
# Section E - Human examples


## change directory to one up from examples, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

library("STITCH")



###
### Section A - Set common variables
###


# General variables - modify as appropriate
if (Sys.info()["sysname"] == "Darwin") {
  operating_system <- "mac"
} else if(Sys.info()["sysname"] == "Linux"){
    operating_system <- "linux"
} else {
  stop("Unable to determine OS - please find this error and manually choose OS to properly download data")
}
tempdir <- tempdir() # try /dev/shm/ or put on local fast disk if possible
server_environment <- "server"
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

###
### Section B - List of examples
###


# ------ Mouse examples
# Mouse example 1 - Run whole chromosome
# Mouse example 2 - Run part of a chromosome (e.g. a few megabases, with buffer)
# Mouse example 3 - Use fewer internal input files (useful for very large N, e.g. > 5000)
# Mouse example 4 - Pseudo-haploid method
# Mouse example 5 - 2-step part 1 - Generate input for a region (useful when BAMs are not accesible on a large computational cluster)
# Mouse example 6 - 2-step part 2 - Run using previously generated input
# Mouse example 7 - Run 2-step imputation on a smaller region
# Mouse example 8 - Output original input in VCF format (useful to run other programs on the input)
# Mouse example 9 - Downsample the number of samples (useful for evaluating how performance varies with number of samples)
# Mouse example 10 - Downsample the read count of samples (useful for evaluating how performance varies with read depth)
# Mouse example 11 - Run on CRAM files
# Mouse example 12 - Change VCF output name
# Mouse example 13 - Use gridWindowSize to speed up analysis
# Mouse example 14 - Write to bgen

# ------ Human examples
# Human example 1 - Run diploid on previously generated input, downsample samples and coverage
# Human example 2 - Run with reference panel









###
### Section C - Download data for mouse and human examples
###


ancillary_http <- "http://mus.well.ox.ac.uk/rwdavies/STITCH/"
ancillary_http <- "http://www.well.ox.ac.uk/~rwdavies/ancillary/"

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

## Download and prepare mouse data
setwd(mouse_datadir)
get_and_untar_if_tgz_file("STITCH_example_2016_05_10.tgz")

## Download mouse cram data
get_and_untar_if_tgz_file("STITCH_mouse_crams_2017_02_05.tgz")


## Download and prepare human data
setwd(human_datadir)
get_and_untar_if_tgz_file("STITCH_human_example_2016_10_18.tgz")

## Get example human reference panel data
## Note - downloaded and extracted from 1000GP_Phase3.tgz at https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html, October 14, 2016
files <- c("1000GP_Phase3.sample", "1000GP_Phase3_chr20.hap.gz", "1000GP_Phase3_chr20.legend.gz")
for(file in files) {
    get_and_untar_if_tgz_file(file)
}

## Get example of human data with reference data at exactly the same sites
setwd(human_matched_to_reference_datadir)
get_and_untar_if_tgz_file("STITCH_human_reference_example_2017_05_24.tgz")


## Get mouse reference genome - need the reference genome to work with CRAM files
setwd(mouse_datadir)
get_and_untar_if_tgz_file("mm10_2016_10_02.fa.gz")
if (file.exists("mm10_2016_10_02.fa") == FALSE | file.exists("mm10_2016_10_02.fa.fai") == FALSE) {
    system("gzip -cd mm10_2016_10_02.fa.gz > mm10_2016_10_02.fa")
    system(paste0("samtools faidx mm10_2016_10_02.fa"))
}


###
### Section D - Mouse Examples
###







# Mouse example 1 - Run whole chromosome
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "whole_region/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)



# Mouse example 2 - Run part of a chromosome
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "subset/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)



# Mouse example 3 - Run part of a chromosome, internally, use fewer input files (useful when N > 5000)
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "bundle/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  inputBundleBlockSize = 100,
  method = "diploid",
  outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen)



# Mouse example 4 - Pseudo-haploid method
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "pseudo_haploid/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  inputBundleBlockSize = 100,
  method = "pseudoHaploid",
  outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen)



# Mouse example 5 - 2-step part 1 - Generate input for a region
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "two_stage/")
outputdir_example5 <- outputdir
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  generateInputOnly = TRUE,
  outputdir = outputdir, method = "diploid", chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)



# Mouse example 6 - 2-step part 2 - Run using previously generated input (useful when BAMs are not accesible from server
setwd(mouse_datadir)
outputdir_moved <- paste0(mouse_resultsdir, "two_stage_moved/")
system(paste0("rm -r ", outputdir_moved), ignore.stderr = TRUE)
system(paste0("mkdir -p ", outputdir_moved))
system(paste0("rsync -a ", outputdir_example5, " ", outputdir_moved))
STITCH(
  outputdir = outputdir_moved,
  method = "diploid",
  originalRegionName = "chr19.10000000.10800000",
  regenerateInput = FALSE,
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)



# Mouse example 7 - Run 2-step imputation on a smaller region
setwd(mouse_datadir)
outputdir_moved <- paste0(mouse_resultsdir, "two_stage_moved_smaller/")
system(paste0("rm -r ", outputdir_moved), ignore.stderr = TRUE)
system(paste0("mkdir -p ", outputdir_moved))
system(paste0("rsync -a ", outputdir_example5, " ", outputdir_moved))
STITCH(
  outputdir = outputdir_moved,
  method = "diploid",
  originalRegionName = "chr19.10000000.10800000",
  regenerateInput = FALSE,
  regionStart = 10300000,
  regionEnd = 10700000,
  buffer = 50000,
  chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)




# Mouse example 8 - Output original input in VCF format (useful to run other programs on the input)
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "output_input_as_vcf/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  outputInputInVCFFormat = TRUE,
  method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)



# Mouse example 9 - Downsample the number of samples (useful for evaluating how performance varies with number of samples)
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "downsample_samples/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  downsampleSamples = 0.25,
  method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen)



# Mouse example 10 - Downsample the read count of samples (useful for evaluating how performance varies with read depth)
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "downsample_reads/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  downsampleFraction = 0.5,
  method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen)



# Mouse example 11 - Run on CRAM files
setwd(mouse_datadir)
### Either un-comment below and convert
### Or use downloaded mouse CRAM files
#system("mkdir crams")
#bams <- as.character(read.table(mouse_bamlist)[, 1])
#print("Note that conversion of BAM to CRAM may be fairly slow")
#crams <- unlist(mclapply(bams, mc.cores = n_cores, function(bam) {
#  i_bam <- match(bam, bams)
#  if ((i_bam %% 100) == 1)
#    print(paste0("Done CRAM conversion of ", i_bam, " / ", length(bams), ", ", date()))
#  cram <- gsub("bam", "cram", gsub("bams", "crams", bam))
#  if (file.exists(cram) == FALSE) {
#    system(paste0("samtools view -T ", mouse_datadir, "mm10_2016_10_02.fa -C -o ", cram, " ", bam))
#    system(paste0("samtools index ", cram))
#  }
#  return(cram)
#}))

cramlist <- paste0(mouse_datadir, "cramlist.txt")
#write.table(matrix(crams, ncol = 1), file = cramlist, row.names = FALSE, col.names = FALSE, quote = FALSE)

setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "use_crams/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  method = "diploid",
  cramlist = cramlist,
  bamlist = "",
  reference = "mm10_2016_10_02.fa",
  chr = mouse_chr,
  K = mouse_K, outputdir = outputdir, posfile = mouse_posfile, genfile = mouse_genfile, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)

# Mouse example 12 - Change VCF output name
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "change_name/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(output_filename = "test.vcf.gz", method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize)

# Mouse example 13 - Use gridWindowSize to speed up analysis
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "gridWindowSize/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(
  regionStart = 10000000,
  regionEnd = 10800000,
  buffer = 50000,
  gridWindowSize = 10000,
  inputBundleBlockSize = 100,  
  method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen)

# Mouse example 14 - Write to bgen
setwd(mouse_datadir)
outputdir <- paste0(mouse_resultsdir, "bgen/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
STITCH(method = "diploid", outputdir = outputdir, chr = mouse_chr, posfile = mouse_posfile, genfile = mouse_genfile, bamlist = mouse_bamlist, K = mouse_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = mouse_nGen, inputBundleBlockSize = inputBundleBlockSize, output_format = "bgen")




###
### Section E - Human examples
###



# Human example 1 - Run diploid on previously generated input, downsample samples and coverage
outputdir <- paste0(human_resultsdir, "downsampling/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
system(paste0("rsync -a ", human_datadir, "/* ", outputdir))
STITCH(
  outputdir = outputdir,
  method = "diploid",
  originalRegionName = "20.1000000.1100000",
  regenerateInput = FALSE,
  regionStart = 1000000,
  regionEnd = 1100000,
  buffer = 10000,
  niterations = 20,
  downsampleSamples = 0.1,
  chr = "20",
  shuffleHaplotypeIterations = NA,
  refillIterations = NA,
  inputBundleBlockSize = 100,
  genfile = human_genfile, posfile = human_posfile, K = human_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = human_nGen)


# Human example 2 - Run with reference panel
outputdir <- paste0(human_resultsdir, "with_ref_panel/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
system(paste0("rsync -a ", human_datadir, "/* ", outputdir))
STITCH(
  outputdir = outputdir,
  method = "diploid",
  originalRegionName = "20.1000000.1100000",
  regenerateInput = FALSE,
  regionStart = 1000000,
  regionEnd = 1100000,
  buffer = 10000,
  niterations = 5,
  downsampleSamples = 1,
  downsampleFraction = 0.1,
  chr = "20",
  inputBundleBlockSize = 100,
  reference_populations = c("CHB", "CHS", "CHD"),
  reference_haplotype_file = human_reference_haplotype_file,
  reference_sample_file = human_reference_sample_file,
  reference_legend_file = human_reference_legend_file,
  shuffleHaplotypeIterations = NA,
  refillIterations = NA,
  genfile = human_genfile, posfile = human_posfile, K = human_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = human_nGen)


# Human example 3 - Run with reference panel with no updating
outputdir <- paste0(human_resultsdir, "matched_to_ref_panel/")
system(paste0("rm -r ", outputdir), ignore.stderr = TRUE)
system(paste0("rsync -a ", human_matched_to_reference_datadir, "/* ", outputdir))
human_matched_to_reference_genfile <- paste0(human_matched_to_reference_datadir, "gen_sequencing.intersect.txt")
human_matched_to_reference_posfile <- paste0(human_matched_to_reference_datadir, "pos.intersect.txt")
human_matched_to_reference_reference_legend_file <- paste0(human_matched_to_reference_datadir, "1000GP_Phase3_20.1000000.1100000.legend.gz")
human_matched_to_reference_reference_haplotype_file <- paste0(human_matched_to_reference_datadir, "1000GP_Phase3_20.1000000.1100000.hap.gz")
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
  reference_haplotype_file = human_matched_to_reference_reference_haplotype_file,
  reference_sample_file = human_reference_sample_file,
  reference_legend_file = human_matched_to_reference_reference_legend_file,
  shuffleHaplotypeIterations = NA,
  refillIterations = NA,
  genfile = human_matched_to_reference_genfile, posfile = human_matched_to_reference_posfile, K = human_K, tempdir = tempdir, environment = server_environment, nCores = n_cores, nGen = human_nGen)
