#!/usr/bin/env Rscript

message("--- Running script to compare VCF and truth ---")
## Script to compare VCF output to bespoke array genotypes
## Could be generalized

args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

##setwd("~/proj/STITCH")
source("scripts/compare_vcf_to_truth_functions.R")
source("STITCH/R/functions.R") ## for function generate_hwe_on_counts

library("data.table")
library("optparse")
library("rrbgen")

option_list <- list(
    make_option(
        "--test-file",
        type = "character",
        help = "VCF output",
        default = NA
    ),
    make_option(
        "--truth-vcf",
        type = "character",
        help = "VCF output",
        default = NA
    ),
    make_option(
        "--chr",
        type = "character",
        help = "chr",
        default = NA
    ),
    make_option(
        "--compare-against",
        type = "character",
        help = "Compare against either affy, megamuga, or truth_vcf",
        dest = "compare_against",
        default = "nothing"
    ),
    make_option(
        "--verbose",
        action = "store_true",
        default = FALSE
    ),
    make_option(
        "--rebuild-mega-save-file",
        action = "store_true",
        default = FALSE,
        dest = "rebuild_mega_save_file"
    ),
    make_option(
        "--mega-save-file",
        type = "character",        
        help = "RData file to save mega-muga processed results",
        dest = "mega_save_file",
        default = NA
    ),
    make_option(
        "--comparison-save-file",
        type = "character",        
        help = "RData file to save final compared results to",
        dest = "comparison_save_file",
        default = NA
    ),
    make_option(
        "--cfw-data-dir",
        type = "character",        
        help = "Where CFW truth files are stored",
        dest = "cfw_data_dir",
        default = "/data/smew1/rdavies/stitch_development/truth/cfw/"
    ),
    make_option(
        "--whose-samples",
        type = "character",        
        help = "Whether to analyze oxford, chicago, or both samples",
        dest = "whose_samples",
        default = "oxford"
    )
)

opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))

if (1 == 0) {

    opt <- list(
        "test-file"="/well/myers/rwdavies/chicago_oxford_2019_08_07/stitch.chr19.vcf.gz",
        chr = "chr19",
        compare_against="megamuga",
        verbose=TRUE,
        whose_samples="both",
        rebuild_mega_save_file=TRUE,
        cfw_data_dir="/well/myers/rwdavies/cfw_truth/",
        mega_save_file="/well/myers/rwdavies/chicago_oxford_2019_08_07/mega.chr19.RData"
    )
    ## opt[["test-file"]] <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/test-results//whole_chr_CFW_1.5.0//stitch.chr19.vcf.gz";
    
}

test_file <- opt[["test-file"]]
truth_vcf_file <- opt$truth_vcf
chr <- opt$chr
whose_samples <- opt$whose_samples
compare_against <- opt$compare_against
verbose <- opt$verbose
mega_save_file <- opt$mega_save_file
comparison_save_file <- opt$comparison_save_file
cfw_data_dir <- opt$cfw_data_dir
rebuild_mega_save_file <- opt$rebuild_mega_save_file
megamugadir <- file.path(cfw_data_dir, "megamuga")
affydir <- file.path(cfw_data_dir, "affy")


if (substr(test_file, nchar(test_file) - 3, nchar(test_file)) == "bgen") {
    output_format <- "bgen"
} else {
    output_format <- "bgvcf"
}


if (is.na(chr))
    stop("specify chr")

if (compare_against == "affy") {

    message("--- Compare Affymetrix to Input VCF ---")
    affy_calls <- get_affymetrix_calls_for_chr(chr)
    affy_calls <- filter_calls_for_chr(affy_calls)
    
    subjects <- colnames(affy_calls)[grep("Q_", colnames(affy_calls))]
    truth_calls <- affy_calls
    
} else if (compare_against == "megamuga") {

    message("--- Compare MegaMUGA to Input VCF ---")
    if (
        !rebuild_mega_save_file &
        !is.na(mega_save_file) &
        file.exists(mega_save_file)
    ) {
        message(paste0("Load previous MegaMuga file for ", chr))
        load(mega_save_file)
    } else {
        megamuga <- get_megamuga_calls_for_chr(chr)
        mega_calls <- convert_megamuga_to_calls(megamuga, chr, whose_samples)
        mega_calls <- filter_calls_for_chr(mega_calls)
        subjects <- setdiff(colnames(mega_calls), c("chr", "pos"))
        if (is.na(mega_save_file) == FALSE) {
            message("Save MegaMUGA calls for chr")
            save(mega_calls, subjects, file = mega_save_file)
        }
    }
    truth_calls <- mega_calls

} else if (compare_against == "truth_vcf") {

    stop("This is not written or tested")
    message("--- Compare Truth VCF to Input VCF ---")
    out <- get_dosages_for_vcf(truth_vcf_file, chr, subjects = NULL, calls = NULL, use_fread = FALSE, tabix = FALSE)
    truth_dosages <- out$dosages
    out <- get_dosages_for_vcf(test_file, chr, subjects = NULL, calls = NULL, use_fread = FALSE, tabix = FALSE)
    dosages <- out$dosages

} else {
    
    stop("incorrect choice of comparison")

}


if (output_format == "bgvcf") {
    out <- get_dosages_for_vcf(test_file, chr, subjects, truth_calls)
    ## out$dosages with row = SNP, col = sample, content is 0-2 dosage
    ## rows labelled with chr-snp-ref-alt
    ## out$dosages_meta = 3 cols with hwe, eaf, info
} else {
    out <- get_dosages_for_bgen(test_file, subjects, truth_calls)
}

dosages <- out$dosages
dosages_meta <- out$dosages_meta
genotypes <- out$genotypes

out <- compare_array_calls_and_dosages(
    calls = truth_calls,
    dosages = dosages,
    dosages_meta = dosages_meta,
    genotypes = genotypes,
    return_things = !is.na(comparison_save_file)
)
if (!is.na(comparison_save_file)) {
    message(paste0("Saving results to:", comparison_save_file))
    save(out, file = comparison_save_file)
}


quit()
