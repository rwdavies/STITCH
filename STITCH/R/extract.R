#' Extract ancestral haplotype dosage to RData cube
#' @param vcf_file path to VCF
#' @param ref path to reference fasta (required by GATK)
#' @param bcftools path to vcftools (or just "bcftools" if in path)
#' @param gatk_jar path to GATK jar
#' @param samples vector of sample names (or NULL for all)
#' @param pos pos matrix, a matrix with at least two columns where the first two columns are chrom and 1-based physical position, respectively. specifies which SNPs to extract (or NULL for all)
#' @param field What to get from the VCF. Default HD for haplotype dosages
#' @return A cube with dimensions of SNPs x samples x ancestral haplotypes
#' @author Robert Davies
#' @export
extract_hd_to_cube <- function(
    vcf_file,
    ref,
    bcftools = "bcftools",
    gatk_jar = "/data/smew1/rdavies/stitch_richard_paper/bin/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar",
    samples = NULL,
    pos = NULL,
    ram = "-Xmx4g",
    field = "HD",
    verbose = TRUE
) {
    ##
    ## bcftools
    ##
    command <- paste0(bcftools, " view ", vcf_file)
    if (!is.null(samples)) {
        samples_file <- tempfile()
        write.table(matrix(samples, ncol = 1), quote = FALSE, col.names =FALSE, row.names = FALSE, file = samples_file)
        command <- paste0(command, " -S ", samples_file)
    }
    if (!is.null(pos)) {
        snps_file <- tempfile()
        write.table(pos[, 1:2], sep = "\t", quote = FALSE, col.names =FALSE, row.names = FALSE, file = snps_file)
        command <- paste0(command, " -R ", snps_file)
    }
    vcf_output_file <- tempfile(fileext = ".vcf")
    if (verbose) {
        message("Extract using bcftools")
    }
    command <- paste0(command, " > ", vcf_output_file)
    system(command)
    system(paste0("bgzip -f ", vcf_output_file))
    system(paste0("tabix -f ", vcf_output_file, ".gz"))
    ##
    ## GATK
    ##
    table_file <- tempfile(fileext = ".table")
    args <- c(
        ram, "-jar", gatk_jar, "-T", "VariantsToTable", "-R", ref,
        "-V", paste0(vcf_output_file, ".gz"),
        "--genotypeFields", field,
        "-o", table_file
    )
    if (verbose) {
        message("Convert using GATK")
    }
    stdout <- tempfile()
    stderr <- tempfile()
    result <- system2("java", args, stdout = stdout, stderr = stderr)
    if (result > 0) {
        system(paste0("cat ", stdout))
        system(paste0("cat ", stderr))
        stop("Failed to extract properly")
    }
    ##
    ## read and format
    ##
    if (verbose) {
        message("Read and format using R")
    }
    data <- data.table::fread(table_file, sep = "\t", data.table = FALSE)
    vcf <- data.table::fread(cmd = paste0("gunzip -c ", vcf_output_file, ".gz | grep -v '^##' | cut -f1-9"), data.table =FALSE)
    ## argh no sep2 yet
    ## turn into cube
    K <- length(strsplit(data[1, 1], ",")[[1]])
    cube <- array(NA, c(nrow(data), ncol(data), K))
    for(i_sample in 1:ncol(data)) {
        cube[, i_sample, ] <- t(sapply(strsplit(data[, i_sample], ","), as.numeric))
    }
    dimnames(cube)[[1]] <- paste(
        vcf[, "#CHROM"],
        vcf[, "POS"],
        vcf[, "REF"],
        vcf[, "ALT"],
        sep = ":"
    )
    dimnames(cube)[[2]] <- dimnames(data)[[2]]
    dimnames(cube)[[3]] <- 1:K
    return(cube)
}

## cube <- extract_hd_to_cube(
##     vcf_file = "STITCH.Chr1.1.100000.vcf.gz",
##     ref = "../ref/tair10_um.fa",
##     bcftools = "bcftools",
##     gatk_jar = "/data/smew1/rdavies/stitch_richard_paper/bin/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar",
##     pos = NULL,
##     ram = "-Xmx4g"
## )

## cube <- extract_hd_to_cube(
##     vcf_file = "STITCH.Chr1.1.100000.vcf.gz",
##     ref = "../ref/tair10_um.fa",
##     bcftools = "bcftools",
##     gatk_jar = "/data/smew1/rdavies/stitch_richard_paper/bin/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar",
##     pos = cbind(
##         c("Chr1", "Chr1"),
##         c(11, 83)
##     ),
##     ram = "-Xmx4g",
##     samples = c("MAGIC.114", "MAGIC.14")
## )

