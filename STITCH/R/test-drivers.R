#' @export
make_unique_tempdir <- function() {
    ## make every folder have a space!
    x <- tempfile(pattern = "folder", fileext = "wer wer2")
    dir.create(x)
    return(x)
}



skip_test_if_TRUE <- function(run_acceptance_tests) {
    if (run_acceptance_tests == FALSE)
        skip("skipping")
}


make_simple_sam_text <- function(
    entries = NULL,
    chr = 1,
    sample_name = "sample_name",
    chr_length = 45,
    include_rg_tag = TRUE,
    include_rg_tag_with_no_sm = FALSE
) {
    out <- paste0(
        "@HD\tVN:1.5\tSO:coordinate\n",
        "@SQ\tSN:", chr, "\tLN:", chr_length, "\n"
    )
    if (include_rg_tag)
        out <- paste0(out, "@RG\tID:7369_8x15\tSM:", sample_name, "\n")
    if (include_rg_tag_with_no_sm)
        out <- paste0(out, "@RG\tID:7369_8x15\n")
    if (length(entries) > 0)
        for(entry in entries)
            out <- paste0(out, paste0(entry, collapse = "\t"), "\n")
    return(out)
}

#' @export
make_simple_bam <- function(
    file_stem,
    sam
) {
    cat(sam, file = paste0(file_stem, ".sam"))
    out <- system2(
        "samtools",
        args = c(
            "view", "-bS", shQuote(paste0(file_stem, ".sam")),
            ">", shQuote(paste0(file_stem, ".bam"))
        ),
        stderr = FALSE
    )
    file.remove(paste0(file_stem, ".sam"))
    system2("samtools", c("index", shQuote(paste0(file_stem, ".bam"))))
    return(paste0(file_stem, ".bam"))
}

make_ref_from_sam <- function(sam, pos = NULL) {
    h <- strsplit(sam, "\n")[[1]]
    b <- strsplit(h[grep("@SQ", h)], "\t")[[1]]
    chr_name <- strsplit(b[grep("SN", b)], "SN:")[[1]][2]
    chr_length <- as.integer(strsplit(b[grep("LN", b)], "LN:")[[1]][2])
    ref <- paste0(rep("A", chr_length), collapse = "")
    if (is.null(pos) == FALSE) {
        for(i_row in 1:nrow(pos)) {
            x <- as.integer(pos[i_row, 2])
            substr(ref, x, x) <- as.character(pos[i_row, "REF"])
        }
    }
    ref_fa <- paste0(">", chr_name, "\n", ref, "\n")
    return(ref_fa)
}

## for now - just make simple ref with As
## later, do more intelligent thing of making sure ref is in agreement
make_simple_cram <- function(
    file_stem,
    sam,
    pos = NULL
) {
    cat(make_ref_from_sam(sam, pos), file = paste0(file_stem, ".fa"))
    cat(sam, file = paste0(file_stem, ".sam"))
    system2(
        "samtools",
        args = c(
            "view", "-T", shQuote(paste0(file_stem, ".fa")), "-C",
            "-o", shQuote(paste0(file_stem, ".cram")),
            shQuote(paste0(file_stem, ".sam"))
        ),
        stderr = FALSE
    )
    file.remove(paste0(file_stem, ".sam"))
    system2("samtools", c("index", shQuote(paste0(file_stem, ".cram"))))
    return(
        list(
            cram_file = paste0(file_stem, ".cram"),
            ref = paste0(file_stem, ".fa")
        )
    )
}

write_names_to_disk <- function(
    bam_names,
    bamlist
) {
    write.table(
        matrix(bam_names, ncol = 1),
        file = bamlist,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}


make_posfile <- function(
    posfile,
    n_snps = 10,
    refs = NA,
    alts = NA,
    L = NA,
    chr = NA,
    seed = 1
) {
    set.seed(seed)
    if (is.na(chr[1])) chr <- rep(1, n_snps)
    if (is.na(L[1])) L <- 1:n_snps
    if (is.na(refs[1])) refs <- rep("A", n_snps)
    if (is.na(alts[1])) alts <- rep("G", n_snps)
    pos <- data.frame(chr, L, refs, alts)
    colnames(pos) <- NULL
    write.table(pos, file = posfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    colnames(pos) <- c("CHR", "POS", "REF", "ALT")
    return(pos)
}

make_genfile <- function(
    genfile,
    n_snps = 10,
    n_samples = 3,
    include_header = TRUE,
    vals = c(0, 1, 2),
    seed = 1,
    header = NULL
) {
    set.seed(seed)
    gen <- array(
        round(sample(vals, size = n_snps * n_samples, replace = TRUE)),
        c(n_snps, n_samples)
    )
    if (is.null(header)) {
        colnames(gen) <- paste0("samp", 1:n_samples)
    } else {
        colnames(gen) <- header
    }
    write.table(gen, file = genfile, sep = "\t", row.names = FALSE, col.names = include_header, quote = FALSE)
    return(gen)
}


make_phasefile <- function(
    phasefile,
    n_snps = 10,
    n_samples = 3,
    include_header = TRUE,
    split_character = "|",
    vals = c(0, 1),
    seed = 1,
    header = NULL,
    K = NA,
    phasemaster = NULL,
    samples_are_inbred = FALSE,
    write_row_as_NA = NULL ## boolean vector when defined
) {
    set.seed(seed)
    if (is.na(K)) {
        phase <- array(
            round(sample(vals, size = n_snps * n_samples * 2, replace = TRUE)),
            c(n_snps, n_samples, 2)
        )
    } else {
        if (is.null(phasemaster))
            phasemaster <- array(round(sample(vals, size = n_snps * K, replace = TRUE)), c(n_snps, K))
        phase <- array(0, c(n_snps, n_samples, 2))
        for(i_sample in 1:n_samples) {
            if (samples_are_inbred) {
                haps_to_sample <- rep(sample(K, 1), 2)
            } else {
                haps_to_sample <- c(sample(K, 1), sample(K, 1))
            }
            for(i_hap in 1:2) {
                phase[, i_sample, i_hap] <- phasemaster[, haps_to_sample[i_hap]]
            }
        }
    }
    out <- sapply(1:n_samples, function(i_samp) {
        paste(phase[, i_samp, 1], phase[, i_samp, 2], sep = split_character)
    })
    if (is.null(header)) {
        colnames(out) <- paste0("samp", 1:n_samples)
    } else {
        colnames(out) <- header
    }
    if (!is.null(write_row_as_NA)) {
        if (nrow(out) != length(write_row_as_NA)) {
            stop("Supplied write_row_as_NA must be the same dimension as phase")
        }
        out[write_row_as_NA, ] <- "NA|NA"
    }
    write.table(out, file = phasefile, sep = "\t", row.names = FALSE, col.names = include_header, quote = FALSE)
    if (length(colnames(out)) == 1) {
        dimnames(phase)[[2]] <- list(colnames(out))
    } else {
        dimnames(phase)[[2]] <- colnames(out)
    }
    return(phase)
}


#' @export
make_acceptance_test_data_package <- function(
    n_samples = 10,
    n_snps = 3,
    seed = 1,
    chr = 10,
    K = 2,
    n_reads = 4,
    L = NA,
    phasemaster = NULL,
    reads_span_n_snps = NULL,
    n_cores = 1,
    tmpdir = tempdir(),
    use_tmpdir_directly = FALSE,
    use_crams = FALSE,
    sample_names = NULL,
    samples_are_inbred = FALSE,
    phred_bq = 25,
    regionName = NA,
    write_row_as_NA = NULL,
    use_bx_tag = FALSE,
    refs = NA,
    alts = NA
) {

    if (length(n_reads) == 1) {
        n_reads <- rep(round(n_reads), n_samples)
    }

    if (is.null(reads_span_n_snps)) {
        reads_span_n_snps <- n_snps
    }

    if (is.na(regionName)) {
        if (is.na(L[1])) {
            regionName <- chr
        } else {
            regionName <- paste0(chr, ".", min(L), ".", max(L))
        }
    }

    if (!use_tmpdir_directly) {
        ## default - build new, specific directory for this
        outputdir <- tempfile(pattern = "dir", tmpdir = tmpdir)
    } else {
        outputdir <- tmpdir
    }
    if (use_crams) {
        rawdir <- file.path(outputdir, "crams")
    } else {
        rawdir <- file.path(outputdir, "bams")
    }
    if (!dir.exists(outputdir)) {
        dir.create(outputdir)
    }

    if (!dir.exists(rawdir)) {
        dir.create(rawdir)
    }

    if (is.na(L[1])) {
        L_is_simple <- TRUE
    } else {
        L_is_simple <- FALSE
    }
    posfile <- file.path(outputdir, paste0("pos.", regionName, ".txt"))
    pos <- make_posfile(
        posfile,
        L = L,
        n_snps = n_snps,
        chr = chr,
        refs = refs,
        alts = alts
    )
    L <- pos[, 2]

    phasefile <- file.path(outputdir, paste0("phase.", regionName, ".txt"))
    
    phase <- make_phasefile(
        phasefile,
        n_snps = n_snps,
        n_samples = n_samples,
        K = K,
        seed = seed,
        phasemaster = phasemaster,
        samples_are_inbred = samples_are_inbred,
        write_row_as_NA = write_row_as_NA
    )

    genfile <- file.path(outputdir, paste0("gen.", regionName, ".file.txt"))
    write.table(
        phase[, , 1] + phase[, , 2] ,
        file = genfile,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    r <- as.character(pos[, "REF"])
    a <- as.character(pos[, "ALT"])
    n <- reads_span_n_snps
    cigar <- paste0(n, "M")
    phred_bq_char <-  rawToChar(as.raw(c(phred_bq + 33)))
    bq <- paste0(rep(phred_bq_char, n), collapse = "")

    if (is.null(sample_names)) {
        sample_names <- sapply(1:n_samples, function(i_sample)
            return(paste0("samp", i_sample)))
    } else if (length(sample_names) != n_samples) {
        stop(paste0("You idiot, length(sample_names) = ", length(sample_names), " is not the same as n_samples which is ", n_samples))
    }

    sample_files <- lapply(1:n_samples, function(i_sample) {
        to_sample <- 1:(n_snps - reads_span_n_snps + 1)
        out <- mclapply(
            1:n_reads[i_sample],
            mc.cores = n_cores,
            simulate_a_read,
            n_snps = n_snps,
            reads_span_n_snps = reads_span_n_snps,
            i_sample = i_sample,
            phase = phase,
            seq = seq,
            r = r,
            a = a,
            phred_bq_char = phred_bq_char,
            L_is_simple = L_is_simple,
            bq = bq,
            cigar = cigar,
            chr = chr,
            pos = pos,
            L = L
        )
        reads <- lapply(out, function(x) x[["sam_bit"]])
        o <- order(as.integer(sapply(reads, function(x) x[[4]])))
        reads <- reads[o]
        if (use_bx_tag) {
            ## combine a few of them
            hap <- sapply(out, function(x) x[["ihap"]])
            hap <- hap[o]
            for(ihap in 1:2) {
                ## hmm, choose 2 groups of >= 4, to have 2 qnames, and 2 bx tagsb
                g <- sample(which(hap == ihap), 8)
                for(j in 1:2) {
                    off <- 4 * (j - 1)
                    bxtag <- paste0("BX:A0", j, "B0", ihap)
                    reads[[g[off + 1]]][1] <- paste0("read", j, "_", ihap)
                    reads[[g[off + 1]]][12] <- bxtag
                    reads[[g[off + 1 + 1]]][1] <- paste0("read", j, "_", ihap)
                    reads[[g[off + 1 + 1]]][12] <- bxtag
                    ##
                    reads[[g[off + 1 + 2]]][1] <- paste0("read", j, "B_", ihap)
                    reads[[g[off + 1 + 2]]][12] <- bxtag
                    reads[[g[off + 1 + 1 + 2]]][1] <- paste0("read", j, "B_", ihap)
                    reads[[g[off + 1 + 1 + 2]]][12] <- bxtag
                }
            }
        }
        sample_name <- sample_names[i_sample]
        key <- round(10e4 * runif(1))
        file_stem <- file.path(rawdir, paste0(sample_name, ".", key))
        if (use_crams) {
            out <- make_simple_cram(
                file_stem = file_stem,
                sam = make_simple_sam_text(
                    reads,
                    chr,
                    sample_name = sample_name,
                    chr_length = max(pos[, 2]) + 1
                ),
                pos = pos
            )
        } else {
            out <- make_simple_bam(
                file_stem = file_stem,
                sam = make_simple_sam_text(
                    reads,
                    chr,
                    sample_name = sample_name,
                    chr_length = max(pos[, 2]) + 1
                )
            )
        }
        return(out)
    })

    if (use_crams) {
        cramlist <- file.path(outputdir, paste0("cramlist.", regionName, ".txt"))
        bamlist <- NULL
        samplelist <- cramlist
        ref <- sample_files[[1]]$ref
        sample_files <- sapply(sample_files, function(x) x$cram_file)
        bam_files <- NULL
    } else {
        ref <- NULL
        cramlist <- NULL
        bamlist <- file.path(outputdir, paste0("bamlist.", regionName, ".txt"))
        samplelist <- bamlist
        sample_files <- sapply(sample_files, function(x) x)
        bam_files <- sample_files
    }
    write.table(
        matrix(sample_files, ncol = 1),
        file = samplelist, row.names = FALSE, col.names = FALSE,
        quote = FALSE
    )

    return(
        list(
            bamlist = bamlist,
            cramlist = cramlist,
            ref = ref,
            posfile = posfile,
            genfile = genfile,
            chr = chr,
            outputdir = outputdir,
            phase = phase,
            L = as.integer(pos[, 2]),
            pos = pos,
            sample_names = sample_names,
            bam_files = bam_files,
            nSNPs = as.integer(nrow(pos)),
            phasefile = phasefile
        )
    )

}


check_output_against_phase <- function(
    file,
    data_package,
    output_format,
    which_snps = NULL,
    tol = 0.2,
    who = NULL,
    min_info = 0.98
) {
    if (is.null(which_snps)) {
        which_snps <- 1:length(data_package$L)
    }
    if (is.null(who)) {
        who <- 1:dim(data_package$phase)[2]
    }
    if (substr(file, nchar(file) - 4, nchar(file)) == ".bgen") {
        out <- rrbgen::rrbgen_load(bgen_file = file)
        ## check what was written to .bgen
        expect_equal(as.character(out$var_info[, "chr"]), as.character(data_package$pos[which_snps, "CHR"]))
        expect_equal(as.character(out$var_info[, "position"]), as.character(data_package$pos[which_snps, "POS"]))
        expect_equal(as.character(out$var_info[, "ref"]), as.character(data_package$pos[which_snps, "REF"]))
        expect_equal(as.character(out$var_info[, "alt"]), as.character(data_package$pos[which_snps, "ALT"]))
        ##
        check_bgen_gp_against_phase(
            gp = out$gp,
            phase = data_package$phase,
            which_snps = which_snps,
            tol = tol,
            who = who
        )
        var_info <- read.table(paste0(file, ".per_snp_stats.txt.gz"), header = TRUE)
        expect_equal(nrow(var_info), nrow(out$gp))
        if (sum((var_info[, "INFO_SCORE"] < min_info)) > 0) {
            print(paste0("min_info = ", min_info))
            print(var_info[, "INFO_SCORE"])
        }
        expect_equal((var_info[, "INFO_SCORE"] >= min_info), rep(TRUE, nrow(var_info)))
    } else {
        ##
        vcf <- read.table(
            file,
            header = FALSE,
            stringsAsFactors = FALSE
        )
        ## check imputation worked here
        check_vcf_against_phase(
            vcf,
            data_package$phase,
            which_snps,
            tol = tol,
            who = who
        )
        check_vcf_info_scores(
            vcf = vcf,
            min_info = min_info
        )
    }
    return(NULL)
}

#' @export
check_vcf_info_scores <- function(vcf, min_info = 0.98) {
    info_scores <- sapply(1:nrow(vcf), function(i_snp) {
        x <- strsplit(vcf[i_snp, 8], ";")[[1]]
        y <- grep("INFO_SCORE", x)
        expect_true(length(y) == 1)
        info_score <- as.numeric(strsplit(x[y], "INFO_SCORE=")[[1]][2])
        return(info_score)
    })
    expect_equal(info_scores >= min_info, rep(TRUE, length(info_scores)))
}

check_vcf_against_phase <- function(
    vcf,
    phase,
    which_snps,
    who = NULL,
    tol = 0.2
) {
    if (length(unique(vcf[, 9])) > 1) {
        stop("not written for this")
    }
    if (is.null(who)) {
        who <- 1:(ncol(vcf) - 9)
    }
    gt_names <- strsplit(vcf[1, 9], ":")[[1]]
    for(i_sample in (9 + who)) {
        vcf_col <- vcf[, i_sample]
        vcf_col_split <- t(sapply(strsplit(vcf_col, ":"), I))
        colnames(vcf_col_split) <- gt_names
        ## check dosage
        dosage <- as.numeric(vcf_col_split[, "DS"])
        truth_dosage <- phase[which_snps, i_sample - 9, 1] + phase[which_snps, i_sample - 9, 2]
        expect_equal(sum(dosage < 0), 0)
        expect_equal(sum(dosage > 2), 0)
        ## print to screen first
        if (sum(abs(dosage - truth_dosage) > tol) > 0) {
            print(paste0("i_sample = ", i_sample))
            print(cbind("dosage" = dosage, "truth_dosage" = truth_dosage))
        }
        expect_equal(max(abs(dosage - truth_dosage)) <= tol, TRUE)
        ## check genotype probability
        genotype_posteriors <- t(sapply(strsplit(vcf_col_split[, "GP"], ","), as.numeric))
        r <- rowSums(genotype_posteriors)
        ## check their sum, up to a tolerance
        expect_equal(sum(abs(r - 1) > 0.00101), 0)
        expect_equal(sum(genotype_posteriors < 0), 0)
        expect_equal(sum(genotype_posteriors > 1), 0)
        d2 <- genotype_posteriors[, 2] + 2 * genotype_posteriors[, 3]
        expect_equal(sum(abs(d2 - dosage) > 0.00101), 0)
        ## check ancestral haplotype dosages (if applicable)
        if ("HD" %in% gt_names) {
            q_t <- sapply(strsplit(vcf_col_split[, "HD"], ","), as.numeric)
            y <- sum(abs(colSums(q_t) - 2) > 0.01)
            if (y > 0) {
                print(q_t)
            }
            expect_equal(y, 0)
        }
    }
}

check_bgen_gp_against_phase <- function(
    gp,
    phase,
    which_snps,
    who = NULL,
    tol = 0.2
) {
    if (is.null(who)) {
        who <- 1:dim(gp)[2]
    }
    for(i_sample in who) {
        ## check genotype probability
        genotype_posteriors <- array(NA, c(dim(gp)[1], 3))
        genotype_posteriors[, ] <- gp[, i_sample, , drop = FALSE]
        r <- rowSums(genotype_posteriors)
        ## check their sum, up to a tolerance
        expect_equal(sum(abs(r - 1) > 0.00101), 0)
        expect_equal(sum(genotype_posteriors + 1e-8 < 0), 0)
        expect_equal(sum(genotype_posteriors - 1e-8 > 1), 0)
        gp_dosage <- genotype_posteriors[, 2] + 2 * genotype_posteriors[, 3]
        truth_dosage <- phase[which_snps, i_sample, 1] + phase[which_snps, i_sample, 2]
        ## pre-check
        if (sum(abs(gp_dosage - truth_dosage) > tol) > 0) {
            print(paste0("i_sample = ", i_sample))
            print(cbind("gp_dosage" = gp_dosage, "truth_dosage" = truth_dosage))
        }
        expect_equal(max(abs(gp_dosage - truth_dosage)) <= tol, TRUE)
    }
}



simple_write <- function(matrix, file, gzip = FALSE, col.names = TRUE) {
    if (gzip)
        file <- gzfile(file, "w")
    write.table(
        matrix,
        file = file,
        col.names = col.names,
        row.names = FALSE,
        sep = " ",
        quote = FALSE
    )
    if (gzip)
        close(file)
}

#' @export
make_reference_package <- function(
    n_snps = 10,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    L = NA,
    chr = 1,
    reference_sample_header = NA,
    reference_genders = c("male", "female"),
    phasemaster = NULL,
    tmpdir = tempdir(),
    use_tmpdir_directly = FALSE,
    regionName = NA,
    expRate = 0.5
) {

    if (is.na(regionName)) {
        if (is.na(L[1])) {
            regionName <- chr
        } else {
            regionName <- paste0(chr, ".", min(L), ".", max(L))
        }
    }
    if (!use_tmpdir_directly) {
        ## default - build new, specific directory for this
        outputdir <- tempfile(pattern = "dir", tmpdir = tmpdir)
    } else {
        outputdir <- tmpdir
    }
    if (!dir.exists(outputdir)) {
        dir.create(outputdir)
    }

    ## add a space in it
    reference_vcf_file <- file.path(outputdir, paste0("ref hap.", regionName, ".vcf.gz"))
    reference_haplotype_file <- file.path(outputdir, paste0("ref hap.", regionName, ".txt.gz"))
    reference_legend_file <- file.path(outputdir, paste0("ref legend.", regionName, ".txt.gz"))
    reference_sample_file <- file.path(outputdir, paste0("ref sample.", regionName, ".txt"))
    reference_genetic_map_file <- file.path(outputdir, paste0("ref gen.", regionName, ".txt.gz"))
    posfile <- file.path(outputdir, paste0("ref.", regionName, ".pos.txt"))
    
    ##     
    n_total_samples <- length(reference_populations) * n_samples_per_pop
    pos <- make_posfile(
        posfile = posfile,
        n_snps = n_snps,
        seed = 1,
        L = L,
        chr = chr
    )

    ##ID POP GROUP SEX
    ##HG00096 GBR EUR male
    reference_samples <- array(0, c(length(reference_populations) * n_samples_per_pop, 4))
    colnames(reference_samples) <- c("ID", "POP", "GROUP", "SEX")
    reference_samples[, "POP"] <- rep(reference_populations, each = n_samples_per_pop)
    reference_samples[, "GROUP"] <- "NOT_USED"
    reference_samples[, "SEX"] <- rep(reference_genders, length.out = nrow(reference_samples))
    reference_samples[, "ID"] <- paste0("SAMPLE", 1:nrow(reference_samples))

    ## override the header for testing
    if (is.na(reference_sample_header[1]) == FALSE) {
        colnames(reference_samples) <- reference_sample_header
    }
    simple_write(reference_samples, reference_sample_file)

    ##id position a0 a1 TYPE AFR AMR EAS EUR SAS ALL
    ##20:60343:G:A 60343 G A Biallelic_SNP 0 0.00144092219020173 0 0 0 0.000199680511182109
    reference_legend <- data.frame(
        id = "NOT_USED",
        position = pos[, "POS"],
        a0 = pos[, "REF"],
        a1 = pos[, "ALT"]
    )
    simple_write(reference_legend, reference_legend_file, gzip = TRUE)

    ## either sample at random
    ## or sample from phasemaster without recomb
    reference_haplotypes <- array(NA, c(n_snps, 2 * n_total_samples))
    for (i_sample in 1:n_total_samples) {
        for (i_hap in 1:2) {
            c <- 2 * (i_sample - 1) + i_hap
            if (is.null(phasemaster )) {
                g <- sample(c(0L, 1L), n_snps, replace = TRUE)
            } else {
                g <- phasemaster[, sample(1:ncol(phasemaster), 1)]
            }
            reference_haplotypes[, c] <- g
        }
    }

    if (chr == "X") {
        ## 2nd hap for each male goes to -
        w <- 2 * which(reference_samples[, "SEX"] == "male")
        reference_haplotypes[, w] <- "-"
    }

    simple_write(reference_haplotypes, reference_haplotype_file, gzip = TRUE, col.names = FALSE)

    ## make genetic map as well
    ## do most SNPs, assume 1 cM / Mbp
    ## assume
    genetic_map <- make_genetic_map_file(L = L, n_snps = n_snps, expRate = expRate)
    simple_write(genetic_map, reference_genetic_map_file, gzip = TRUE, col.names = TRUE)

    ## not sure how important long term    
    if (is.na(reference_sample_header[1])) {
        colClasses <- get_reference_colClasses(
            reference_sample_file = reference_sample_file,
            reference_populations = reference_populations,
            chr = chr
        )
    } else {
        ## not NA, these have been set, and will fail the above
        colClasses <- NULL
    }
    
    ## do reference VCF here
    make_and_write_reference_vcf(pos, reference_haplotypes, reference_samples, reference_vcf_file) 
    
    return(
        list(
            reference_vcf_file = reference_vcf_file,
            reference_haplotype_file = reference_haplotype_file,
            reference_sample_file = reference_sample_file,
            reference_legend_file = reference_legend_file,
            reference_genetic_map_file = reference_genetic_map_file,
            reference_populations = reference_populations,
            pos = pos,
            reference_haplotypes = reference_haplotypes,
            reference_legend = reference_legend,
            reference_samples = reference_samples,
            colClasses = colClasses
        )
    )

}



make_and_write_reference_vcf <- function(pos, reference_haplotypes, reference_samples, reference_vcf_file) {
    
    reference_vcf_file_unzipped <- gsub(".gz", "", reference_vcf_file)
    sampleNames <- reference_samples[, 1]
    header <- paste0(
        '##fileformat=VCFv4.1\n',
        '##FILTER=<ID=PASS,Description="All filters passed">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n',
        paste0('##contig=<ID=', pos[1, 1], '>', sep = ""), "\n"
    )
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    cat(header, header2, "\n", sep="", file = reference_vcf_file_unzipped)

    ## make rest
    N <- length(sampleNames)
    vcf_matrix_to_out <- data.frame(matrix(
        data = NA,
        nrow = nrow(pos),
        ncol = N + 9
    ))
    
    vcf_matrix_to_out[, 1] <- pos[, 1]
    vcf_matrix_to_out[, 2] <- pos[, 2]
    vcf_matrix_to_out[, 3] <- "."
    vcf_matrix_to_out[, 4] <- pos[, 3]
    vcf_matrix_to_out[, 5] <- pos[, 4]
    vcf_matrix_to_out[, 6] <- "."
    vcf_matrix_to_out[, 7] <- "PASS"
    vcf_matrix_to_out[, 8] <- "."
    vcf_matrix_to_out[, 9] <- "GT"
    
    for(iSample in 1:length(sampleNames)) {
        i <- 2 * iSample - 1
        vcf_matrix_to_out[, 9 + iSample] <- paste0(reference_haplotypes[, i], "|", reference_haplotypes[, i + 1])
    }

    write.table(
        vcf_matrix_to_out,
        file = reference_vcf_file_unzipped,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        append = TRUE,
        quote = FALSE
    )

    system(paste0("bgzip -f ", shQuote(reference_vcf_file_unzipped)))
    system(paste0("tabix ", shQuote(reference_vcf_file)))

    NULL
    
}




##position COMBINED_rate(cM/Mb) Genetic_Map(cM)
##150118 1.13462264157027 0
##154675 1.12962782559127 0.00517047537763574
##154753 1.13654510133156 0.00525858634803186
##168567 1.58657526542862 0.0209588203778261

#' @export
make_genetic_map_file <- function(L, n_snps, expRate = 0.5) {
    if (is.na(L[1])) {
        L <- 1:n_snps
    }
    genetic_map <- array(NA, c(n_snps, 3))
    colnames(genetic_map) <- c("position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM.")
    genetic_map[, "position"] <- L
    genetic_map[, "COMBINED_rate.cM.Mb."] <- expRate
    genetic_map[n_snps, "COMBINED_rate.cM.Mb."] <- 0
    ## fill in - be simple!
    genetic_map <- fill_in_genetic_map_cm_column(genetic_map)
    return(genetic_map)
}



simulate_a_read <- function(
    i_read,
    n_snps,
    reads_span_n_snps,
    i_sample,
    phase,
    seq,
    r,
    a,
    phred_bq_char,
    L_is_simple,
    bq,
    cigar,
    chr,
    pos,
    L
) {
    ## choose which SNP to intersect, then choose a position
    if (L_is_simple) {
        w <- sample(1:(n_snps - reads_span_n_snps + 1), 1) + 0:(reads_span_n_snps - 1)
        ihap <- (i_read %% 2) + 1
        h <- phase[w, i_sample, ihap]
        seq <- r[w]
        seq[h == 1] <- a[w][h == 1]
        seq <- paste0(seq, collapse = "")
        seq <- r[w]
        seq[h == 1] <- a[w][h == 1]
        seq <- paste0(seq, collapse = "")
        local_bq <- bq
        local_cigar <- cigar
    } else {
        snp_to_intersect <- sample(n_snps, 1)
        which_in_intercept <- sample(reads_span_n_snps, 1) ## 1-based
        ##print("------------------------")
        ##print(snp_to_intersect)
        ##print(which_in_intercept)
        w <- snp_to_intersect + (0:(reads_span_n_snps - 1)) - (which_in_intercept - 1)
        ## now, if out of bounds, nudge back in
        if (sum(w <= 0) > 0) {
            w <- w + -min(w)+ 1
        }
        if (sum(w > n_snps) > 0) {
            w <- w - (max(w) - n_snps)
        }
        ## w is 1-based sampling of start to end
        ## w <- sample(to_sample, 1) + 0:(reads_span_n_snps - 1)
        ihap <- sample(2, 1)
        h <- phase[w, i_sample, ihap]
        ## make "A" otherwise
        seq <- rep("A", tail(L[w], 1) - head(L[w], 1) + 1)
        seq_w <- L[w] - min(L[w]) + 1
        seq[seq_w] <- r[w]
        seq[seq_w][h == 1] <- a[w][h == 1]
        n <- length(seq)
        seq <- paste0(seq, collapse = "")
        local_bq <- paste0(rep(phred_bq_char, n), collapse = "")
        local_cigar <- paste0(n, "M")
    }
    return(
        list(
            sam_bit = c(paste0("r00", i_read), "0", chr, pos[w[1], 2], "60",
              local_cigar, "*", "0", "0",
              seq, local_bq),
            ihap = ihap
        )
    )
}


#' @export
make_fb_test_package <- function(
    K = 4,
    nReads = 8,
    nSNPs = 10,
    S = 2,
    gridWindowSize = 3
) {
    set.seed(4916)
    L <- 1:nSNPs
    maxDifferenceBetweenReads <- 1e10
    maxEmissionMatrixDifference <- 1e50
    Jmax_local <- 1e3
    ##
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    snps_in_grid_1_based <- out$snps_in_grid_1_based
    ## make almost the same
    eHapsCurrent_tc <- array(NA, c(K, nSNPs, S))
    alphaMatCurrent_tc <- array(NA, c(K, nGrids - 1, S))
    eHapsCurrent_t <-
        0.1 * array(runif(K * nSNPs), c(K, nSNPs)) +
        0.9 * array(sample(c(0, 1), K * nSNPs, replace = TRUE), c(K, nSNPs))
    for(s in 1:S) {
        eHapsCurrent_tc[, , s] <- 0.90 * eHapsCurrent_t + 0.10 * runif(nSNPs * K)
        m <- array(runif(K * (nGrids - 1)), c(K, (nGrids - 1)))
        alphaMatCurrent_tc[, , s] <- t(t(m) / colSums(m))
    }
    sigmaCurrent_m <- array(0.9 + 0.1 * runif((nGrids - 1) * S), c(nGrids - 1, S))
    priorCurrent_m <- array(1 / K, c(K, S))
    ##
    transMatRate_tc_D <- get_transMatRate_m("diploid", sigmaCurrent_m)
    transMatRate_tc_H <- get_transMatRate_m("diploid-inbred", sigmaCurrent_m)
    ## sample reads
    cr <- sort(sample(1:nGrids, nReads, replace = TRUE)) - 1
    true_H <- sample(c(1, 2), nReads, replace = TRUE)
    ## choose them from haps 1 and 2
    sampleReads <- lapply(
        1:length(cr),
        function(ii) {
        i <- cr[ii] ## 0-based
        w <- which(grid == i)
        g <- w[sample(1:length(w), size = 1)] - 1
        x <- c(g - 2, g - 1, g, g + 1, g + 2)
        x <- x[x %in% 0:(nSNPs - 1)] ## 0-based
        i_hap <- true_H[ii]
        bq <- matrix(round(30 * (eHapsCurrent_tc[i_hap, x + 1, 1] - 0.5)), ncol = 1)
        list(
            length(x) - 1,
            i,
            bq,
            matrix(x,ncol=1,nrow=length(x))
        )
    })
    N <- 1
    ## almost certainly want to pre-declare
    list_of_eMatRead_t <- lapply(0:(S - 1), function(s) {
        eMatRead_t <- array(1, c(K, nReads))
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t,
            sampleReads = sampleReads,
            eHapsCurrent_tc = eHapsCurrent_tc,
            s = s,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = Jmax_local,
            eMatHapOri_t = array(0, c(1, 1)), ## ugh
            pRgivenH1 = array(0),
            pRgivenH2 = array(0),
            run_pseudo_haploid = FALSE,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text"
        )
        return(eMatRead_t)
    })
    ##
    list_of_eMatGrid_t <- lapply(0:(S - 1), function(s) {
        eMatGrid_t <- array(1, c(K, nGrids))
        rcpp_make_eMatGrid_t(
            eMatGrid_t = eMatGrid_t,
            eMatRead_t = list_of_eMatRead_t[[s + 1]],
            H = 1,
            sampleReads = sampleReads,
            hap = 1,
            nGrids = nGrids,
            run_fb_grid_offset = 0,
            use_all_reads = TRUE,
            bound = TRUE,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text"
        )
        return(eMatGrid_t)
    })
    ##
    return(
        list(
            eHapsCurrent_tc = eHapsCurrent_tc,
            alphaMatCurrent_tc = alphaMatCurrent_tc,
            sigmaCurrent_m = sigmaCurrent_m,
            priorCurrent_m = priorCurrent_m,
            transMatRate_tc_D = transMatRate_tc_D,
            transMatRate_tc_H = transMatRate_tc_H,
            S = S,
            K = K,
            gridWindowSize = gridWindowSize,
            nReads = nReads,
            nSNPs = nSNPs,
            L = L,
            N = N,
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based,
            sampleReads = sampleReads,
            true_H = true_H,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax_local = Jmax_local,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            list_of_eMatRead_t = list_of_eMatRead_t,
            list_of_eMatGrid_t = list_of_eMatGrid_t
        )
    )
}


random_R_version_checker <- function() {
    minor1 <- as.numeric(strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1])
    minor2 <- as.numeric(strsplit(R.version$minor, ".", fixed = TRUE)[[1]][2])    
    if (as.numeric(R.version$major) == 3) {
        if (minor1 <= 5) {
            return(1)
        } else {
            return(2)
        }
    }
}


