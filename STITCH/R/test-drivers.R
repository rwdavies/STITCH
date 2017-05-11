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

make_simple_bam <- function(
    file_stem,
    sam
) {
    cat(sam, file = paste0(file_stem, ".sam"))
    system2(
        "samtools",
        args = c(
            "view", "-bS", paste0(file_stem, ".sam"),
            ">", paste0(file_stem, ".bam")
        ),
        stderr = FALSE
    )
    file.remove(paste0(file_stem, ".sam"))
    system2("samtools", c("index", paste0(file_stem, ".bam")))
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
            "view", "-T", paste0(file_stem, ".fa"), "-C",
            "-o", paste0(file_stem, ".cram"),
            paste0(file_stem, ".sam")
        ),
        stderr = FALSE
    )
    file.remove(paste0(file_stem, ".sam"))
    system2("samtools", c("index", paste0(file_stem, ".cram")))
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
    phasemaster = NULL
) {
    set.seed(seed)
    if (is.na(K)) {
        phase <- array(
            round(sample(vals, size = n_snps * n_samples * 2, replace = TRUE)),
            c(n_snps, n_samples, 2)
        )
    } else {
        if (is.null(phasemaster))
            phasemaster <- array(round(sample(vals, size = n_snps * K, replace = TRUE)), c(n_snps, 2))
        phase <- array(0, c(n_snps, n_samples, 2))
        for(i_sample in 1:n_samples)
            for(i_hap in 1:2)
                phase[, i_sample, i_hap] <- phasemaster[, sample(K, 1)]
    }
    out <- sapply(1:n_samples, function(i_samp) {
        paste(phase[, i_samp, 1], phase[, i_samp, 2], sep = split_character)
    })
    if (is.null(header)) {
        colnames(out) <- paste0("samp", 1:n_samples)
    } else {
        colnames(out) <- header
    }
    write.table(out, file = phasefile, sep = "\t", row.names = FALSE, col.names = include_header, quote = FALSE)
    if (length(colnames(out)) == 1) {
        dimnames(phase)[[2]] <- list(colnames(out))
    } else {
        dimnames(phase)[[2]] <- colnames(out)
    }
    return(phase)
}



make_acceptance_test_data_package <- function(
    n_samples = 10,
    n_snps = 3,
    seed = 1,
    chr = 10,
    K = 2,
    n_reads = 4,
    phasemaster = NULL,
    reads_span_n_snps = NULL,
    n_cores = 1,
    tmpdir = tempdir(),
    use_crams = FALSE
) {

    if (is.null(reads_span_n_snps))
        reads_span_n_snps <- n_snps

    outputdir <- tempfile(pattern = "dir", tmpdir = tmpdir)
    if (use_crams) {
        rawdir <- file.path(outputdir, "crams")
    } else {
        rawdir <- file.path(outputdir, "bams")
    }
    dir.create(outputdir)
    dir.create(rawdir)

    posfile <- file.path(outputdir, "posfile.txt")
    pos <- make_posfile(posfile, n_snps = n_snps, chr = chr)

    phasefile <- file.path(outputdir, "phasefile.txt")
    phase <- make_phasefile(
        phasefile, n_snps = n_snps, n_samples = n_samples, K = K, seed = seed, phasemaster = phasemaster
    )

    genfile <- file.path(outputdir, "genfile.txt")
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
    read_names <- paste0("r00", 1:n_reads)
    bq <- paste0(rep(":", n), collapse = "")

    sample_names <- sapply(1:n_samples, function(i_sample)
        return(paste0("samp", i_sample)))

    sample_files <- lapply(1:n_samples, function(i_sample) {
        reads <- mclapply(
            1:n_reads,
            mc.cores = n_cores,
            function(i_read) {
            w <- sample(1:(n_snps - reads_span_n_snps + 1), 1) + 0:(reads_span_n_snps - 1)
            h <- phase[w, i_sample, (i_read %% 2) + 1]
            seq <- r[w]
            seq[h == 1] <- a[w][h == 1]
            seq <- paste0(seq, collapse = "")
            return(
                c(read_names[i_read], "0", chr, pos[w[1], 2], "60",
                  cigar, "*", "0", "0",
                  seq, bq)
            )
        })
        reads <- reads[order(as.integer(sapply(reads, function(x) x[[4]])))]
        sample_name <- sample_names[i_sample]
        if (use_crams) {
            out <- make_simple_cram(
                file_stem = file.path(rawdir, sample_name),
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
                file_stem = file.path(rawdir, sample_name),
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
        cramlist <- file.path(outputdir, "cramlist.txt")
        bamlist <- NULL
        samplelist <- cramlist
        ref <- sample_files[[1]]$ref        
        sample_files <- sapply(sample_files, function(x) x$cram_file)
        bam_files <- NULL
    } else {
        ref <- NULL
        cramlist <- NULL
        bamlist <- file.path(outputdir, "bamlist.txt")
        samplelist <- bamlist
        sample_files <- sapply(sample_files, function(x) x)
        bam_files <- sample_files
    }
    write.table(
        matrix(sample_files, ncol = 1),
        file = samplelist, row.names = FALSE, col.names = FALSE,
        quote = FALSE
    )

    return(list(
        bamlist = bamlist,
        cramlist = cramlist,
        ref = ref,
        posfile = posfile,
        genfile = genfile,
        bamlist = bamlist,
        chr = chr,
        outputdir = outputdir,
        phase = phase,
        L = as.integer(pos[, 2]),
        pos = pos,
        sample_names = sample_names,
        bam_files = bam_files,
        T = as.integer(nrow(pos)),
        phasefile = phasefile
    ))

}


check_vcf_against_phase <- function(
    vcf,
    phase,
    tol = 0.2
) {
    if (length(unique(vcf[, 9])) > 1)
        stop("not written for this")
    gt_names <- strsplit(vcf[1, 9], ":")[[1]]
    for(i_sample in 10:ncol(vcf)) {
        vcf_col <- vcf[, i_sample]
        vcf_col_split <- t(sapply(strsplit(vcf_col, ":"), I))
        colnames(vcf_col_split) <- gt_names
        ## check dosage
        dosage <- as.numeric(vcf_col_split[, "DS"])
        truth_dosage <- phase[, i_sample - 9, 1] + phase[, i_sample - 9, 2]
        expect_equal(sum(dosage < 0), 0)
        expect_equal(sum(dosage > 2), 0)
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

make_reference_package <- function(
    regionName = "test",
    n_snps = 10,
    n_samples_per_pop = 4,
    reference_populations = c("CEU", "GBR", "CHB"),
    L = NA,
    chr = 1,
    reference_sample_header = NA,
    reference_genders = c("male", "female")
) {

    n_total_samples <- length(reference_populations) * n_samples_per_pop
    posfile <- tempfile()
    pos <- make_posfile(
        posfile = posfile,
        n_snps = n_snps,
        seed = 1,
        L = L
    )

    reference_haplotype_file <- tempfile(fileext = ".gz")
    reference_legend_file <- tempfile(fileext = ".gz")
    reference_sample_file <- tempfile()

    ##ID POP GROUP SEX
    ##HG00096 GBR EUR male
    reference_samples <- array(0, c(length(reference_populations) * n_samples_per_pop, 4))
    colnames(reference_samples) <- c("ID", "POP", "GROUP", "SEX")
    reference_samples[, "POP"] <- rep(reference_populations, each = n_samples_per_pop)
    reference_samples[, "GROUP"] <- "NOT_USED"
    reference_samples[, "SEX"] <- rep(reference_genders, length.out = nrow(reference_samples))
    reference_samples[, "ID"] <- paste0("NOT_USED", 1:nrow(reference_samples))

    ## override the header for testing
    if (is.na(reference_sample_header[1]) == FALSE)
        colnames(reference_samples) <- reference_sample_header
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

    ## 0 0 1 1 - 2 entries per sample, row = SNP, one per sample
    reference_haplotypes <- array(
        sample(c(0, 1), n_snps * 2 * n_total_samples, replace = TRUE),
        c(n_snps, 2 * n_total_samples)
    )
    if (chr == "X") {
        ## 2nd hap for each male goes to -
        w <- 2 * which(reference_samples[, "SEX"] == "male")
        reference_haplotypes[, w] <- "-"
    }
    simple_write(reference_haplotypes, reference_haplotype_file, gzip = TRUE, col.names = FALSE)

    return(
        list(
            reference_haplotype_file = reference_haplotype_file,
            reference_sample_file = reference_sample_file,
            reference_legend_file = reference_legend_file,
            reference_populations = reference_populations,
            pos = pos,
            reference_haplotypes = reference_haplotypes,
            reference_legend = reference_legend,
            reference_samples = reference_samples
        )
    )

}
