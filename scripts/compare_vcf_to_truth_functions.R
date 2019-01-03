## get all MegaMuga samples for a chromosome
get_megamuga_calls_for_chr <- function(chr) {
    if (verbose)
        message("Get megamuga calls for chr")

    samples <- read.table(file.path(megamugadir, "Sample_Map.txt"), sep="\t", header = TRUE)
    oxford <- is.na(as.numeric(as.character(samples[,"Name"])))==TRUE
    oxford[32]=TRUE # it takes the e as a number!
    ## re-do the oxford ones
    samples <- cbind(samples,"")
    colnames(samples)[ncol(samples)]="betterName"
    samples[,"betterName"]=as.character(samples[,"Name"])
    for(w in (1:length(oxford))[oxford]) {
        x=as.character(samples[w,"betterName"])
        y=unlist(strsplit(x,""))
        z=c(y[1:(length(y))-1],".0",y[length(y)])
        e=paste(z,collapse="")
        samples[w,"betterName"]=e
    }
    N <- nrow(samples)
    snps <- read.table(file.path(megamugadir, "SNP_Map.txt"), sep = "\t", header = TRUE)
    nSNP <- nrow(snps)

    ##
    ## load in all samples
    ##
    in_file <- file.path(megamugadir, "Univ_of_Chicago_St_Pierre_MEGMUGV01_20150505_FinalReport.txt")
    data <- fread(
        in_file,
        data.table = FALSE,
        showProgress = FALSE,
        skip = 9
    )

    ##
    ## turn into a matrix
    ##
    head <- read.table(
        file.path(megamugadir, "Univ_of_Chicago_St_Pierre_MEGMUGV01_20150505_FinalReport.txt"),
        skip = 9, nrow = 1, sep = "\t"
    )
    pos <- array(0, c(nSNP, 4))
    colnames(pos) <- c("CHR","POS","REF","ALT")
    options(scipen=999)
    pos[,1:2] <- cbind(as.character(snps[,"Chromosome"]),as.numeric(as.character(snps[,"Position"])))
    pos[,3] <- unlist(lapply(snps[,"SNP"],function(x) substr(x,2,2)))
    pos[,4] <- unlist(lapply(snps[,"SNP"],function(x) substr(x,4,4)))
    g <- array(0,c(nSNP,N))
    gp <- array(0,c(nSNP,N,3)) # the cluster probabilities
    colnames(g) <- samples[,"Name"]
    rownames(g) <- snps[,"Name"]

    ## try filling it in
    for(i in 1:N) {
        samp <- colnames(g)[i]
        x <- data[data[, "Sample ID"] == samp, ]
        y <- match(x[, 1], rownames(g))
        g[y, i] <- as.integer(x[,7]=="B") + as.integer(x[,8]=="B")
        g[y, i][x[,9]<0.15] <- NA
        for(j in 1:3)
            gp[y, i, j] <- x[,8 + j]
    }

    ## rename
    megaG <- g
    megaPos <- pos
    megaSNPs <- snps
    megaSamples <- samples


    ##
    ## fix mm9/mm10 here
    ##
    mm10 <- read.table(file.path(megamugadir, "snps.megamuga.mm10.bed"))
    mm9 <- read.table(file.path(megamugadir, "snps.megamuga.mm9.bed"))
    un <- read.table(file.path(megamugadir, "snps.megamuga.unlifted.bed"))
    ## add to mm10 first
    mm9=cbind(mm9,-1,-1,-1,-1)
    mm9[match(un[,4],mm9[,4]),5:8]=NA
    for(i in 5:8)
        mm9[is.na(mm9[,5])==FALSE,i]=as.character(mm10[,i-4])

    ##
    megaSNPs=cbind(megaSNPs,NA,NA)
    colnames(megaSNPs)[ncol(megaSNPs)+-1:0]=c("Chromosome10","Position10")
    ##
    x=match(mm9[,4],megaSNPs[,"Name"])
    megaSNPs[x,c("Chromosome10","Position10")]=mm9[,c(5,7)]
    ## also - redo pos
    megaPos[,1]=megaSNPs[,"Chromosome10"]
    megaPos[,2]=megaSNPs[,"Position10"]

    nSNPsRemoved=nrow(mm9)-nrow(mm10)
    if (verbose)
        message(paste0(nSNPsRemoved," SNPs removed due to liftOver"))

    ##
    ## get information about the run quality
    ##
    ##GenCall Version = 6.3.0,Low GenCall Score Cutoff = 0.1500
    run <- read.table(
        file.path(megamugadir, "Univ_of_Chicago_St_Pierre_MEGMUGV01_20150505_DNAReport.csv"),
        sep=",",skip=2,header=TRUE,comment.char="@"
    )

    return(
        list(
            megaG = megaG,
            megaPos = megaPos,
            meagSNPs = megaSNPs,
            megaSamples = megaSamples
        )
    )

}


## take the list of megamuga objects
## convert to the "calls" format for affy
convert_megamuga_to_calls <- function(megamuga, chr, whose_samples) {

    if (verbose)
        message("Convert megamuga genotypes to calls")

    ## calls requires both chr and pos
    ## calls also has entries for SNPs of  -1, 0, 1, 2
    ## rownames chr-pos-ref-alt
    ## colnames of Q_*
    megaG <- megamuga$megaG
    megaPos <- megamuga$megaPos
    meagSNPs <- megamuga$megaSNPs
    megaSamples <- megamuga$megaSamples

    ## re-size and re-name megaG
    ## Q_CFW-SW___49.0e_recal
    ## subset both megaG and megaSamples to Oxford CFW samples
    if (whose_samples == "oxford") {
        samples <- megaSamples[grep(".", megaSamples[, "betterName"], fixed = TRUE), "Name"]
    } else if (whose_samples == "chicago") {
        samples <- megaSamples[-grep(".", megaSamples[, "betterName"], fixed = TRUE), "Name"]
        ## remove missing samples
        ##samples <- setdiff(samples, c("43919", "26658", "26351"))
    }
    megaG <- megaG[, match(samples, colnames(megaG))]
    megaSamples <- megaSamples[match(samples, megaSamples[, "Name"]), ]

    if (whose_samples == "oxford") {
        colnames(megaG) <- paste0("Q_CFW-SW___", megaSamples[, "betterName"], "_recal")
    }
    ## re-name mis-labelled samples
    if (whose_samples == "chicago") {
        ## re-name mis-labelled samples
        colnames(megaG)[match(c("42401", "42937"), colnames(megaG))] <- c("42937", "42401")
    }


    
    ## these files failed QC - see email from Sat Jun 13, 2015, 6:57 AM
    ## ~/proj/outbred/megaMugaQC_v0.0.0.R
    ## remove samples that are of poor quality
    ## from this failQC=run[,"X10._GC_Score"]<0.50
    remove1 <- c('42856', 'Q_CFW-SW___121.0b_recal', '42857', '43343', '42921', '42634', 'Q_CFW-SW___139.0g_recal', '42628', '45693', '46818', 'Q_CFW-SW___116.0m_recal', '46831', 'Q_CFW-SW___91.0d_recal', '45846', '46879', '43751', '26351', '47018') ## fail QC
    remove2 <- c('43919', '26658','26351') ## samples that do not exists
    remove <- union(remove1, remove2)
    to_remove <- match(remove, colnames(megaG))
    if (sum(!is.na(to_remove)) > 0) {
        megaG <- megaG[, -to_remove[!is.na(to_remove)]]
    }
    
    message(paste0("There are ", ncol(megaG), " samples available from the truth megamuga data"))

    ## re-size and add rownames
    keep <- megaPos[, "CHR"] == chr & is.na(megaPos[, "CHR"]) == FALSE
    megaG <- megaG[keep, ]
    megaPos <- megaPos[keep, ]
    ## remove duplicates
    keep <- match(unique(megaPos[, "POS"]), megaPos[, "POS"])
    megaG <- megaG[keep, ]
    megaPos <- megaPos[keep, ]
    temp <- colnames(megaG)
    megaG <- data.frame(
        chr = megaPos[, "CHR"],
        pos = megaPos[, "POS"],
        megaG
    )
    colnames(megaG)[-c(1:2)] <- temp
    rownames(megaG) <- NULL
    rownames(megaG) <- paste(
        megaPos[, "CHR"], megaPos[, "POS"], megaPos[, "REF"], megaPos[, "ALT"], sep = "-"
    )

    return(megaG)

}


## get all CFW samples for a chromosome
## return in a matrix with row = SNP, column = individual
get_affymetrix_calls_for_chr <- function(chr) {
    if (verbose)
        message("Get Affymetrix calls for chr")

    calls <- fread(
        paste0("grep -v '^#' ", file.path(affydir, "quant-norm.pm-only.brlmm-p.calls.txt")),
        data.table = FALSE,
        showProgress = FALSE
    )
    annot <- fread(
        paste0("grep -v '^#' ", file.path(affydir, "MOUSEDIVm520650.na33.annot.csv")),
        data.table = FALSE,
        showProgress = FALSE
    )
    annot <- annot[, c(1:4, 7:8)]
    colnames(annot) <- c("probeset_id", "rsid", "chr", "pos", "allele_A", "allele_B")
    annot[, "chr"] <- paste0("chr", annot[, "chr"])
    annot <- annot[annot[, "chr"] == chr, ]
    remove <- nchar(annot[, "allele_A"]) != 1 | nchar(annot[, "allele_B"]) != 1
    annot <- annot[remove == FALSE, ]
    calls <- calls[is.na(match(calls[, "probeset_id"], annot[, "probeset_id"])) == FALSE, ]
    both <- merge(annot, calls, by = "probeset_id")
    rownames(both) <- paste(both[, "chr"], both[, "pos"], both[, "allele_A"], both[, "allele_B"], sep = "-")

    info <- read.table(file.path(affydir, "Animal_info.txt"), header = FALSE)
    colnames(info) <- c("cel_name","name","cel_files","location","batch1","batch2","gender")
    report2 <- read.table(file.path(affydir, "quant-norm.pm-only.brlmm-p.report.txt"),header=TRUE)
    keep_samples <- is.na(match(report2[, "cel_files"], info[, "cel_files"])) == FALSE
    old_names <- info[, "cel_files"]
    new_names <- paste0(as.character(info[, "cel_name"]), "_recal")

    ## subset both to include necessary columns + desired columns
    keep_cols <- c(
        colnames(annot),
        colnames(both)[is.na(match(colnames(both), old_names)) == FALSE]
    )
    both <- both[, keep_cols]
    a <- match(colnames(both), old_names)
    colnames(both)[is.na(a) == FALSE] <- new_names[a[is.na(a) == FALSE]]

    return(both)
}



## remove SNPs that perform poorly
## not bothering with sites close to the array for now
filter_calls_for_chr <- function(calls) {
    if (verbose)
        message("Filter array calls for chr")
    arrayCounts <- cbind(
        rowSums(calls == 0, na.rm=TRUE),
        rowSums(calls == 1, na.rm=TRUE),
        rowSums(calls == 2, na.rm=TRUE)
    )
    arrayHWE <- apply(arrayCounts, 1, calculate_hwe_p)
    ##    N <- length(grep("Q_", colnames(calls)))
    N <- ncol(mega_calls) - 2 ## not ideal!
    missing <- rowSums(is.na(calls)) / N
    ## (distance < 12) |
    remove <-
        (missing > 0.05) |
        (arrayHWE < 1e-10)
    if (verbose)
        message(paste0("Removing ", sum(remove), " out of ", length(remove), " SNPs"))
    calls <- calls[remove == FALSE, ]
    return(calls)
}


## for a VCF
## return calls as appropriate
## possibly filter on subjects and SNPs
get_dosages_for_vcf <- function(
    vcf_file,
    chr,
    subjects = NULL,
    calls = NULL,
    use_fread = TRUE,
    use_chr = TRUE,
    tabix = TRUE,
    nCores = 1
) {

    if (verbose)
        message("Get dosages for VCF")

    if (tabix) {
        if (!file.exists(paste0(vcf_file, ".tbi"))) {
            system(paste0("tabix -f ", vcf_file))
        }
    }

    command <- paste0("bcftools  view ")

    if (is.null(subjects) == FALSE) {
        subjects_file <- tempfile()
        command <- paste0(command, "-S ", subjects_file, " ")
        if (length(subjects) == 0) {
            stop("0 length subject vector observed. Possible lack of intersection between truth and test data")
        }
        write.table(
            matrix(subjects, ncol = 1), row.names = FALSE, col.names = FALSE, quote = FALSE,
            file = subjects_file
        )
        col_names <- strsplit(system(paste0("bcftools view -h ", vcf_file, " -S ", subjects_file, "| grep '^#CHROM' | head -n1"), intern = TRUE), "\t")[[1]]
    } else {
        col_names <- strsplit(system(paste0("bcftools view -h ", vcf_file, " | grep '^#CHROM' | head -n1"), intern = TRUE), "\t")[[1]]
    }
    if (is.null(calls) == FALSE) {
        regions_file <- tempfile()
        command <- paste0(command, "-R ", regions_file, " ")
        regions <- calls[, c("chr", "pos")]
        write.table(
            regions, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t",
            file = regions_file
        )
    }

    command <- paste0(command, " -H ", vcf_file)    
    if (use_chr)
        command <- paste0(command, " ", chr)

    if (verbose)
        message("Extract and load VCF")
    if (use_fread) {
        vcf <- fread(cmd = command, data.table = FALSE, sep = "\t")
    }  else {
        vcf <- read.table(vcf_file)
    }
    colnames(vcf) <- col_names

    n_snps <- nrow(vcf)
    format <- strsplit(as.character(vcf[, "FORMAT"]), ":")
    gt_spot <- lapply(format, function(x) which(x == "GT"))
    gp_spot <- lapply(format, function(x) which(x == "GP"))
    ds_spot <- lapply(format, function(x) which(x == "DS"))

    if (verbose)
        message("Extract dosages from loaded VCF")
    dosages <- array(NA, c(nrow(vcf), ncol(vcf) - 9))
    colnames(dosages) <- colnames(vcf)[(10):ncol(vcf)]
    rownames(dosages) <- paste(vcf[, "#CHROM"], vcf[, "POS"], vcf[, "REF"], vcf[, "ALT"], sep = "-")
    genotypes <- dosages
    
    for(i_samp in 1:ncol(dosages)) {
        samp <- colnames(dosages)[i_samp]
        ## first, get DS if it exists
        a <- strsplit(as.character(vcf[, 9 + i_samp]), ":")
        for(i_snp in 1:nrow(vcf)) {
            b <- a[[i_snp]]
            x1 <- ds_spot[[i_snp]]
            x2 <- gp_spot[[i_snp]]
            x3 <- gt_spot[[i_snp]]
            if (length(x1) == 1) {
                d <- as.numeric(b[x1])
            } else if (length(x2) == 1) {
                c <- as.numeric(strsplit(b[x2], ",")[[1]])
                d <- sum(c[2] + c[3] * 2)
            } else if (length(x3) == 1) {
                d <- c(NA, 0, 1, 2)[
                    match(b[x3], c("./.", "0/0", "0/1", "1/1"))
                ]
            }
            dosages[i_snp, i_samp] <- d
            if (length(x2) == 1) {
                c <- as.numeric(strsplit(b[x2], ",")[[1]])
                genotypes[i_snp, i_samp] <- which.max(c) - 1
            }            
        }

    }

    if (verbose)
        message("Extract meta information from VCF")
    info <- strsplit(as.character(vcf[, "INFO"]), ":")
    dosages_meta <- t(sapply(info, function(x) {
        y <- strsplit(x, ";")[[1]]
        return(c(
            hwe = get_col_from_info(y, col = "HWE="),
            eaf = get_col_from_info(y, col = "EAF="),
            info = get_col_from_info(y, col = "INFO_SCORE=")
        ))
    }))
    rownames(dosages_meta) <- rownames(dosages)

    return(
        list(
            dosages = dosages,
            dosages_meta = dosages_meta,
            genotypes = genotypes
        )
    )

}

get_col_from_info <- function(y, col = "EAF=") {
    z <- substr(y, 1, nchar(col)) == col
    if (sum(z) == 1) {
        return(as.numeric(strsplit(y[z], "=")[[1]][2]))
    } else {
        return(NA)
    }
}


## for a VCF
## return calls as appropriate
## possibly filter on subjects and SNPs
get_dosages_for_bgen <- function(
    test_file,
    subjects,
    calls
) {

    if (verbose)
        message("Get dosages for bgen")

    ## pre-intersect
    var_info <- rrbgen_load_variant_info(test_file)
    var_info[, "varid"] <- gsub("-", ":", var_info[, "varid"])
    snps1 <- gsub("-", ":", rownames(calls))
    x <- t(sapply(strsplit(snps1, ":"), I))
    snps2 <- paste0(x[, 1], ":", x[, 2], ":", x[, 4], ":", x[, 3])
    vars_to_get <- intersect(var_info[, "varid"], c(snps1, snps2))

    ## if only at positions
    ## var_info <- rrbgen_load_variant_info(test_file)
    ## var12 <- paste(var_info[, "chr"], var_info[, "position"], sep = ":")
    ## m <- t(sapply(strsplit(rownames(mega_calls), "-"), I))
    ## calls12 <- paste(m[, 1], m[, 2], sep = ":")
    ## ##
    ## to_get <- intersect(var12, calls12)
    ## vars_to_get <- var_info[match(to_get, var12), "varid"]
    
    outX <- rrbgen_load(
        bgen_file = test_file,
        samples_to_get = subjects,
        vars_to_get = vars_to_get
    )

    ##dosages <- cbind(
    ##    chr = outX$var_info[, "chr"],
    ##    pos = outX$var_info[, "position"],
    ##    outX$gp[, , 2] + 2 * outX$gp[, , 3]
    ## )

    dosages <- outX$gp[, , 2] + 2 * outX$gp[, , 3]
    rownames(dosages) <- gsub(":", "-", rownames(dosages))
    
    per_snp_stats <- read.table(paste0(test_file, ".per_snp_stats.txt.gz"), header = TRUE)

    dosages_meta <- per_snp_stats[match(vars_to_get, per_snp_stats[, "varid"]), c("HWE", "EAF", "INFO_SCORE")]
    colnames(dosages_meta) <- c("hwe", "eaf", "info")

    ## out$dosages with row = SNP, col = sample, content is 0-2 dosage
    ## rows labelled with chr-snp-ref-alt
    ## out$dosages_meta = 3 cols with hwe, eaf, info

    return(
        list(
            dosages = dosages,
            dosages_meta = dosages_meta
        )
    )

}


get_col_from_info <- function(y, col = "EAF=") {
    z <- substr(y, 1, nchar(col)) == col
    if (sum(z) == 1) {
        return(as.numeric(strsplit(y[z], "=")[[1]][2]))
    } else {
        return(NA)
    }
}


compare_array_calls_and_dosages <- function(calls, dosages, dosages_meta, genotypes, return_things = FALSE) {
    
    if (verbose) {
        message("Compare array calls and dosages")
    }

    remove_cols <- c("chr", "pos", "probeset_id", "rsid", "allele_A", "allele_B")
    t <- match(remove_cols, colnames(calls))
    to_remove <- t[is.na(t) == FALSE]
    if (length(to_remove) > 0) {
        calls <- calls[, -to_remove]
    }
    if (ncol(calls) != ncol(dosages)) {
        stop("Bad removal of unnecessary columns")
    }
                     
    ## match on both alleles
    ## if match on flip, then flip dosages
    ## actually, only merge on first two items
    ## then reject if alleles mismatch
    c1A <- t(sapply(strsplit(rownames(calls), "-"), I))
    c1B <- c1A[, c(1, 2, 4, 3)]
    c1A <- apply(c1A, 1, paste, collapse = "-")
    c1B <- apply(c1B, 1, paste, collapse = "-")
    
    ## which ones required flips? flip now, then re-align
    ## I think I was just matching on position before?
    t1 <- match(c1B, rownames(dosages))
    to_flip_from_calls <- which(is.na(t1) == FALSE)
    if (length(to_flip_from_calls) > 0) {
        ## crap, different number of things in them?
        calls[to_flip_from_calls, ] <- 2 - calls[to_flip_from_calls, ]
        rownames(calls)[to_flip_from_calls] <- c1B[to_flip_from_calls]
    }

    joint <- intersect(rownames(calls), rownames(dosages))
    callsS <- calls[match(joint, rownames(calls)), ]
    genotypesS <- genotypes[match(joint, rownames(genotypes)), ]    
    dosagesS <- dosages[match(joint, rownames(dosages)), ]
    dosages_metaS <- dosages_meta[match(joint, rownames(dosages)), ]

    c1 <- colnames(callsS)
    c2 <- colnames(dosagesS)
    ##    subjects <- intersect(c1[grep("Q", c1)], c2[grep("Q", c2)])
    subjects <- intersect(c1, c2)
    callsS <- callsS[, match(subjects, c1)]
    dosagesS <- round(dosagesS[, match(subjects, c2)], 2) ## why are these rounded at all
    genotypesS <- genotypesS[, match(subjects, c2)]
    callsS[callsS == -1] <- NA

    r2 <- sapply(1:nrow(callsS), function(i_row) {
        cor(
            unlist(callsS[i_row, ]),
            dosagesS[i_row, ],
            use = "complete.obs"
        ) ** 2
    })
    pass_qc <- dosages_metaS[, "info"] > 0.4 & dosages_metaS[, "hwe"] > 1e-6

    message(paste0(sum(pass_qc), " / ", length(pass_qc), " SNPs pass QC"))
    message(
        paste0(
            "Average r2 of all SNPs:",
            round(mean(r2, na.rm = TRUE), 3)
        )
    )
    message(
        paste0(
            "Average r2 of SNPs that pass QC:",
            round(mean(r2[pass_qc], na.rm = TRUE), 3)
        )
    )
    message(
        paste0(
            "Average r2 of SNPs that fail QC:",
            round(mean(r2[pass_qc == FALSE], na.rm = TRUE), 3)
        )
    )

    maf <- dosages_metaS[, "eaf"]
    maf[maf > 0.5] <- 1 - maf[maf > 0.5]

    message("r2 for various MAF ranges for SNPs that pass QC")
    print(tapply(r2[pass_qc], cut(maf[pass_qc], c(0, 0.01, 0.05, 0.1, 0.2, 0.5)), mean, na.rm=TRUE))

    r <- c(
        agree = sum(genotypesS == callsS, na.rm = TRUE),
        disagree = sum(genotypesS != callsS, na.rm = TRUE),
        missing = sum(is.na(genotypesS == callsS))
    )
    conc <- round(r["agree"] / (r["agree"] + r["disagree"]), 3)
    message(paste0(100 * conc, "% concordance among non-missing sites"))
    miss <- round(r["missing"] / sum(r), 3)
    message(paste0(100 * miss, "% concordance among non-missing sites"))        

    if (return_things) {
        return(list(callsS = callsS, dosagesS = dosagesS, maf = maf, pass_qc = pass_qc, dosages_metaS = dosages_metaS))
    } else {
        return(NULL)
    }

}


