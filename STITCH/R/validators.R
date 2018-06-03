validate_output_format <- function(output_format) {
    if (class(output_format) != "character") {
        stop(paste0("output_format must have class character but you have selected:", class(output_format)))
    }
    if ((output_format != "bgen") & (output_format != "gvcf")) {
        stop(paste0("output_format must be either gvcf of bgen"))
    }
    return(NULL)
}


validate_output_filename <- function(
    output_filename,
    output_format
) {
    if (output_format == "gvcf") {
        extension <- ".vcf.gz"
    } else if (output_format == "bgen") {
        extension <- ".bgen"
    } else {
        stop("internal error")
    }
    min_chars <- (nchar(extension))        
    if (is.null(output_filename) == FALSE) {
        err_msg <- paste0("output_filename must have at least ", min_chars + 1, " characters and end with ", extension, ", and you have supplied:", output_filename)
        if (nchar(output_filename) < min_chars) {
            stop(err_msg)
        }
        if (substr(
            output_filename,
            nchar(output_filename) - (min_chars - 1),
            nchar(output_filename)
        ) != extension) {
            stop(err_msg)
        }
    }
    return(NULL)
}


validate_pos_and_legend_snps_for_niterations_equals_1 <- function(legend_snps, pos_snps, niterations) {
    if (niterations == 1) {
        x <- is.na(match(legend_snps, pos_snps))
        if (sum(x) > 0)
            stop(paste0("You have selected to use reference haplotypes with niterations=1, which requires exact matching of reference legend SNPs and posfile SNPs. However, reference legend SNP with pos-ref-alt ", legend_snps[which.max(x)], " was not found in posfile"))
        x <- is.na(match(pos_snps, legend_snps))
        if (sum(x) > 0)
            stop(paste0("You have selected to use reference haplotypes with niterations=1, which requires exact matching of reference legend SNPs and posfile SNPs. However, posfile SNP with pos-ref-alt ", pos_snps[which.max(x)], " was not found in reference legend"))
    }
    return(NULL)
}


print_and_validate_reference_snp_stats <- function(
    pos_snps,
    legend_snps,
    both_snps
) {
    ## print out stats
    print_message(paste0("In the region to be imputed plus buffer there are the following number of SNPs:"))
    print_message(paste0(length(pos_snps), " from the posfile"))
    print_message(paste0(length(legend_snps), " from the reference haplotypes"))
    print_message(paste0(length(both_snps), " in the intersection of the posfile and the reference haplotypes"))
    if (length(both_snps) == 0)
        stop("There are 0 SNPs that intersect between the posfile and the reference haplotypes. Please troubleshoot to see if there is an error, or re-run without reference haplotypes")
}


validate_reference_sample_file <- function(samples, reference_sample_file) {
    for(col in c("POP", "SEX")) {
        if (is.na(match(col, colnames(samples))))
            stop(
                paste0(
                    "Cannot find column ", col, " in reference_sample_file:",
                    reference_sample_file
                )
            )
    }
    x <- samples[, "SEX"] != "male" & samples[, "SEX"] != "female"
    if (sum(x) > 0)
        stop(
            "Reference samples must have SEX male or female ",
            "and in row ", which.max(x), " there is entry ",
            samples[which.max(x), "SEX"], " from file:", reference_sample_file
        )
}



validate_phase <- function(phase) {
    if (length(dim(phase)) != 3)
        stop("The phasefile does not have the right number of dimensions")
    if (sum(phase[, , 1] != 0 & phase[, , 1] != 1 & phase[, , 2] != 0 & phase[, , 2] != 1) > 0) {
        for(i_col in 1:dim(phase)[2]) {
            p <- phase[, i_col, ]
            i_row <- which.max(rowSums(p != 0 & p != 1) > 0)
            p2 <- phase[i_row, i_col, ]
            stop(paste0(
                "The phasefile contains entries other than 0 or 1. ",
                "One such entry is in column ", i_col, " and row ", i_row, " ",
                " with value ", paste(p2, collapse = "|")
            ))
        }
    }
}


## validate pos matrix shortly after loading
validate_pos <- function(
    pos,
    chr,
    stop_file_name = "pos file"
) {
    if (sum(pos[,1] != chr) > 0) {
        m <- which.max(is.na(match(pos[,1], chr)))
        stop(paste0(
            stop_file_name, " needs to be unique to chromosome ",
            "and you supplied chr=", chr, " but ",
            "chromosome ", pos[m, 1], " was observed in row ", m
        ))
    }
    x <- suppressWarnings(is.na(as.integer(as.character(pos[, 2]))))
    if (sum(x) > 0) {
        stop(paste0(
            stop_file_name, " column 2 needs to be integer valued ",
            "between 1 and ", .Machine$integer.max, " but ",
            "in row ", which.max(x), " ", pos[which.max(x), 2],
            " was observed"
        ))
    }
    if (sum(diff(pos[, 2]) <= 0) > 0) {
        w <- which.max(diff(pos[, 2]) <= 0)
        stop(paste0(
            stop_file_name, " column 2 needs to be sorted on position ",
            "with increasing positions between rows ",
            "but row number ", w, " has position ", pos[w, 2], " ",
            "and row number ", w + 1, " has position ", pos[w + 1, 2]
        ))
    }
    for(col in 3:4) {
        x <- is.na(match(pos[, col], c("A", "C", "G", "T")))
        if (sum(x) > 0) {
            m <- which.max(x)
            y <- c("", "", "ref", "alt")[col]
            stop(paste0(
                stop_file_name, " ", y, " column entry ", pos[m, col],
                " in row ", m, " ",
                "contains is not one or A, C, G or T. STITCH is ",
                "only supported for bi-allelic SNPs"
            ))
        }
    }
    return(NULL)
}



validate_vcf_output_name <- function(
    vcf_output_name
) {
    if (is.null(vcf_output_name) == FALSE) {
        if (nchar(vcf_output_name) < 8)
            stop(paste0("vcf_output_name must have at least 8 characters and end with .vcf.gz, and you have supplied vcf_output_name:", vcf_output_name))
        if (substr(vcf_output_name, nchar(vcf_output_name) - 6, nchar(vcf_output_name)) != ".vcf.gz")
            stop(paste0("vcf_output_name must end with .vcf.gz, and you have supplied vcf_output_name:", vcf_output_name))
    }
    return(NULL)
}




validate_phase_header <- function(phasefile) {
    first_row <- as.character(unlist(read.table(phasefile,sep="\t",nrows=1)))
    ## check none are 0|0, 0|1, 1|0, 1|1
    m <- match(first_row, c("0|0", "0|1", "1|0", "1|1"))
    if (sum(is.na(m) == FALSE) > 0) {
        m2 <- which.max(is.na(m) == FALSE)
        stop(paste0(
            "The header for the phasefile is either invalid or missing. ",
            "The first invalid entry is in position ", m2, " and is entry ", first_row[m2]
        ))
    }
}

validate_phase_col <- function(col, i_samp) {
    x <- sapply(col, length)
    if (sum(x != 2) > 0) {
        m <- which.max(x)
        stop(paste0(
            "Unable to split column ", i_samp, " of phasefile at position ", m,
            " with entry ", col[m], " due to lack of field separator |"
        ))
    }
}




get_and_validate_pos_gen_and_phase <- function(
    posfile,
    genfile = "",
    phasefile = "",
    chr,
    verbose = FALSE
) {
    if (verbose)
        print_message("Get and validate pos and gen")
    pos <- get_and_validate_pos(posfile, chr)
    gen <- get_and_validate_gen(genfile)
    phase <- get_and_validate_phase(phasefile)
    if (is.null(gen) == FALSE)
        if(nrow(gen) != nrow(pos))
            stop(paste0(
                "posfile and genfile must have the same number of SNPs ",
                "but posfile has ", nrow(pos), " SNPs and genfile has ",
                nrow(gen), " SNPs (and 1 header row)"
            ))
    if (is.null(phase) == FALSE)
        if(dim(phase)[1] != nrow(pos))
            stop(paste0(
                "posfile and phasefile must have the same number of SNPs ",
                "but posfile has ", nrow(pos), " SNPs and phasefile has ",
                dim(phase)[1], " SNPs (and 1 header row)"
            ))
    L <- as.integer(as.character(pos[, "POS"]))
    nSNPs <- as.integer(length(L))
    if (verbose)
        print_message("Done get and validate pos and gen")
    return(
        list(
            pos = pos,
            gen = gen,
            phase = phase,
            nSNPs = nSNPs,
            L = L
        )
    )
}

get_and_validate_gen <- function(genfile) {
    if (genfile == "")
        return(NULL)
    gen2 <- read.table(genfile, header=TRUE, sep="\t")
    gen <- array(0, c(nrow(gen2), ncol(gen2)))
    for(i_col in 1:ncol(gen)) {
        col <- gen2[, i_col]
        validate_gen_col(col, i_col)
        gen[, i_col] <- as.integer(as.character(gen2[, i_col]))
    }
    header <- as.character(unlist(read.table(genfile, sep="\t", nrows = 1)))
    validate_gen_header(header)
    colnames(gen) <- header
    return(gen)
}


validate_gen_header <- function(header) {
    ## only allowed are 0, 1, 2, NA
    m <- match(header, c(0, 1, 2, NA))
    if (sum(is.na(m) == FALSE) > 0) {
        m2 <- which.max(is.na(m) == FALSE)
        stop(paste0(
            "The header for the genfile is either invalid or missing. ",
            "The first invalid entry is in position ", m2, " and is entry ", header[m2]
        ))
    }
}

## can only be 0, 1, 2, NA
validate_gen_col <- function(col, i_col) {
    m <- match(col, c(0, 1, 2, NA))
    if (sum(is.na(m)) > 0) {
        m2 <- which.max(is.na(m))
        stop(paste0(
            "Ineligible entry for genfile column ", i_col, " at position ",
            m2, " with entry ", col[m2], ". Acceptable entries are 0, 1, 2 or NA"
        ))
    }
}

get_and_validate_pos <- function(posfile, chr) {
    pos <- read.table(posfile, header = FALSE, sep = "\t")
    colnames(pos) <- c("CHR", "POS", "REF", "ALT")
    validate_pos(pos, chr)
    return(pos)
}

get_and_validate_phase <- function(
    phasefile
) {
    if (phasefile == "")
        return(NULL)
    phaseX <- read.table(phasefile, header = TRUE, stringsAsFactors = FALSE)
    validate_phase_header(phasefile)
    phase <- array(0, c(nrow(phaseX), ncol(phaseX), 2))
    for(i_samp in 1:ncol(phase)) {
        col <- strsplit(phaseX[, i_samp], "|", fixed = TRUE)
        validate_phase_col(col, i_samp)
        for(j in 1:2)
            phase[, i_samp, j] <- as.integer(sapply(col, function(x) x[j]))
    }
    if (length(colnames(phaseX)) == 1) {
        dimnames(phase)[[2]] <- list(colnames(phaseX))
    } else {
        dimnames(phase)[[2]] <- colnames(phaseX)
    }
    validate_phase(phase)
    return(phase)
}





validate_gridWindowSize <- function(gridWindowSize) {
    if (is.na(gridWindowSize)) {
        return(NULL)
    } else if (is.integer(gridWindowSize) | is.numeric(gridWindowSize)) {
        if (gridWindowSize < 2) {
            stop(paste0("gridWindowSize must be an integer > 1 but you have selected:", gridWindowSize))
        } else if (round(gridWindowSize) != gridWindowSize) {
            stop(paste0("gridWindowSize must be an integer but you have selected:", gridWindowSize))
        }  else {
            return(NULL)
        }
    } else {
        stop(paste0("gridWindowSize must be an integer or numeric or null but you have selected an object of class ", class(gridWindowSize), " with entry:", gridWindowSize))
    }
}


validate_chr <- function(chr) {
    if (length(chr) == 0)
        stop("Please specify chr, the chromosome to impute")
    if(chr=="")
        stop("Please specify chr, the chromosome to impute")
}


validate_nCores <- function(nCores) {
    if (is.numeric(nCores) == FALSE) {
        stop("nCores must be an integer")
    }
    if (nCores < 1)
        stop("nCores must be greater than or equal to 1")
    if (round(nCores) != nCores)
        stop("nCores must be an integer")
}

validate_nGen <- function(nGen) {
    nGen_error <- paste0(
        "Please specify nGen, the estimate of the number of generations ",
        "since founding. Can be approximated by 4 * N_e / K, where",
        "N_e ~= 20,000 for humans"
    )
    if (is.na(nGen))
        stop(nGen_error)
    if( nGen == "")
        stop(nGen_error)
    if (is.numeric(nGen) == FALSE)
        stop("nGen needs to be a numeric")
    if (nGen < 0)
        stop("nGen must be greater than 0")
}

validate_posfile <- function(posfile) {
    if(posfile=="")
        stop("Please specify posfile, the file with sites to impute over")
    if (file.exists(posfile) == FALSE)
        stop(paste0("Cannot find supplied posfile:", posfile))
}

validate_genfile <- function(genfile) {
    if(genfile!="")
        if (file.exists(genfile) == FALSE)
            stop(paste0("Cannot find supplied genfile:", genfile))
}

validate_K <- function(K) {
    if(K=="") {
        stop("Pleasy specify K")
    }
    if (is.numeric(K) == FALSE) {
        stop(paste0("K must be numeric but class(K)=", class(K)))
    }
    if (round(K) != K) {
        stop(paste0("K must be an integer but you have selected K=", K))
    }
    return(NULL)
}

validate_K_subset <- function(method, K, K_subset) {
    if (method == "diploid_subset") {
        if (is.na(K_subset)) {
            stop(paste0("K_subset must be specified if method=diploid_subset"))
        }
        if (is.na(K)) {
            stop(paste0("K must not be NA"))
        }
        if (is.numeric(K_subset) == FALSE) {
            stop(paste0("K_subset must be an integer but you have selected class(K_subset)=", class(K_subset)))
        }
        if (round(K_subset) != K_subset) {
            stop(paste0("K_subset must be an integer but you have selected K_subset=", K_subset))
        }
        if ((K_subset %% 2) != 0) {
            stop(paste0("K_subset must an even integer but you have selected K_subset=", K_subset))
        }
        if (K_subset > K) {
            stop(paste0("K_subset must an even integer less than K but you have selected K_subset=", K_subset, ", and K=", K))
        }
    }
    return(NULL)
}

validate_outputdir <- function(outputdir)
    if(outputdir=="")
        stop("Please specify outputdir")

validate_tempdir <- function(tempdir) {
    if (is.na(tempdir) == FALSE) {
        if (tempdir=="")
            stop(
                paste0(
                    "Please specify tempdir. If in doubt, use tempdir = tempdir()"
                )
            )
    }
}

validate_outputBlockSize <- function(outputBlockSize) {
    if(outputBlockSize < 1)
        stop("outputBlockSize must be greater than 1")
    if(round(outputBlockSize) != outputBlockSize)
        stop("outputBlockSize must be an integer")
    return(NULL)
}


validate_downsampleFraction <- function(downsampleFraction) {
    if(downsampleFraction < 0 | downsampleFraction > 1)
        stop("downsampleFraction must be between 0 and 1")
}

validate_downsampleSamples <- function(downsampleSamples) {
    if (is.numeric(downsampleSamples) == FALSE)
        stop("downsampleSamples must be a integer > 0 and <=1")
    if (downsampleSamples <= 0 | downsampleSamples >1)
        stop("downsampleSamples must be a integer > 0 and <=1")
    return(NULL)
}

validate_plotHapSumDuringIterations <- function(plotHapSumDuringIterations)
    if (is.logical(plotHapSumDuringIterations) == FALSE | is.na(plotHapSumDuringIterations))
        stop("plotHapSumDuringIterations must be either TRUE or FALSE")


validate_regionStart_regionEnd_and_buffer <- function(regionStart, regionEnd, buffer) {
    w <- as.integer(is.na(regionStart))+
        as.integer(is.na(regionEnd))+
        as.integer(is.na(buffer))
    if ( w !=0 & w!=3)
        stop("either all three of regionStart, regionEnd and buffer are NA, or none of them are")
    if (is.na(regionStart)==FALSE) {
        if (regionStart<1) stop("regionStart must be >0")
        if (round(regionStart)!=regionStart) stop("regionEnd must be an integer")
    }
    if (is.na(regionEnd)==FALSE) {
        if (regionEnd<1) stop("regionEnd must be >0")
        if (round(regionEnd)!=regionEnd) stop("regionEnd must be an integer")
    }
    if (is.na(buffer)==FALSE) {
        if (buffer < 0)
            stop("buffer must be >=0")
        if (round(buffer)!=buffer) stop("buffer must be an integer")
    }
    if (is.na(regionStart)==FALSE & is.na(regionEnd)==FALSE) {
        if (regionStart>regionEnd) stop("regionStart must be before regionEnd")
    }
}

validate_bamlist_and_cramlist_for_input_generation <- function(
    regenerateInput = TRUE,
    originalRegionName = NA,
    bamlist = "",
    cramlist = "",
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    regenerateInputWithDefaultValues = FALSE,
    reference = NULL
) {
    if (regenerateInput == FALSE) {
        if (is.na(originalRegionName) == TRUE)
            stop("if regenerateInput is FALSE (i.e. using existing data), you must supply the original region name (must not be NA) to load the input properly. Also don't forget to supply the position file used to make the original input data, as well as to use the same regionStart, regionEnd and buffer values")
        if (bamlist != "" || cramlist != "")
            stop("If not regenerating input, please do not supply bamlist")
        if ((is.na(regionStart) | is.na(regionEnd) | is.na(buffer)) & regenerateInputWithDefaultValues == FALSE)
            stop("If regenerateInput is FALSE, please supply original regionStart, regionEnd and buffer used to make the inputs. This is used to subset the pos file appropriately. If you used the default values of regionStart, regionEnd and buffer (NA), please set regenerateInputWithDefaultValues to TRUE")
    } else {
        if (bamlist == "" & cramlist == "")
            stop("If regenerateInput is TRUE, please supply either bamlist or cramlist")
        if (cramlist == "") {
            if (file.exists(bamlist) == FALSE) {
                stop(paste0("Cannot find bamlist:", bamlist))
            }
        }
        if (bamlist == "") {
            if (file.exists(cramlist) == FALSE) {
                stop(paste0("Cannot find cramlist:", cramlist))
            }
        }
    }
    if (cramlist != "") {
        if (bamlist != "")
            stop("Please specify only one or cramlist or bamlist")
        if (reference == "")
            stop("If you supply a list of CRAM files, you need to supply the reference they were aligned against")
    }
}


validate_inputBundleBlockSize <- function(inputBundleBlockSize, readAware) {
    if (is.na(inputBundleBlockSize) == FALSE) {
        if (readAware == FALSE)
            stop("readAware is FALSE is not currently supported when bundling inputs together")
        if (inputBundleBlockSize < 1)
            stop("inputBundleBlockSize must be an integer greater than 0")
        if (round(inputBundleBlockSize) != inputBundleBlockSize)
            stop("inputBundleBlockSize must be an integer")
    }
}




validate_reference_files <- function(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations) {
    if (reference_haplotype_file == "" &
        reference_legend_file != "")
        stop("If you specify reference_legend_file you must also specify reference_haplotype_file")
    if (reference_haplotype_file != "" &
        reference_legend_file == "")
        stop("If you specify reference_haplotype_file you must also specify reference_legend_file")
    if (reference_sample_file != "" & is.na(reference_populations[1]) == TRUE)
        stop("The reference sample file has been provided but you are not specifying any populations to retain. If you want to use all reference haplotypes, please omit this file")
    if (reference_sample_file != "" & file.exists(reference_sample_file) == FALSE)
        stop(paste0("Cannot find reference_sample_file:", reference_sample_file))
    if (reference_legend_file != "" & file.exists(reference_legend_file) == FALSE)
        stop(paste0("Cannot find reference_legend_file:", reference_legend_file))
    if (reference_haplotype_file != "" & file.exists(reference_haplotype_file) == FALSE)
        stop(paste0("Cannot find reference_haplotype_file:", reference_haplotype_file))
}


validate_refillIterations <- function(refillIterations, niterations) {
    if (is.na(refillIterations[1]) == FALSE) {
        x <- niterations - refillIterations
        if (sum (x < 5) > 0)
            stop(paste0("The parameter niterations is set to ", niterations, " and refillIterations is set to c(", paste(refillIterations, collapse = ","), "). refillIterations (trying to refill infrequently used ancestral haplotypes) has a small but beneficial improvement on imputation performance, but due to methods involve, will initially lead to a decrease in imputation performance for the next few EM iterations. Please set niterations > (max(refillIterations) + 4), or if you want to use a small number of iterations, it is recommended you disable refillIterations by setting this parameter equal to NA"))
    }
}

validate_shuffleHaplotypeIterations <- function(shuffleHaplotypeIterations, niterations) {
    if (is.na(shuffleHaplotypeIterations[1]) == FALSE) {
        x <- niterations - shuffleHaplotypeIterations
        if (sum (x < 5) > 0)
            stop(paste0("The parameter niterations is set to ", niterations, " and shuffleHaplotypeIterations is set to c(", paste(shuffleHaplotypeIterations, collapse = ","), "). shuffleHaplotypeIterations (trying to reshuffle the ancestral haplotypes to minimize switching among ancestral haplotypes in analyzed samples) has a small but beneficial improvement on imputation performance, but due to methods involve, will initially lead to a decrease in imputation performance for the next few EM iterations. Please set niterations > (max(shuffleHaplotypeIterations) + 4), or if you want to use a small number of iterations, it is recommended you disable shuffleHaplotypeIterations by setting this parameter equal to NA"))
    }
}


validate_hapProb <- function(initial_min_hapProb, initial_max_hapProb) {
    if (is.numeric(initial_min_hapProb) == FALSE)
        stop("initial_min_hapProb must be numeric")
    if (is.numeric(initial_max_hapProb) == FALSE)
        stop("initial_max_hapProb must be numeric")
    if (initial_min_hapProb < 0 | 1 < initial_min_hapProb)
        stop("initial_min_hapProb must be between 0 and 1")
    if (initial_max_hapProb < 0 | 1 < initial_max_hapProb)
        stop("initial_max_hapProb must be between 0 and 1")
}
