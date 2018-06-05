## functions to do with writing results to disk
## either VCF or bgen


## if want to output genotype probabilities in VCF format
outputInputInVCFFunction <- function(
    outputdir,
    pos,
    nSNPs,
    tempdir,
    N,
    nCores,
    regionName,
    environment,
    sampleNames,
    outputBlockSize,
    bundling_info,
    output_filename,
    vcf.piece_unique,
    output_format
) {
    ##
    ##
    ## 1 - build the data
    ##
    ##
    x3 <- getSampleRange(N, nCores)
    f <- function(sampleRange, N, outputBlockSize, tempdir, regionName, nSNPs, bundling_info, vcf.piece_unique) {
        i_core <- match(sampleRange[1], sapply(x3, function(x) x[[1]]))
        outputBlockRange <- getOutputBlockRange(sampleRange,outputBlockSize )
        bundledSampleReads <- NULL
        outputBlockRange <- getOutputBlockRange(
            sampleRange,
            outputBlockSize
        )
        vcf_matrix_to_out <- data.frame(array(
            NA,
            c(nSNPs, outputBlockRange[2] - outputBlockRange[1])
        ))
        vcf_matrix_to_out_offset <- outputBlockRange[1]
        for (iSample in sampleRange[1]:sampleRange[2]) {
            out <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info,
                bundledSampleReads = bundledSampleReads
            )
            sampleReads <- out$sampleReads
            bundledSampleReads <- out$bundledSampleReads
            ##
            o <- get_pl_and_rd(sampleReads, nSNPs)
            pl <- o$pl
            rd <- o$rd
            ## start saving
            vcf_matrix_to_out[
               ,
                iSample - vcf_matrix_to_out_offset
            ] <- make_column_of_vcf_from_pl_and_rd(pl, rd)
            ## if appropriate - write to disk!
            iBlock <- match(iSample, outputBlockRange)
            if (iBlock > 1 & is.na(iBlock)==FALSE) {
                write_block_of_vcf(
                    i_core = i_core,
                    iBlock = iBlock,
                    vcf_matrix_to_out = vcf_matrix_to_out,
                    outputdir = outputdir,
                    regionName = regionName,
                    outputBlockRange = outputBlockRange,
                    vcf.piece_unique = vcf.piece_unique
                )
                ## initialize new matrix if not last one
                if (iBlock < length(outputBlockRange)) {
                    vcf_matrix_to_out <- data.frame(array(
                        NA,
                        c(nSNPs, outputBlockRange[iBlock] - outputBlockRange[iBlock - 1])
                    ))
                    vcf_matrix_to_out_offset <- outputBlockRange[iBlock - 1]
                }
            }
        }
    }
    ##
    ## loop over samples - write to disk!
    ##
    print_message("Prepare data to use to build vcf from input")
    if (environment == "server") {
        out2 <- mclapply(x3,mc.cores=nCores,FUN=f, N = N, outputBlockSize = outputBlockSize, tempdir = tempdir, regionName = regionName, nSNPs = nSNPs, bundling_info = bundling_info, vcf.piece_unique = vcf.piece_unique)
    }
    if (environment == "cluster") {
        cl <- makeCluster(nCores, type = "FORK")
        out2 <- parLapply(cl, x3, fun=f,N = N, outputBlockSize = outputBlockSize, tempdir = tempdir, regionName = regionName, nSNPs = nSNPs, bundling_info = bundling_info, vcf.piece_unique = vcf.piece_unique)
        stopCluster(cl)
    }
    error_check <- sapply(out2, class) == "try-error"
    if (sum(error_check) > 0) {
        print_message(out2[[which(error_check)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }
    ##
    ##
    ## step 2 - build the VCF itself
    ##
    ##
    ##
    ## build header
    ## build left side
    ## stitch everything together
    print_message("Build VCF from input")
    output_vcf <- get_output_filename(
        output_filename = output_filename,
        outputdir = outputdir,
        regionName = regionName,
        output_format = output_format,
        prefix = "stitch.input"
    )
    output_vcf_header <- paste0(output_vcf, ".header.gz")
    output_vcf_left <- paste0(output_vcf, ".left.gz")
    ##
    ## do header lines
    ##
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n'
    )
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    ##
    ## add output
    ##
    cat(header, header2, "\n", sep="", file = output_vcf_header)
    ##
    ## build first part of file
    ##
    options(scipen=6)
    INFO <- "."
    write.table(
        matrix(paste(pos[,1], pos[,2], ".", pos[,3], pos[,4], ".", "PASS", INFO , "GT:AD:PL", sep="\t"), ncol=1),
        file = gzfile(output_vcf_left),
        row.names = FALSE,
        col.names = FALSE,
        sep = "",
        quote = FALSE
    )
    ##
    ## get number of cores and blocks
    ##
    files_to_paste <- get_files_to_paste(
        N, nCores, outputBlockSize, outputdir, vcf.piece_unique, regionName
    )
    unlisted_files_to_paste <- paste(unlist(files_to_paste), collapse = " ")    
    ##
    ## merge together
    ##
    command <- paste0(
        'bash -c "',
        'paste -d ',shQuote("\t"),' ',
        '<(gunzip -c ', shQuote(output_vcf_left), ' | cat ) ',
        unlisted_files_to_paste, ' | ',
        ' cat ', shQuote(output_vcf_header) ,' - | bgzip -c > ',
        shQuote(output_vcf)
       ,'"'
    )
    system(command)
    unlink(output_vcf_header)
    unlink(output_vcf_left)
    ##
    ## remove
    ##
    files_to_remove <- get_files_to_paste(
        N, nCores, outputBlockSize, outputdir, vcf.piece_unique, regionName, as_command = FALSE
    )
    for(file in unlist(files_to_remove)) {
        unlink(file)
    }
    print_message("Done building VCF from input")
    return(NULL)
}



## write a block of samples to disk
## don't write header or row information
write_block_of_vcf <- function(
    i_core,
    iBlock,
    vcf_matrix_to_out,
    outputdir,
    regionName,
    outputBlockRange,
    vcf.piece_unique
) {
    print_message(paste0(
        "Write block of VCF for samples:",
        (outputBlockRange[iBlock - 1] + 1), "-", outputBlockRange[iBlock]
    ))
    file <- file.path(
        outputdir,
        paste0(
            "vcf.piece",
            vcf.piece_unique,
            i_core, ".", iBlock, ".", regionName, ".txt"
        )
    )
    data.table::fwrite(
        vcf_matrix_to_out,
        file = file,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    system(paste0("gzip -1 '", file, "'"))
    return(NULL)
}



# write out a single samples worth of VCF entry
make_column_of_vcf_from_pl_and_rd <- function(
  pl,
  rd
) {
  # write out genotype, genotype likelihood, and dosage
  # gt, ad, pl
  # GT:AD:PL
  #chr19   3126522 .       G       A       7661.22 .       AC=464;AF=0.935;AN=496;BaseQRankSum=4.232;DP=333;Dels=0.00;FS=2.773;GC=39.90;HRun=0;HaplotypeScore=0.0738;InbreedingCoeff=-0.0083;MLEAC=465;MLEAF=0.938;MQ=42.81;MQ0=23;MQRankSum=-0.576;QD=28.69;ReadPosRankSum=-0.673      GT:AD:DP:GQ:PL       ./.     ./.     ./.     1/1:0,1:1:3:32,3,0      ./.     ./.
  str <- paste0(
    "./.", ":",
    rd[,1], ",",
    rd[,2], ":" ,
    pl[,1], ",",
    pl[,2], "," ,
    pl[,3]
  )
  str[(rd[,1] + rd[,2])==0] = "./."
  return(str)
}




## now works on transpose, i.e. 3 x nSNPs, find max id per-col
get_max_gen_rapid <- function(x) {
    ## assume matrix 3 columns >1 row
    z <- rep(1, ncol(x))
    y <- x[1, ]
    for(i in 2:3) {
        w <- x[i, ] > y
        z[w] <- i
        y[w] <- x[i, w]
    }
    return(cbind(z, 1:ncol(x)))
}



## write out a single samples worth of VCF entry
## note - C++ version does not have "read_proportions" functionality
## I think this was for the phasing things
make_column_of_vcf <- function(
    gp_t,
    read_proportions = NULL
) {
    ## write out genotype, genotype likelihood, and dosage
    ##GT:GL:DS
    ##FORMAT=<ID=GT:,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
    ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Posterior probability of 0/0, 0/1, and 1/1">
    ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
    ## 1/1:0,0.054,0.946:1.946
    ## add one samples worth of info to a VCF
    z <- get_max_gen_rapid(gp_t)
    gt <- c("0/0","0/1","1/1")[z[, 1]]
    gt[gp_t[z] < 0.9] <- "./."
    precision <- 3
    format_string <- paste0(":%.", precision, "f,%.", precision, "f,%.", precision, "f:%.", precision, "f")
    str <- paste0(
        gt, 
        sprintf(format_string, gp_t[1, ], gp_t[2, ], gp_t[3, ], gp_t[2, ] + 2 * gp_t[3, ])
    )
    if (is.null(read_proportions) == FALSE) {
        format_string <- paste0(":%.", precision, "f,%.", precision, "f,%.", precision, "f,%.", precision, "f")
        str <- paste0(
            str, 
            sprintf(
                format_string,
                read_proportions[, 1],
                read_proportions[, 2],
                read_proportions[, 3],
                read_proportions[, 4]
            )
        )
    }
    return(str)
}





## specify here the output file from the process
get_output_filename <- function(
    output_filename,
    outputdir,
    regionName,
    output_format,
    prefix = "stitch"
) {
    if (output_format == "gvcf") {
        extension <- ".vcf.gz"
    } else {
        extension <- "bgen"
    }
    if (is.null(output_filename)) {
        return(
            file.path(
                outputdir,
                paste0(prefix, ".", regionName, extension)
            )
        )
    } else {
        if (basename(output_filename) == output_filename) {
            return(file.path(
                outputdir,
                output_filename
            ))
        } else {
            return(output_filename)
        }
    }
}


## build the header
## build the left part of the VCF
## paste together the interim VCF files
## gzip the whole thing together
write_vcf_after_EM <- function(
    output_filename,
    outputdir,
    regionName,
    sampleNames,
    tempdir,
    nCores,
    info,
    hwe,
    estimatedAlleleFrequency,
    pos,
    N,
    outputBlockSize,
    reference_panel_SNPs,
    method,
    output_format,
    vcf.piece_unique = "."
) {
    ## set up file names
    print_message("Build final VCF")
    output_vcf <- get_output_filename(
        output_filename = output_filename,
        outputdir = outputdir,
        regionName = regionName,
        output_format = output_format
    )
    output_vcf_header <- paste0(output_vcf, ".header.gz")
    output_vcf_left <- paste0(output_vcf, ".left.gz")
    ##
    ## do header lines
    ##
    header_line_if_ref_panel_snps <- NULL
    if (sum(reference_panel_SNPs) > 0)
        header_line_if_ref_panel_snps <- '##INFO=<ID=REF_PANEL,Number=.,Type=Integer,Description="Whether a SNP was (1) or was not (0) found in the reference panel during imputation">\n'
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score">\n',
        '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
        '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">\n',
        header_line_if_ref_panel_snps,
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">\n',
        '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">\n'
    )
    if (method == "pseudoHaploid" && 1 == 0) ## disable for now
        header <- paste0(
            header,
            '##FORMAT=<ID=ER1,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 1">\n',
            '##FORMAT=<ID=EA1,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 1">\n',
            '##FORMAT=<ID=ER2,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 2">\n',
            '##FORMAT=<ID=EA2,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 2">\n'
        )
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    ##
    ## add output
    ##
    cat(header, header2, "\n", sep="", file = output_vcf_header)
    ##
    ## build first part of file
    ##
    options(scipen=6)
    INFO <- paste0(
        "EAF=", round(estimatedAlleleFrequency, 5), ";",
        "INFO_SCORE=", round(info, 5), ";",
        "HWE=", hwe
    )
    if (sum(reference_panel_SNPs) > 0) {
        INFO <- paste0(INFO, ";REF_PANEL=", as.integer(reference_panel_SNPs))
    }
    write.table(
        matrix(paste(pos[,1], pos[,2], ".", pos[,3], pos[,4], ".", "PASS", INFO, "GT:GP:DS", sep="\t"), ncol=1),
        file = gzfile(output_vcf_left),
        row.names = FALSE,
        col.names = FALSE,
        sep = "",
        quote = FALSE
    )
    ##
    ## get number of cores and blocks
    ##
    files_to_paste <- get_files_to_paste(
        N, nCores, outputBlockSize, outputdir, vcf.piece_unique, regionName
    )
    unlisted_files_to_paste <- paste(unlist(files_to_paste), collapse = " ")    
    ##
    ## this command is crazy - write to tempfile then run
    ## merge together
    ##
    command <- paste0(
        'bash -c "',
        'paste -d ',shQuote("\t"),' ',
        '<(gunzip -c ', shQuote(output_vcf_left), ' | cat ) ',
        unlisted_files_to_paste, ' | ',
        ' cat ', shQuote(output_vcf_header) ,' - | bgzip -c > ',
        shQuote(output_vcf)
       ,'"'
    )
    ##cat(command, file = temp_runfile)
    ##system(paste0("bash ", temp_runfile))
    ## temp_runfile <- tempfile()
    system(command)
    unlink(output_vcf_header)
    unlink(output_vcf_left)
    ##
    ## remove files
    ##
    files_to_remove <- get_files_to_paste(
        N, nCores, outputBlockSize, outputdir, vcf.piece_unique, regionName, as_command = FALSE
    )
    for(file in unlist(files_to_remove)) {
        unlink(file)
    }
    return(NULL)
}



get_files_to_paste <- function(N, nCores, outputBlockSize, outputdir, vcf.piece_unique, regionName, as_command = TRUE) {
    x3 <- getSampleRange(N, nCores)
    files_to_paste <- lapply(1:length(x3), function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        files_x <- NULL
        if (nBlocks > 1) {
            for(iBlock in 2:nBlocks) {
                file_x <- file.path(
                    outputdir,
                    paste0(
                        "vcf.piece",
                        vcf.piece_unique,
                        i_core, "." , iBlock, ".",
                        regionName, ".txt.gz"
                    )
                )
                if (as_command) {
                    files_x <- c(files_x, paste0('<(gunzip -c ', shQuote(file_x), ' | cat ) '))
                } else {
                    files_x <- c(files_x, file_x)
                }
            }
            return(files_x)
        }
    })
    return(files_to_paste)
}
