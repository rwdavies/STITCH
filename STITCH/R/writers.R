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
    vcf_output_name,
    vcf.piece_unique
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
        vcf_matrix_to_out <- array(
            NA,
            c(nSNPs, outputBlockRange[2] - outputBlockRange[1])
        )
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
                    vcf_matrix_to_out <- array(
                        NA,
                        c(nSNPs, outputBlockRange[iBlock] - outputBlockRange[iBlock - 1])
                    )
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
    output_vcf <- get_output_vcf(vcf_output_name,  outputdir, regionName, prefix = "stitch.input")
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
    x3 <- getSampleRange(N, nCores)
    files_to_paste <- lapply(1:nCores, function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(
            paste0(
                '<(gunzip -c ', outputdir, "vcf.piece",
                vcf.piece_unique, i_core,
                ".", 2:nBlocks, ".", regionName,
                ".txt.gz", ' | cat ) '
            )
        )
    })
    files_to_paste <- paste(unlist(files_to_paste), collapse = " ")
    ##
    ## merge together
    ##
    command <- paste0(
        'bash -c "',
        'paste -d ',shQuote("\t"),' ',
        '<(gunzip -c ',output_vcf_left,' | cat ) ',
        files_to_paste, " | ",
        " cat ",output_vcf_header," - | bgzip -c > ",
        output_vcf
       ,'"'
    )
    system(command)
    system(paste0("rm ",output_vcf_header))
    system(paste0("rm ",output_vcf_left))
    ##
    ## remove
    ##
    x3 <- getSampleRange(N, nCores)
    files_to_remove <- unlist(lapply(1:nCores, function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(
            paste0(
                outputdir, "vcf.piece", vcf.piece_unique,
                i_core,
                ".", 2:nBlocks, ".", regionName, ".txt.gz"
            )
        )
    }))
    system(paste0("rm ",paste0(files_to_remove, collapse=" ")))
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
    ##m <- matrix(
    ##    sapply(list_of_vcf_columns_to_out[iSample_list], I),
    ##    ncol = length(iSample_list)
    ##)
    write.table(
        vcf_matrix_to_out,
        file = gzfile(
            paste0(
                outputdir, "vcf.piece",
                vcf.piece_unique,
                i_core, ".", iBlock, ".", regionName, ".txt.gz"
            )
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
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




get_max_gen_rapid <- function(x) {
  # assume matrix 3 columns >1 row
  z <- rep(1, nrow(x))
  y <- x[,1]
  for(i in 2:3) {
    w <- x[,i]>y
    z[w] <- i
    y[w] <- x[w,i]
  }
  return(cbind(1:nrow(x), z))
}



## write out a single samples worth of VCF entry
make_column_of_vcf <- function(
    gp,
    read_proportions
) {
    ## write out genotype, genotype likelihood, and dosage
    ##GT:GL:DS
    ##FORMAT=<ID=GT:,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
    ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Posterior probability of 0/0, 0/1, and 1/1">
    ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
    ## 1/1:0,0.054,0.946:1.946
    ## add one samples worth of info to a VCF
    z <- get_max_gen_rapid(gp)
    gt <- c("0/0","0/1","1/1")[z[,2]]
    gt[gp[z] < 0.9] <- "./."
    precision <- 3
    str <- paste0(
        gt,":",
        round(gp[,1], precision), ",",
        round(gp[,2], precision), ",",
        round(gp[,3], precision), ":",
        round(gp[,2] + 2 * gp[,3], precision)
    )
    if (is.null(read_proportions) == FALSE)
        str <- paste0(
            str, ":",
            round(read_proportions[, 1], precision), ",",
            round(read_proportions[, 2], precision), ",",
            round(read_proportions[, 3], precision), ",",
            round(read_proportions[, 4], precision)
        )
    return(str)
}



## specify here the output VCF from the process
get_output_vcf <- function(
    vcf_output_name,
    outputdir,
    regionName,
    prefix = NULL
) {
    if (is.null(vcf_output_name)) {
        if (is.null(prefix)) {
            return(file.path(
                outputdir,
                paste0("stitch.", regionName, ".vcf.gz")
            ))
        } else {
            return(file.path(
                outputdir,
                paste0(prefix, ".", regionName, ".vcf.gz")
            ))
        }
    } else {
        if (basename(vcf_output_name) == vcf_output_name) {
            return(file.path(
                outputdir,
                vcf_output_name
            ))
        } else {
            return(vcf_output_name)
        }
    }
}

# build the header
# build the left part of the VCF
# paste together the interim VCF files
# gzip the whole thing together
write_vcf_after_EM <- function(
  vcf_output_name,
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
  vcf.piece_unique = "."
) {
    ## set up file names
    print_message("Build final VCF")
    output_vcf <- get_output_vcf(vcf_output_name,  outputdir, regionName)
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
    x3 <- getSampleRange(N, nCores)
    files_to_paste <- lapply(1:length(x3), function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(
            paste0(
                '<(gunzip -c ', outputdir, "vcf.piece",
                vcf.piece_unique,
                i_core, "." , 2:nBlocks, ".",
                regionName, ".txt.gz", ' | cat ) '
            )
        )
    })
    files_to_paste <- paste(unlist(files_to_paste), collapse = " ")
    ##
    ## merge together
    ##
    command <- paste0(
        'bash -c "',
        'paste -d ',shQuote("\t"),' ',
        '<(gunzip -c ',output_vcf_left,' | cat ) ',
        files_to_paste, " | ",
        " cat ",output_vcf_header," - | bgzip -c > ",
        output_vcf
       ,'"'
    )
    system(
        command
    )
    system(paste0("rm ",output_vcf_header))
    system(paste0("rm ",output_vcf_left))
    ##
    ## remove
    ##
    x3 <- getSampleRange(N, nCores)
    files_to_remove <- unlist(lapply(1:length(x3), function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(
            paste0(
                outputdir, "vcf.piece",
                vcf.piece_unique, i_core,
                ".", 2:nBlocks, ".", regionName, ".txt.gz"
            )
        )
    }))
    system(paste0("rm ",paste0(files_to_remove, collapse=" ")))
    return(NULL)
}



