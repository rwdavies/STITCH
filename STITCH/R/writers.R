## the function called by STITCH
make_and_write_output_file <- function(
    output_filename,
    outputdir,
    regionName,
    output_format,
    blocks_for_output,
    alphaBetaBlockList,
    reference_panel_SNPs,
    priorCurrent,    
    sigmaCurrent,
    alphaMatCurrent_t,
    eHapsCurrent_t,
    N,
    method,
    sampleNames,
    nSNPs,
    nGrids,
    nCores,
    B_bit_prob,
    bundling_info,
    tempdir,
    grid,
    alleleCount,
    pos,
    K,
    highCovInLow,
    start_and_end_minus_buffer,
    allSampleReads
) {

    print_message("Begin making and writing output file")
    to_use_output_filename <- get_output_filename(
        output_filename = output_filename,
        outputdir = outputdir,
        regionName = regionName,
        output_format = output_format
    )
    annot_header <- get_per_snp_annot(
        output_format = output_format,
        reference_panel_SNPs = reference_panel_SNPs,
        return_annotation_only = TRUE
    )$annot_header

    read_starts_and_ends <- determine_reads_in_output_blocks(
        N = N,
        blocks_for_output = blocks_for_output,
        nCores = nCores,
        tempdir = tempdir,
        regionName = regionName,
        bundling_info = bundling_info,
        nGrids = nGrids,
        allSampleReads = allSampleReads
    )     

    print_message("Initialize output file")
    if (output_format == "bgvcf") {
        ## put into temp file, then later bgzip (future: can I stream this?)
        output_unbgzipped <- paste0(to_use_output_filename, ".building.vcf")
        unlink(output_unbgzipped)
        make_and_write_vcf_header(
            output_vcf_header = output_unbgzipped,
            annot_header = annot_header,
            method = method,
            sampleNames = sampleNames
        )
    } else if (output_format == "bgen") {
        ## bgen header here
        out <- rrbgen_write(
            to_use_output_filename,
            sample_names = sampleNames,
            B_bit_prob = B_bit_prob,
            close_bgen_file = FALSE,
            header_M = start_and_end_minus_buffer[2] - start_and_end_minus_buffer[1] + 1
        )
        bgen_file_connection <- out$bgen_file_connection
        previous_offset <- out$final_binary_length
    }
    print_message("Done initializing output file")    

    sampleRanges <- getSampleRange(N, nCores)    
    ## write blocks
    transMatRate_t <- get_transMatRate(
        method = method,
        sigmaCurrent = sigmaCurrent
    )
    info <- array(NA, nSNPs)
    hwe <- array(NA, nSNPs)
    hweCount_total <- array(NA, c(nSNPs, 3))    
    estimatedAlleleFrequency <- array(NA, nSNPs)
    if (length(highCovInLow) > 0) {
        gen_imp <- array(NA, c(nSNPs, length(highCovInLow)))
    } else {
        gen_imp <- NULL
    }

    print_message("Loop over and write output file")
    print_i_output_block <- round(seq(1, nrow(blocks_for_output), length.out = 10))    
    for(i_output_block in 1:nrow(blocks_for_output)) {

        ## print out no more than 10 messages
        if (i_output_block %in% print_i_output_block) {
            print_message(paste0(
                "Making output piece ",
                i_output_block, " / ", nrow(blocks_for_output)
            ))
        }
        
        first_snp_in_region <- blocks_for_output[i_output_block, "snp_start_1_based"]
        last_snp_in_region <- blocks_for_output[i_output_block, "snp_end_1_based"]        
        snps_in_output_block <- first_snp_in_region:last_snp_in_region
        
        first_grid_in_region <- blocks_for_output[i_output_block, "grid_start_0_based"]
        last_grid_in_region <- blocks_for_output[i_output_block, "grid_end_0_based"]                
        grids_in_output_block <- first_grid_in_region:last_grid_in_region
        
        nSNPsInOutputBlock <- length(snps_in_output_block)

        ## shrink params
        if (first_grid_in_region < last_grid_in_region) {
            grids_to_use <- first_grid_in_region:(last_grid_in_region - 1)
            ## what? is this right?
            alphaMatCurrentLocal_t <- alphaMatCurrent_t[, 1 + grids_to_use, drop = FALSE]
            transMatRateLocal_t <- transMatRate_t[, 1 + grids_to_use, drop = FALSE]
        }
        
        ##
        out <- mclapply(
            sampleRanges,
            mc.cores = nCores,
            B_bit_prob = B_bit_prob,
            tempdir = tempdir,
            regionName = regionName,
            bundling_info = bundling_info,
            K = K,
            alphaBetaBlockList = alphaBetaBlockList,
            alphaMatCurrentLocal_t = alphaMatCurrentLocal_t,
            eHapsCurrent_t = eHapsCurrent_t,
            transMatRateLocal_t = transMatRateLocal_t,
            first_grid_in_region = first_grid_in_region,
            last_grid_in_region = last_grid_in_region,
            i_output_block = i_output_block,
            read_starts_and_ends = read_starts_and_ends,
            first_snp_in_region = first_snp_in_region,
            last_snp_in_region = last_snp_in_region,
            nSNPsInOutputBlock = nSNPsInOutputBlock,
            output_format = output_format,
            method = method,
            grid = grid,
            highCovInLow = highCovInLow,
            allSampleReads = allSampleReads,
            FUN = per_core_get_results
        )
        
        check_mclapply_OK(out)
        ## rebuild / re-assemble
        hweCount <- array(0, c(nSNPsInOutputBlock, 3))
        infoCount <- array(0, c(nSNPsInOutputBlock, 2))
        afCount <- array(0, nSNPsInOutputBlock)
        for(i in 1:length(out)) {
            infoCount <- infoCount + out[[i]]$infoCount
            hweCount <- hweCount + out[[i]]$hweCount
            afCount <- afCount + out[[i]]$afCount
        }
        if (output_format == "bgen") {
            list_of_gp_raw_t <- lapply(1:length(out), function(i_core) {
                gp_raw_t <- out[[i_core]]$gp_raw_t
                return(gp_raw_t)
            })
        } else if (output_format == "bgvcf") {
            ## make the WHOLE THING here
            vcf_matrix_to_out <- data.frame(matrix(
                data = NA,
                nrow = nSNPsInOutputBlock,
                ncol = N + 9
            ))
            for(i in 1:length(out)) {
                sampleRange <- sampleRanges[[i]]
                vcf_matrix_to_out[, 9 + sampleRange[1]:sampleRange[2]] <- out[[i]]$vcf_matrix_to_out
            }
        }
        ##
        rm(out)
        if (N > 10000) {
            gc(reset = TRUE); gc(reset = TRUE)
        }

        if (length(highCovInLow) > 0) {
            for(j in 1:length(highCovInLow)) {
                ## get gp_t
                load(file = file_dosages(tempdir, highCovInLow[j], regionName, "piece.gp_t"))
                gen_imp[snps_in_output_block, j] <- 0.5 * gp_t[2, ] + gp_t[3, ]
            }
        }
        
        ## now make HWE etc
        thetaHat <- infoCount[, 1] / 2 / N
        denom <- 2 * N * thetaHat * (1-thetaHat)
        info[snps_in_output_block] <- 1 - infoCount[, 2] / denom
        info[snps_in_output_block][(thetaHat == 0) | (thetaHat == 1)] <- 1
        estimatedAlleleFrequency[snps_in_output_block] <- afCount / N
        hwe[snps_in_output_block] <- generate_hwe_on_counts(hweCount, nSNPsInOutputBlock, nCores)
        hweCount_total[snps_in_output_block, ] <- hweCount

        ## now, can output!
        out <- get_per_snp_annot(
            output_format = output_format,
            reference_panel_SNPs = reference_panel_SNPs[snps_in_output_block],
            estimatedAlleleFrequency = estimatedAlleleFrequency[snps_in_output_block],
            info = info[snps_in_output_block],
            hwe = hwe[snps_in_output_block],
            alleleCount = alleleCount[snps_in_output_block, , drop = FALSE]
        )
        
        INFO <- out$INFO

        ## now write block for
        if (output_format == "bgvcf") {
            
            vcf_matrix_to_out[, 1] <- pos[snps_in_output_block, 1]
            vcf_matrix_to_out[, 2] <- pos[snps_in_output_block, 2]
            vcf_matrix_to_out[, 3] <- "."
            vcf_matrix_to_out[, 4] <- pos[snps_in_output_block, 3]
            vcf_matrix_to_out[, 5] <- pos[snps_in_output_block, 4]
            vcf_matrix_to_out[, 6] <- "."
            vcf_matrix_to_out[, 7] <- "PASS"
            vcf_matrix_to_out[, 8] <- INFO
            vcf_matrix_to_out[, 9] <- "GT:GP:DS"
            data.table::fwrite(
                vcf_matrix_to_out,
                file = output_unbgzipped,
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                append = TRUE,
                nThread = nCores
            )
            
        } else if (output_format == "bgen") {

            var_info <- make_var_info(pos, c(first_snp_in_region, last_snp_in_region))
            out <- rrbgen_write(
                bgen_file_connection = bgen_file_connection,
                previous_offset = previous_offset,
                add_to_bgen_connection = TRUE,
                close_bgen_file = FALSE,
                list_of_gp_raw_t = list_of_gp_raw_t,
                sample_names = sampleNames,
                var_info = var_info,
                B_bit_prob = B_bit_prob,
                nCores = nCores
            )
            bgen_file_connection <- out$bgen_file_connection
            previous_offset <- out$final_binary_length

        }

    }
    print_message("Done looping over and writing output file")
    
    ## 
    if (output_format == "bgvcf") {
        print_message("bgzip output file and move to final location")
        system(paste0("bgzip --threads ", nCores, " -f ", shQuote(output_unbgzipped)))
        system(paste0("mv ", shQuote(paste0(output_unbgzipped, ".gz")), " ", shQuote(to_use_output_filename)))
    } else {
        close(bgen_file_connection)
        print_message("Write out variant statistics to accompany bgen file")
        var_info <- make_var_info(pos, start_and_end_minus_buffer)
        make_and_write_bgen_per_snp_annot_file(
            to_use_output_filename = to_use_output_filename,
            estimatedAlleleFrequency = estimatedAlleleFrequency,
            info = info,
            hwe = hwe,
            reference_panel_SNPs = reference_panel_SNPs,
            start_and_end_minus_buffer = start_and_end_minus_buffer,
            nSNPs = nSNPs,
            alleleCount = alleleCount,
            var_info = var_info
        )
        print_message("Done writing out variant statistics to accompany bgen file")
    }
    
    print_message("Done making and writing output file")

    return(
        list(
            gen_imp = gen_imp,
            estimatedAlleleFrequency = estimatedAlleleFrequency,
            info = info,
            hwe = hwe,
            hweCount = hweCount_total
        )
    )
    
}

per_core_get_results <- function(
    sampleRange,
    B_bit_prob,
    tempdir,
    regionName,
    bundling_info,
    K,
    alphaBetaBlockList,
    alphaMatCurrentLocal_t,
    eHapsCurrent_t,
    transMatRateLocal_t,
    first_grid_in_region,
    last_grid_in_region,    
    i_output_block,
    read_starts_and_ends,
    first_snp_in_region,
    last_snp_in_region,
    nSNPsInOutputBlock,
    output_format,
    method,
    grid,
    highCovInLow,
    allSampleReads
) {

    bundledSampleReads <- NULL
    bundledSampleProbs <- NULL
    
    ## load sample
    ## load pRgivenH1
    who_to_run <- sampleRange[1]:sampleRange[2]
    N_core <- sampleRange[2] - sampleRange[1] + 1 ## number in this core
    hweCount <- array(0, c(nSNPsInOutputBlock, 3))
    infoCount <- array(0, c(nSNPsInOutputBlock, 2))
    afCount <- array(0, nSNPsInOutputBlock)
    
    if (output_format == "bgvcf") {
        vcf_matrix_to_out <- data.frame(matrix(
            data = NA,
            nrow = nSNPsInOutputBlock,
            ncol = N_core
        ))
        gp_raw_t <- NULL
    } else if (output_format == "bgen") {
        gp_raw_t <- matrix(
            data = raw(0),
            nrow = N_core * 2 * (B_bit_prob / 8),
            ncol = nSNPsInOutputBlock
        )
        vcf_matrix_to_out <- NULL
    }
    pRgivenH1 <- NULL
    pRgivenH2 <- NULL

    for (iiSample in 1:(length(who_to_run))) {
        
        iSample <- who_to_run[iiSample]
        ## get these from list?
        out <- get_sampleReads_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = iSample,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads,
            allSampleReads = allSampleReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads
        
        if (method == "pseudoHaploid") {
            out <- get_sampleProbs_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info,
                bundledSampleProbs = bundledSampleProbs
            )
            pRgivenH1 <- out$pRgivenH1
            pRgivenH2 <- out$pRgivenH2
            srp <- out$srp
            bundledSampleProbs <- out$bundledSampleProbs
        }

        ## run forward backwards
        s <- read_starts_and_ends[iSample, i_output_block, 1]
        e <- read_starts_and_ends[iSample, i_output_block, 2]
        if (is.na(s) | is.na(e)) {
            which_reads <- NULL
        } else {
            which_reads <- s:e
        }

        if (first_grid_in_region == last_grid_in_region) {
            ## only 1 grid being imputed (rare!)
            ## also, imputation more or less meaningless
            ## but anyway, make valid output
            abSmall <- alphaBetaBlockList[[iSample]]
            fbsoL <- lapply(1:length(abSmall), function(i) {
                a <- abSmall[[i]]$alphaHatBlocks_t
                b <- abSmall[[i]]$betaHatBlocks_t
                gamma_t <- a[, i_output_block, drop = FALSE] * b[, i_output_block, drop = FALSE]
                gamma_t <- gamma_t / sum(gamma_t) ## don't have c here
                return(list(gamma_t = gamma_t))
            })
        } else {    

            fbsoL <- run_forward_backwards(
                sampleReads = sampleReads[which_reads],
                pRgivenH1 = pRgivenH1[which_reads],
                pRgivenH2 = pRgivenH2[which_reads],
                method = method,
                K = K,
                priorCurrent = array(-1, K), ## irrelevant here
                alphaMatCurrent_t = alphaMatCurrentLocal_t,
                eHapsCurrent_t = eHapsCurrent_t,
                transMatRate_t = transMatRateLocal_t,
                alphaBetaBlock = alphaBetaBlockList[[iSample]],
                i_snp_block_for_alpha_beta = i_output_block,
                run_fb_grid_offset = first_grid_in_region,
                run_fb_subset = TRUE
            )$fbsoL
            
        }

        gp_t <- calculate_gp_t_from_fbsoL(
            eHapsCurrent_t = eHapsCurrent_t,            
            grid = grid,
            method = method,
            fbsoL = fbsoL,
            snp_start_1_based = first_snp_in_region,
            snp_end_1_based = last_snp_in_region,
            grid_offset_0_based = first_grid_in_region
        )

        if (iSample %in% highCovInLow) {
            save(gp_t, file = file_dosages(tempdir, iSample, regionName, "piece.gp_t"))
        }
        
        ## 
        eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
        fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##

        infoCount[, 1] <- infoCount[, 1, drop = FALSE] + eij
        infoCount[, 2] <- infoCount[, 2, drop = FALSE] + (fij - eij**2)
        ## this returns un-transposed results
        max_gen <- get_max_gen_rapid(gp_t)
        ## hweCount is NOT transposed!
        hweCount[max_gen] <- hweCount[max_gen] + 1 ## hmmmmm not ideal
        afCount <- afCount + (eij) / 2
        ##
        if (output_format == "bgvcf") {
            vcf_matrix_to_out[, iiSample] <- rcpp_make_column_of_vcf(gp_t, 0, matrix())
        } else if (output_format == "bgen") {
            rrbgen::rcpp_place_gp_t_into_output(
                gp_t,
                gp_raw_t, ## storage matrix
                iiSample, ## relative position
                nSNPs = ncol(gp_t), ## ncol(gp_t) = nSNPsInRegion here
                B_bit_prob
            )
            ##
        }
        
    }
    
    return(
        list(
            gp_raw_t = gp_raw_t,
            vcf_matrix_to_out = vcf_matrix_to_out,
            infoCount = infoCount,
            hweCount = hweCount,
            afCount = afCount
        )
    )
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
    output_format,
    allSampleReads
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
                bundledSampleReads = bundledSampleReads,
                allSampleReads = allSampleReads
            ) ## not yet loaded into RAM
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
        ' cat ', shQuote(output_vcf_header) ,' - | bgzip --threads ', nCores, '-c > ',
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



## write out a single samples worth of VCF entry
make_column_of_vcf_from_pl_and_rd <- function(
    pl,
    rd
) {
    ## write out genotype, genotype likelihood, and dosage
    ## gt, ad, pl
    ##GT:AD:PL
    ##chr19   3126522 .       G       A       7661.22 .       AC=464;AF=0.935;AN=496;BaseQRankSum=4.232;DP=333;Dels=0.00;FS=2.773;GC=39.90;HRun=0;HaplotypeScore=0.0738;InbreedingCoeff=-0.0083;MLEAC=465;MLEAF=0.938;MQ=42.81;MQ0=23;MQRankSum=-0.576;QD=28.69;ReadPosRankSum=-0.673      GT:AD:DP:GQ:PL       ./.     ./.     ./.     1/1:0,1:1:3:32,3,0      ./.     ./.
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
    return(cbind(1:ncol(x), z))
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
    gt <- c("0/0","0/1","1/1")[z[, 2]]
    z2 <- cbind(z[, 2], z[, 1]) ## flip - argh
    gt[gp_t[z2] < 0.9] <- "./."
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
    if (output_format == "bgvcf") {
        extension <- ".vcf.gz"
    } else {
        extension <- ".bgen"
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
            return(path.expand(output_filename))
        }
    }
}


make_var_info <- function(pos, start_and_end_minus_buffer) {
    pos_local <- pos[start_and_end_minus_buffer[1]:start_and_end_minus_buffer[2], , drop = FALSE]
    var_info <- matrix("", nrow = nrow(pos_local), ncol = 6)
    colnames(var_info) <- c("chr", "varid", "rsid", "position", "ref", "alt")
    var_info[, "chr"] <- as.character(pos_local[, "CHR"])
    var_info[, "varid"] <- paste0(pos_local[, "CHR"], ":", pos_local[, "POS"], ":", pos_local[, "REF"], ":", pos_local[, "ALT"])
    var_info[, "rsid"] <- "."
    var_info[, "position"] <- pos_local[, "POS"]
    var_info[, "ref"] <- as.character(pos_local[, "REF"])
    var_info[, "alt"] <- as.character(pos_local[, "ALT"])
    return(var_info)
}

make_and_write_vcf_header <- function(output_vcf_header, annot_header, method, sampleNames) {
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">\n',
        annot_header,
        '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">\n'
    )
    ## disable for now    
    if (method == "pseudoHaploid" && 1 == 0) {
        header <- paste0(
            header,
            '##FORMAT=<ID=ER1,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 1">\n',
            '##FORMAT=<ID=EA1,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 1">\n',
            '##FORMAT=<ID=ER2,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 2">\n',
            '##FORMAT=<ID=EA2,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 2">\n'
        )
    }
    header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
    ##
    ## add output
    ##
    cat(header, header2, "\n", sep="", file = output_vcf_header)
    return(NULL)
}




## centralize this to assure that whatever the output method, the same annotations are written to both files
get_per_snp_annot <- function(
    output_format,
    reference_panel_SNPs,                              
    estimatedAlleleFrequency = NULL,
    info = NULL,
    hwe = NULL,
    alleleCount = NULL,
    return_annotation_only = FALSE    
) {
    INFO <- NULL
    if (output_format == "bgvcf") {
        if (!return_annotation_only) {
            INFO <- paste0(
                "EAF=", round(estimatedAlleleFrequency, 5), ";",
                "INFO_SCORE=", round(info, 5), ";",
                "HWE=", formatC(hwe, format = "e", digits = 2), ";",
                "ERC=", round(alleleCount[, 1], 5), ";",
                "EAC=", round(alleleCount[, 2] - alleleCount[, 1], 5), ";",
                "PAF=", round(alleleCount[, 3], 5)
            )
        }
        annot_header <- paste0(
            '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score">\n',
            '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
            '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">\n',
            '##INFO=<ID=ERC,Number=.,Type=Float,Description="Estimated number of copies of the reference allele from the pileup">\n',
            '##INFO=<ID=EAC,Number=.,Type=Float,Description="Estimated number of copies of the alternate allele from the pileup">\n',
            '##INFO=<ID=PAF,Number=.,Type=Float,Description="Estimated allele frequency using the pileup of reference and alternate alleles">\n'
        )
        if (sum(reference_panel_SNPs) > 0) {
            if (!return_annotation_only) {
                INFO <- paste0(INFO, ";REF_PANEL=", as.integer(reference_panel_SNPs))
            }
            annot_header <- paste0(
                annot_header,
                '##INFO=<ID=REF_PANEL,Number=.,Type=Integer,Description="Whether a SNP was (1) or was not (0) found in the reference panel during imputation">\n'
            )
        }
    } else {
        if (!return_annotation_only) {
            INFO <- data.frame(
                "EAF" =  round(estimatedAlleleFrequency, 5),
                "INFO_SCORE" = round(info, 5),
                "HWE" = formatC(hwe, format = "e", digits = 2),
                "ERC" = round(alleleCount[, 1], 5),
                "EAC" = round(alleleCount[, 2] - alleleCount[, 1], 5),
                "PAF" = round(alleleCount[, 3], 5)
            )
        }
        annot_header <- c(
            '#INFO_SCORE Imputation info score, same as IMPUTE info measure I_A',
            '#EAF Estimated allele frequency after imputation',
            '#HWE Hardy-Weinberg p-value',
            '#ERC Estimated number of copies of the reference allele from the pileup',
            '#EAC Estimated number of copies of the alternate allele from the pileup',
            '#PAF Estimated allele frequency using the pileup of reference and alternate alleles'
        )
        if (sum(reference_panel_SNPs) > 0) {
            if (!return_annotation_only) {
                INFO$REF_PANEL <- as.integer(reference_panel_SNPs)
            }
            annot_header <- c(
                annot_header,
                '#REF_PANEL Whether a SNP was (1) or was not (0) found in the reference panel during imputation'
            )
        }
    }
    return(
        list(
            INFO = INFO,
            annot_header = annot_header
        )
    )
}


make_and_write_bgen_per_snp_annot_file <- function(
    to_use_output_filename,
    estimatedAlleleFrequency,
    info,
    hwe,
    reference_panel_SNPs,
    start_and_end_minus_buffer,
    nSNPs,
    alleleCount,
    var_info
) {
    ## first six columns match gen spec
    ## next columns give more
    ## this returns a data.frame for bgen
    out <- get_per_snp_annot(
        output_format = "bgen",
        reference_panel_SNPs = reference_panel_SNPs,
        estimatedAlleleFrequency = estimatedAlleleFrequency,
        info = info,
        hwe = hwe,
        alleleCount = alleleCount
    )
    INFO <- out$INFO[start_and_end_minus_buffer[1]:start_and_end_minus_buffer[2], , drop = FALSE]    
    annot_header <- out$annot_header
    out_file <- paste0(to_use_output_filename, ".per_snp_stats.txt")
    annot_header <- c(
        paste0("#STITCH ", sessionInfo()$otherPkgs$STITCH$Version, " per-SNP imputation statistics and other annotations"),
        paste0("#First 6 columns are from .gen specification, where VAR_ID=<chr>:<position>:<ref>:<alt>. The rest are as follows"),
        annot_header
    )
    cat(annot_header, file = out_file, sep = "\n")
    ## add in POS and VARID defined as above
    suppressWarnings(write.table(
        cbind(var_info, INFO),
        file = out_file,
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE,
        append = TRUE
    ))
    system(paste0("gzip -f -1 '", out_file, "'"))
    return(NULL)
}
    



## this assumes reads are sorted but that need not be the case for old STITCH
## 
determine_reads_in_output_blocks <- function(
    N,
    blocks_for_output,
    nCores,
    tempdir,
    regionName,
    bundling_info,
    nGrids,
    allSampleReads
) {
    print_message("Determine reads in output blocks")
    ##
    blocks_in_vector_form <- array(NA, nGrids) ## want 0, 1, 2, etc, depending on which output block
    for(i in 1:nrow(blocks_for_output)) {
        s2 <- blocks_for_output[i, "grid_start_0_based"]
        e2 <- blocks_for_output[i, "grid_end_0_based"]
        blocks_in_vector_form[1 + s2:e2] <- i
    }
    sampleRanges <- getSampleRange(N = N, nCores = nCores)    
    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        tempdir = tempdir,
        regionName = regionName,
        bundling_info = bundling_info,
        blocks_for_output = blocks_for_output,
        blocks_in_vector_form = blocks_in_vector_form,
        allSampleReads = allSampleReads,
        FUN = determine_reads_in_output_blocks_subfunction
    )
    check_mclapply_OK(out)
    ## rebuild
    read_starts_and_ends <- array(NA, c(N, nrow(blocks_for_output), 2))
    for(i_core in 1:length(out)) {
        who_to_run <- sampleRanges[[i_core]][1]:sampleRanges[[i_core]][2]
        outL <- out[[i_core]]
        for (iiSample in 1:(length(who_to_run))) {
            read_starts_and_ends[who_to_run[iiSample], , 1] <- outL[iiSample, , 1, drop = FALSE]
            read_starts_and_ends[who_to_run[iiSample], , 2] <- outL[iiSample, , 2, drop = FALSE]
        }
    }
    print_message("Done determining reads in output blocks")    
    return(read_starts_and_ends)
}

determine_reads_in_output_blocks_subfunction <- function(
    sampleRange,
    tempdir,
    regionName,
    bundling_info,
    blocks_in_vector_form,
    blocks_for_output,
    allSampleReads
) {
    ## 
    bundledSampleReads <- NULL
    who_to_run <- sampleRange[1]:sampleRange[2]
    N_core <- length(who_to_run)
    ##
    read_starts_and_ends <- array(NA, c(N_core, nrow(blocks_for_output), 2))
    ##
    for (iiSample in 1:(length(who_to_run))) {
        ##
        iSample <- who_to_run[iiSample]
        ## get these from list?
        out <- get_sampleReads_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = iSample,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads,
            allSampleReads = allSampleReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads
        ##
        central_snp_in_read <- sapply(sampleReads, function(x) x[[2]]) + 1
        ## if not sorted, throw error
        if (sum(diff(central_snp_in_read) < 0) > 0) {
            x <- which.max(diff(central_snp_in_read) < 0)            
            stop(paste0("In an internal format used by STITCH, the reads are not sorted by the central SNP in the read. For example, sample number ", iSample, " has read number ", x, " with central SNP ", sampleReads[[x]][[2]] + 1, " and read number ", x + 1, " with central SNP ", sampleReads[[x + 1]][[2]] + 1, ". This is probably caused by using previously generated input data from an older version of STITCH. Please regenerate input and try again. Otherwise, one can manually convert by ordering the sampleReads by the 2nd entry of each list member and saving again. If this does not seem correct please get in touch"))
        }
        what_output_block_read_is_in <- blocks_in_vector_form[central_snp_in_read]
        s <- diff(what_output_block_read_is_in) != 0 ## where they switch
        ## start and end of regions
        not_NA <- which(is.na(what_output_block_read_is_in) == FALSE)
        starts <- c(head(not_NA, 1), which(s) + 1)
        ends <- c(which(s), tail(not_NA, 1))
        whats <- what_output_block_read_is_in[starts]
        ## put into matrix?
        read_starts_and_ends[iiSample, whats, 1] <- starts
        read_starts_and_ends[iiSample, whats, 2] <- ends
    }
    return(read_starts_and_ends)
}


    
