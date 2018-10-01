get_sampleProbs_from_dir_for_sample <- function(
  dir,
  regionName,
  iSample,
  bundling_info,
  bundledSampleProbs = NULL
) {
    if (length(bundling_info) == 0) {
        load(file_sampleProbs(dir, iSample, regionName))
    } else {
        x <- bundling_info$matrix[iSample, ]
        iC <- x["iCore"]
        iB <- x["iBundle"]
        iP <- x["iPosition"]
        y <- bundling_info$list[[iC]][[iB]]
        s <- y[1]
        e <- y[2]
        file <- file_bundledSampleProbs(dir, s, e, regionName)
        if (length(bundledSampleProbs) == 0) {
            load(file)
        } else if (bundledSampleProbs$name != paste0(s, "_", e)) {
            load(file)
        }
        sampleProbs <- bundledSampleProbs[[iP]]
        pRgivenH1 <- sampleProbs$pRgivenH1
        pRgivenH2 <- sampleProbs$pRgivenH2
        srp <- sampleProbs$srp
    }
    return(
        list(
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,
            srp = srp,
            bundledSampleProbs = bundledSampleProbs
        )
    )
}

get_sampleReads_from_dir_for_sample <- function(
    dir,
    regionName,
    iSample,
    bundling_info,
    bundledSampleReads = NULL,
    what = "sampleReads",
    allSampleReads = NULL
) {
    if (is.null(allSampleReads) == FALSE) {
        return(
            list(
                sampleReads = allSampleReads[[iSample]],
                bundledSampleReads = NULL
            )
        )
    }
    if (what == "referenceSampleReads") {
        file_loader <- file_referenceSampleReads
        file_bundled_loader <- file_bundledReferenceSampleReads
    } else if (what == "sampleReads") {
        file_loader <- file_sampleReads
        file_bundled_loader <- file_bundledSampleReads
    }
    if (length(bundling_info) == 0) {
        load(file_loader(dir, iSample, regionName))
    } else {
        x <- bundling_info$matrix[iSample, ]
        iC <- x["iCore"]
        iB <- x["iBundle"]
        iP <- x["iPosition"]
        y <- bundling_info$list[[iC]][[iB]]
        s <- y[1]
        e <- y[2]
        file <- file_bundled_loader(dir, s, e, regionName)
        if (length(bundledSampleReads) == 0) {
            load(file)
        } else {
            if (bundledSampleReads$name != paste0(s, "_", e))
                load(file)
        }
        sampleReads <- bundledSampleReads[[iP]]
    }
    return(
        list(
            sampleReads = sampleReads,
            bundledSampleReads = bundledSampleReads
        )
    )
}


get_alphaBetaBlocks_from_dir_for_sample <- function(
    dir,
    regionName,
    iSample,
    bundling_info,
    bundledAlphaBetaBlocks = NULL,
    what = "alphaBetaBlocks",
    allAlphaBetaBlocks = NULL
) {
    if (is.null(allAlphaBetaBlocks) == FALSE) {
        return(
            list(
                alphaBetaBlocks = allAlphaBetaBlocks[[iSample]],
                bundledAlphaBetaBlocks = NULL
            )
        )
    }
    file_loader <- file_alphaBetaBlocks
    file_bundled_loader <- file_bundledAlphaBetaBlocks
    if (length(bundling_info) == 0) {
        load(file_loader(dir, iSample, regionName))
    } else {
        x <- bundling_info$matrix[iSample, ]
        iC <- x["iCore"]
        iB <- x["iBundle"]
        iP <- x["iPosition"]
        y <- bundling_info$list[[iC]][[iB]]
        s <- y[1]
        e <- y[2]
        file <- file_bundled_loader(dir, s, e, regionName)
        if (length(bundledAlphaBetaBlocks) == 0) {
            load(file)
        } else {
            if (bundledAlphaBetaBlocks$name != paste0(s, "_", e)) {
                load(file)
            }
        }
        alphaBetaBlocks <- bundledAlphaBetaBlocks[[iP]]
    }
    return(
        list(
            alphaBetaBlocks = alphaBetaBlocks,
            bundledAlphaBetaBlocks = bundledAlphaBetaBlocks
        )
    )
}


bundle_inputs_after_generation <- function(
    bundling_info,
    iBam,
    dir,
    regionName,
    what = "sampleReads",
    remove_files = TRUE
) {
    if (what == "sampleReads") {
        file_samples <- file_sampleReads
    } else if (what == "sampleProbs") {
        file_samples <- file_sampleProbs
    } else if (what == "referenceSampleReads") {
        file_samples <- file_referenceSampleReads
    } else if (what == "alphaBetaBlocks") {
        file_samples <- file_alphaBetaBlocks
    } else {
        stop("bad bundle choice")
    }
    iC <- bundling_info$matrix[iBam, "iCore"]
    iB <- bundling_info$matrix[iBam, "iBundle"]
    y <- bundling_info$list[[iC]][[iB]]
    s <- y[1]
    e <- y[2]
    ##print_message(paste0("Bundle what=", what, ", iBam=", iBam, ", iCore=", iC, ", iBundle=", iB, ", samples=", s, "-", e))
    ## now, get these samples and write to disk
    samples_in_core_bundle <- s:e
    bundledSampleObject <- lapply(
        samples_in_core_bundle,
        function(iBam) {
            if (what == "sampleReads") {
                file <- file_sampleReads(dir, iBam, regionName)
                load(file = file)
                return(sampleReads)
            } else if (what == "sampleProbs") {
                file <- file_sampleProbs(dir, iBam, regionName)
                load(file = file)
                return(list(pRgivenH1 = pRgivenH1, pRgivenH2 = pRgivenH2, srp = srp))
            } else if (what == "referenceSampleReads") {
                file <- file_referenceSampleReads(dir, iBam, regionName)
                load(file = file)
                return(sampleReads)
            } else if (what == "alphaBetaBlocks") {
                file <- file_alphaBetaBlocks(dir, iBam, regionName)
                load(file = file)
                    return(alphaBetaBlocks)
            }
        }
    )
    bundledSampleObject$name <- paste0(s, "_", e)
    if (what == "sampleReads") {
        bundledSampleReads <- bundledSampleObject
        save(
            bundledSampleReads,
            file = file_bundledSampleReads(dir, s, e, regionName),
            compress = FALSE            
        )
    } else if (what == "sampleProbs") {
        bundledSampleProbs <- bundledSampleObject
        save(
            bundledSampleProbs,
            file = file_bundledSampleProbs(dir, s, e, regionName),
            compress = FALSE            
        )
    } else if (what == "referenceSampleReads") {
        bundledSampleReads <- bundledSampleObject
        save(
            bundledSampleReads,
            file = file_bundledReferenceSampleReads(dir, s, e, regionName),
            compress = FALSE
        )
    } else if (what == "alphaBetaBlocks") {
        bundledAlphaBetaBlocks <- bundledSampleObject
        save(
            bundledAlphaBetaBlocks,
            file = file_bundledAlphaBetaBlocks(dir, s, e, regionName),
            compress = TRUE
        )
    }
    if (remove_files == TRUE) {
        for(iBam in samples_in_core_bundle) {
            unlink(file_samples(dir, iBam, regionName))
        }
    }
    return(NULL)
}



get_rebundled_files <- function(inputdir, regionName, outputdir) {
    command <- paste0(
        'cd ', shQuote(inputdir),
        ' && find . -name "',
        'bundledSamples.*-*.', regionName, '.RData', '"'
    )
    files <- file.path(inputdir, system(command, intern = TRUE))
    ## files <- system(paste0("ls ", file_bundledSampleReads(inputdir, "*", "*", regionName)), intern = TRUE)
    ## get start and end of bundledSampleReads
    a <- unlist(strsplit(files, paste0(".", regionName, ".RData")))
    x <- strsplit(a, paste0( "bundledSamples."))
    if (length(grep("bundledSamples", outputdir)) > 0)
        stop ("Please choose an outputdir name that does not include the term bundledSamples")
    ranges <-  sapply(x, function(x) x[2])
    ranges <- t(sapply(strsplit(ranges, "-"), function(x) as.integer(x)))
    files <- files[order(ranges[, 1])]
    ranges <- ranges[order(ranges[, 1]), ]
    ## also check - range is exact and fully spans
    return(
        list(
            files = files,
            ranges = ranges
        )
    )
}

## if the files already exist, good to go
## otherwise, loop through files using 1 core
## when doesn't exist, load next bunch
## every once in a while, write to disk
##
rebundle_input <- function(
    inputdir,
    regionName,
    bundling_info,
    N,
    tempdir,
    nCores,
    outputdir
) {

    print_message("Rebundle inputs")
    out <- get_rebundled_files(inputdir, regionName, outputdir)
    files <- out$files
    ranges <- out$ranges
    
    if (sum(unlist(apply(ranges, 1, function(x) x[1]:x[2])) != 1:N) > 0) {
        stop ("You are rebundling old input files, however, the originally bundled files do not span the number of files you have as inferred from sampleNames")
    }
    
    if (ranges[nrow(ranges), 2] != N) {
        stop ("You are rebundling old input files, however, the number of samples as inferred from the input bundles is not the same as sampleNames")
    }
    
    files_do_not_exist <- unlist(lapply(bundling_info$list, function(m) {
        sapply(m, function(a) {
            s <- a[1]
            e <- a[2]
            file.exists(file_bundledSampleReads(inputdir, s, e, regionName)) == FALSE
        })
    }))
    if (sum(files_do_not_exist) ==0) {
        print_message("The previously bundled files are the right size. No need to rebundle")
        print_message("Done rebundling inputs")
        return(NULL)
    }
    
    newBundle <- NULL
    bundledSampleReads <- NULL
    
    ## multi-core re-bundling
    sampleRanges <- getSampleRange(N, nCores)

    out2 <- mclapply(
        sampleRanges,
        mc.cores=nCores,
        FUN = rebundle_input_subfunction,
        ranges = ranges,
        files = files,
        tempdir = tempdir,
        regionName = regionName,
        bundling_info = bundling_info,
        inputdir = inputdir    
    )
    check_mclapply_OK(out2)
    
    print_message("Done rebundling inputs")
    return(NULL)
}

rebundle_input_subfunction <- function(
    sampleRange,
    ranges,
    files,
    tempdir,
    regionName,
    bundling_info,
    inputdir
) {

    ## start with the file relevant to the first sample
    iSample <- sampleRange[1]
    i_file <- which(iSample >= ranges[, 1] & iSample <= ranges[, 2])
    load(files[i_file])    
    cor <- ranges[i_file, ]
    cor <- cor[1]:cor[2]
    
    for(iSample in sampleRange[1]:sampleRange[2]) {
        m <- match(iSample, cor)
        if (is.na(m) == FALSE) {
            sampleReads <- bundledSampleReads[[m]]
        } else {
            i_file <- which(iSample >= ranges[, 1] & iSample <= ranges[, 2])
            load(files[i_file])
            cor <- ranges[i_file, ]
            cor <- cor[1]:cor[2]
            m <- match(iSample, cor)
            if (is.na(m)) {
                stop("There is faulty logic in rebundle_input is wrong. Please submit bug report")
            }
            sampleReads <- bundledSampleReads[[m]]
        }
        save(
            sampleReads,
            file = file_sampleReads(inputdir, iSample, regionName),
            compress = FALSE
        )
        last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
        if (last_in_bundle == 1) {
            bundle_inputs_after_generation(
                bundling_info = bundling_info,
                iBam = iSample,
                dir = inputdir,
                regionName = regionName
            )
        }
    }
    return(NULL)
}




load_all_sampleReads_into_memory <- function(
    N,
    nCores,
    tempdir,
    regionName,
    bundling_info
) {
    ##
    print_message("Begin loading all sample reads into memory")
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = function(sampleRange) {
        allSampleReads <- as.list(1:N)        
        bundledSampleReads <- NULL
        for(iSample in sampleRange[1]:sampleRange[2]) {
            out <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info,
                bundledSampleReads = bundledSampleReads
            )
            sampleReads <- out$sampleReads
            allSampleReads[[iSample]] <- sampleReads
            bundledSampleReads <- out$bundledSampleReads
        }
        return(allSampleReads)
    })
    check_mclapply_OK(out)
    ##
    allSampleReads <- as.list(1:N)
    for(i_core in 1:length(out)) {
        sampleRange <- sampleRanges[[i_core]]
        x <- out[[i_core]]        
        for(iSample in sampleRange[1]:sampleRange[2]) {        
            allSampleReads[[iSample]] <- x[[iSample]]
        }
    }
    print_message("Done loading all sample reads into memory")    
    return(allSampleReads)

}



split_reads_completely <- function(
    N,
    nCores,
    tempdir,
    regionName,
    bundling_info,
    allSampleReads
) {

    print_message("Split reads")    
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = split_read_subfunc,
        tempdir = tempdir,
        regionName = regionName,
        bundling_info = bundling_info,
        allSampleReads = allSampleReads,
        N = N
    )
    check_mclapply_OK(out)

    if (is.null(allSampleReads) == FALSE) {
        allSampleReads <- as.list(1:N)
        for(i_core in 1:length(out)) {
            sampleRange <- sampleRanges[[i_core]]
            x <- out[[i_core]]        
            for(iSample in sampleRange[1]:sampleRange[2]) {        
                allSampleReads[[iSample]] <- x[[iSample]]
            }
        }
    }

    print_message("Done splitting reads")

    return(allSampleReads)
}



split_read_subfunc <- function(
    sampleRange,
    tempdir,
    regionName,
    bundling_info,
    allSampleReads,
    N
) {

    bundledSampleReads <- NULL

    who_to_run <- sampleRange[1]:sampleRange[2]
    allSampleReadsOutput <- as.list(1:N)
    
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

        x3 <- unlist(sapply(sampleReads,function(x) x[[3]]))
        x4 <- unlist(sapply(sampleReads,function(x) x[[4]]))
        sr <- lapply(1:length(x3),function(i) {
            list(0, x4[i], x3[i], x4[i])
        })
        sampleReads <- sr[order(x4)]

        if (is.null(allSampleReads) == FALSE) {
            allSampleReadsOutput[[iSample]] <- sampleReads
        } else {
            save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName))
            if (length(bundling_info) > 0) {
                last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
                if (last_in_bundle == 1) {
                    bundle_inputs_after_generation(
                        bundling_info = bundling_info,
                        iBam = iSample,
                        dir = tempdir,
                        regionName = regionName,
                        what = "sampleReads"
                    )
                }
            }
        }

    }
    return(allSampleReadsOutput)
}




get_sampleRead_from_SNP_i_to_SNP_j <- function(
    sampleRead,
    i,
    j,
    L,
    grid
) {
    bq <- sampleRead[[3]][i:j]
    u <- sampleRead[[4]][i:j]
    ## central SNP in read
    cr <- grid[u[getCentral(L[u + 1])] + 1]
    return(
        list(
            j - i,
            cr,
            matrix(bq, ncol = 1),
            matrix(u, ncol = 1)
        )
    )
}

