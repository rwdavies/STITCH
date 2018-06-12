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
    what = "sampleReads"
) {
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



bundle_inputs_after_generation <- function(
    bundling_info,
    iBam,
    dir,
    regionName,
    what = "sampleReads",
    remove_files = TRUE
    ) {
    if (what == "sampleReads")
        file_samples <- file_sampleReads
    if (what == "sampleProbs")
        file_samples <- file_sampleProbs
    if (what == "referenceSampleReads")
        file_samples <- file_referenceSampleReads
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
    }
    if (remove_files == TRUE) {
        for(iBam in samples_in_core_bundle) {
            unlink(file_samples(dir, iBam, regionName))
        }
    }
    return(NULL)
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
    environment
) {

    print_message("Rebundle inputs")
    out <- get_rebundled_files(inputdir, regionName)
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
    
    rebundle_input
    newBundle <- NULL
    bundledSampleReads <- NULL
    
    
    ## multi-core re-bundling
    sampleRanges <- getSampleRange(N, nCores)

    f <- function(sampleRange, ranges, files, tempdir, regionName, bundling_info) {
        ## start with the first file
        i_file <- 1
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
                if (is.na(m)) stop("There is faulty logic in rebundle_input is wrong. Please submit bug report")
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
    }


    if(environment=="server") {
        out2 <- mclapply(sampleRanges, mc.cores=nCores,FUN=f, ranges = ranges, files  = files, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info)
    }
    if(environment=="cluster") {
        cl = makeCluster(nCores, type = "FORK")
        out2 = parLapply(cl, sampleRanges, fun=f, ranges = ranges, files  = files, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info)
        stopCluster(cl)
    }
    error_check <- sapply(out2, class) == "try-error"
    if (sum(error_check) > 0) {
        print_message(out2[[which(error_check)[1]]])
        stop("There has been an error rebundling the input. Please see error message above")
    }
    
    print_message("Done rebundling inputs")
    return(NULL)
}
