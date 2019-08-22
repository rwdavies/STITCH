## can I compile with logical matrix and do operations?

test_that("can profile", {

    skip("turn off for now")
    ## so what if I want to do large scale comparison later
    ## e.g. compare numeric hap to big matrix (scaled?)
    
    if (1 == 0) {

        ## here, each
        plot_shuffle_haplotype_attempts <- TRUE
        sampleRanges <- list(c(1, 416))
        tempdir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/temp/"
        system(paste0("rsync -av ", tempdir(), "/* ", tempdir))
        ## reference_haps,
        ## non_NA_cols,
        non_NA_cols <- which(is.na(reference_haps[ , 1]) == FALSE)
        raw_reference_haps <- make_raw_reference_haps(reference_haps)
        
        save(
            raw_reference_haps,
            reference_haps, non_NA_cols,
            tempdir, regionName, reference_bundling_info,plot_shuffle_haplotype_attempts,
            eHapsCurrent_tc,
            alphaMatCurrent_tc,
            transMatRate_tc_H,
            priorCurrent_m,
            grid,
            list_of_break_results,
            reference_phred,
            nCores,
            sampleRanges,
            nbreaks,
            iteration,
            file = "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/temp_speed.RData")

    }

    load("/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/temp_speed.RData")

    plot_shuffle_haplotype_attempts <- FALSE
    suppressOutput <- 1
    print("Starting new")
    a1 <- Sys.time()
    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = new_subset_of_single_reference_iteration,
        tempdir = tempdir,
        regionName = regionName,
        reference_bundling_info = reference_bundling_info,
        eHapsCurrent_tc =eHapsCurrent_tc,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        priorCurrent_m = priorCurrent_m,
        grid = grid,
        list_of_break_results = list_of_break_results,
        reference_phred = reference_phred,
        nbreaks = nbreaks,
        sampleRanges = sampleRanges,
        reference_haps = reference_haps,        
        non_NA_cols = non_NA_cols,
        suppressOutput = suppressOutput,
        plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts,
        iteration = iteration
    )
    b1 <- Sys.time()
    print(difftime(b1, a1))
    print("Done new")

    if (plot_shuffle_haplotype_attempts) {
        load(file_fbdStore(tempdir, regionName, iteration))
        new_list_of_fbd_store <- list_of_fbd_store
        unlink(file_fbdStore(tempdir, regionName, iteration))
    }
    
    print("Starting old")
    a2 <- Sys.time()
    out2 <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = subset_of_single_reference_iteration,
        tempdir = tempdir,
        regionName = regionName,
        reference_bundling_info = reference_bundling_info,
        eHapsCurrent_tc =eHapsCurrent_tc,
        alphaMatCurrent_tc = alphaMatCurrent_tc,
        transMatRate_tc_H = transMatRate_tc_H,
        priorCurrent_m = priorCurrent_m,
        grid = grid,
        list_of_break_results = list_of_break_results,
        reference_phred = reference_phred,
        nbreaks = nbreaks,
        sampleRanges = sampleRanges,
        plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts,
        suppressOutput = suppressOutput,
        iteration = iteration
    )
    b2 <- Sys.time()
    print(difftime(b2, a2))
    print("Done old")
    print(names(out2[[1]]))

    if (plot_shuffle_haplotype_attempts) {    
        load(file_fbdStore(tempdir, regionName, iteration))
        old_list_of_fbd_store <- list_of_fbd_store
        unlink(file_fbdStore(tempdir, regionName, iteration))
    }
    
    print("--------------")
    print(paste0("New time:", difftime(b1, a1)))
    print(paste0("Old time:", difftime(b2, a2)))
    
    expect_equal(out, out2)
    if (plot_shuffle_haplotype_attempts) {    
        expect_equal(new_list_of_fbd_store, old_list_of_fbd_store)
    }

    save(out, out2, file = "~/temp.RData")
    ## new_list_of_fbd_store, old_list_of_fbd_store, 

    if (1 == 0) {

        load(file = "~/temp.RData")
        i <- 4
        expect_equal(
            out[[1]]$list_of_fromMat[[1]][i, , ],
            out2[[1]]$list_of_fromMat[[1]][i, , ]
        )
        
        ## check out as well
        s <- 1
        iSample <- 3
        new_list_of_fbd_store[[s]][[iSample]][["gammaK_t"]][1:4, 1:4]
        old_list_of_fbd_store[[s]][[iSample]][["gammaK_t"]][1:4, 1:4]

        old_list_of_fbd_store
        
        out
        ## OK!
        load("~/temp.old.RData")
        a <- fbsoL$eMatGrid_t
        load("~/temp.new.RData")
        b1 <- fbsoL$eMatGrid_t_new
        b2 <- fbsoL$eMatGrid_t_old

        ## OK - difference here! quite severe. understand...
        a[1:4, 28 + -2:2]
        b1[1:4, 28 + -2:2]
        b2[1:4, 28 + -2:2]        
        eMatGrid_t1[1:4, 28 + -2:2]
        
        eHapsCurrent_tc[1:4, 28, 2]        
        range(a[, 30])
        range(a[, 28])
        range(a[, 30])        
        
    }
    
})
    




if (1 == 0) {

    library("testthat"); library("STITCH"); library("rrbgen")
    dir <- "/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH"
    setwd(paste0(dir, "/STITCH/R"))
    a <- dir(pattern = "*.R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)
    setwd(dir)
    Sys.setenv(PATH = paste0(getwd(), ":", Sys.getenv("PATH")))

    
    load("/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH/temp_speed.RData")
    K <- dim(eHapsCurrent_tc)[1]
    nGrids <- dim(eHapsCurrent_tc)[2]
            eMatGrid_t1 <- array(1, c(K, nGrids))
            eMatGrid_t2 <- array(1, c(K, nGrids))        
            s <- 1
            iSample <- 0
            reference_phred <- 20
            maxEmissionMatrixDifference <- 1e10
    
            rcpp_ref_make_eMatGrid_t(
                eMatGrid_t = eMatGrid_t1,
                reference_haps = reference_haps,
                non_NA_cols = non_NA_cols,
                eHapsCurrent_tc = eHapsCurrent_tc,
                grid = grid,
                reference_phred = reference_phred,
                s = s,
                iSample = iSample,
                maxEmissionMatrixDifference = maxEmissionMatrixDifference,
                rescale = TRUE,
                bound = TRUE
            )
    eMatGrid_t1[1:4, 28]
    
            ## compare against old way
            sampleReads <- rcpp_make_sampleReads_from_hap(
                non_NA_cols = non_NA_cols,
                reference_phred = reference_phred,
                reference_hap = reference_haps[, iSample + 1]
            )
            ## 
            sampleReads <- snap_sampleReads_to_grid(
                sampleReads = sampleReads,
                grid = grid
            )
            ## 
            eMatRead_t <- array(1, c(K, length(sampleReads)))
            rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads = 1e10, Jmax = 100, eMatHapOri_t = array(0, c(1, 1)), pRgivenH1 = array(0, 1), pRgivenH2 = array(0, 1), prev = 0, suppressOutput = 1, prev_section ="", next_section = "", run_pseudo_haploid = FALSE);
            ## 
            rcpp_make_eMatGrid_t(eMatGrid_t2, eMatRead_t, 1, sampleReads, 1, nGrids, prev = 0, suppressOutput = 1, prev_section = "",next_section = "", run_fb_grid_offset = 0, TRUE, TRUE, maxEmissionMatrixDifference = maxEmissionMatrixDifference, rescale = TRUE)

            expect_equal(eMatGrid_t1, eMatGrid_t2)

}
