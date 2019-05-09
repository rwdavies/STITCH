#' @export
file_sampleReads <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("sample.",iBam,".input.",regionName,".RData")
    ))
}

#' @export
file_sampleReadsInfo <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("sampleReadsInfo.",iBam,".input.",regionName,".RData")
    ))
}

file_referenceSampleReads <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("referenceSampleReads.",iBam,".input.",regionName,".RData")
    ))
}


file_sampleProbs <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("sampleProbs.",iBam,".input.",regionName,".RData")
    ))
}

file_alphaBetaBlocks <- function(
    dir,
    iBam,
    regionName
) {
    return(
        file.path(
            dir,
            paste0("alphaBetaBlocks.", iBam , ".", regionName, ".RData")
        )
    )
}


file_bundledSampleReads <- function(
  dir,
  start,
  end,
  regionName
) {
    return(
        file.path(
            dir,
            paste0("bundledSamples.", start, "-", end, ".", regionName, ".RData")
        )
    )
}


file_bundledReferenceSampleReads <- function(
  dir,
  start,
  end,
  regionName
) {
  return(file.path(
      dir,
      paste0("bundledReferenceSampleReads.", start, "-", end, ".", regionName, ".RData")
  ))
}


file_bundledSampleProbs <- function(
  dir,
  start,
  end,
  regionName
) {
    return(file.path(
        dir,
        paste0("bundledSampleProbs.", start, "-", end, ".", regionName, ".RData")
    ))
}

file_bundledAlphaBetaBlocks <- function(
    dir,
    start,
    end,
    regionName
) {
    return(
        file.path(
            dir,
            paste0("bundledAlphaBetaBlocks.", start, "-", end, ".", regionName, ".RData")
        )
    )
}


file_besthaps <- function(tempdir, iSample, regionName) {
    return(file.path(tempdir, paste0("sample.",iSample,".haps.",regionName,".RData")))
}

file_dosages <- function(
  dir,
  iBam,
  regionName,
  what = "dosage"
) {
    return(file.path(
        dir,
        paste0("sample.", iBam, ".input.", regionName, ".", what, ".RData")
    ))
}

file_haplotypes <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("sample.", iBam, ".input.",regionName,".haplotypes.RData")
    ))
}

file_break_results <- function(tempdir, regionName, iteration = NULL) {
    return(file.path(
        tempdir,
        paste0(
            "nbreaks.", regionName,
            c(".", "")[as.integer(is.null(iteration)) + 1],
            iteration, ".RData"
        )
    ))
}

file_fbdStore <- function(tempdir, regionName, iteration) {
    return(file.path(
        tempdir,
        paste0("fbdStore.",regionName,".", iteration, ".RData")
    ))
}

#' @export
file_date <- function(outputdir, regionName, what) {
    if (!(what %in% c("start", "startEM", "endEM", "end"))) {
        stop("bad date file selection")
    }
    file.path(outputdir, "RData", paste0(what, ".", regionName,".RData"))
}
