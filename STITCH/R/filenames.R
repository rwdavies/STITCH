## centralize the writing of the sampleRead filename
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

file_break_results <- function(tempdir, regionName) {
    return(file.path(
        tempdir,
        paste0("nbreaks.",regionName,".RData")
    ))
}

file_fbdStore <- function(tempdir, regionName, iteration) {
    return(file.path(
        tempdir,
        paste0("fbdStore.",regionName,".", iteration, ".RData")
    ))
}
