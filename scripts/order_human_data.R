## order human data by sampleReads[[2]]
human_matched_to_reference_datadir <- file.path(getwd(), "test-data/human_data_matched_to_reference/")
get_and_untar_if_tgz_file <- function(file) {
    if (file.exists(file) == FALSE) {
        if (operating_system == "mac") {
            system(paste0('curl ', ancillary_http, file, ' -o "', file, '"'))
        } else {
            system(paste0("wget ", ancillary_http, file))
        }
    }
    if (substr(file, nchar(file) - 3, nchar(file)) == ".tgz")
        system(paste0("tar -xzvf ", file))
}
ancillary_http <- "http://www.well.ox.ac.uk/~rwdavies/ancillary/"
setwd("/data/smew1/rdavies/stitch_development/STITCH_github_latest/STITCH")
dir.create("order_human_data")
setwd("order_human_data")
get_and_untar_if_tgz_file("STITCH_human_reference_example_2017_05_24.tgz")
setwd("input")
files <- dir()
file_bundle <- files[1]

for(file_bundle in files) {
    message(file_bundle)
    load(file_bundle)
    for(i_sample in which(names(bundledSampleReads) == "")){
        sampleReads <- bundledSampleReads[[i_sample]]
        central_snp <- sapply(sampleReads, function(x) x[[2]])
        sampleReads <- sampleReads[order(central_snp)]
        bundledSampleReads[[i_sample]] <- sampleReads
    }
    save(bundledSampleReads, file = file_bundle)
}

setwd("../")
unlink("STITCH_human_reference_example_2017_05_24.tgz")
system("tar -czvf STITCH_human_reference_example_2018_07_11.tgz *")
