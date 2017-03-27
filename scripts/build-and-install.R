#!/usr/bin/env Rscript

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}

clean_build <- TRUE
if (Sys.getenv("DIRTY_BUILD") == "TRUE")
    clean_build <- FALSE

## specify package
pkg <- "STITCH"
## Sys.setenv(SEQLIB_ROOT = file.path(getwd(), "SeqLib"))
root_dir <- getwd()

## documentation
devtools::document(pkg = pkg)

## make the tarball
print("Build")
if (clean_build == TRUE) {
    SeqLib_dir <- file.path(pkg, "src", "SeqLib")
    setwd(SeqLib_dir)
    system(paste0("make clean"))
    setwd(root_dir)
}
package_tarball <- devtools::build(pkg = pkg, manual = TRUE)
version <- unlist(
    strsplit(unlist(strsplit(basename(package_tarball), "STITCH_")), ".tar.gz")
)

## move tarball to releases
release_package_tarball <- file.path("releases", paste0("STITCH_", version, ".tar.gz"))
system(paste0("mv ", package_tarball, " ", release_package_tarball))
package_tarball <- release_package_tarball

## install from tarball
install.packages(package_tarball)

## build PDF
pdf_name <- file.path("releases", paste0("STITCH_", version, ".pdf"))
args = c(
    "CMD", "Rd2pdf", pkg, "--batch", "--force", "--no-preview",
    paste0("--output=", pdf_name)
)
system2("R", args)
