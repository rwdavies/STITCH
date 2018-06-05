#' @title Make STITCH command line interface
#' @param function_file to main STITCH function file
#' @param stitch_cli_file where output goes
make_STITCH_cli <- function(
    function_file,
    cli_output_file,
    integer_vectors = c("shuffleHaplotypeIterations", "splitReadIterations", "refillIterations"),
    character_vectors = c("reference_populations"),
    other_logical_params = NULL,
    other_integer_params = NULL,
    other_double_params = NULL,
    other_character_params = NULL,
    function_name = "STITCH"
) {

    a <- readLines(function_file, n = 200)
    params <- a[grep("param", a)]
    x <- grep("initialize_parameters", params)
    if (length(x) > 0)
        params <- params[-x]
    x <- grep("##", params)
    if (length(x) > 0)
        params <- params[-x]
    x <- sapply(strsplit(params, "#' @param "), function(x) x[2])
    x <- hack_for_1_2_5(x)
    param_names <- sapply(strsplit(x, " "), function(x) x[1])
    help_per_param <- sapply(strsplit(x, " "), function(x) paste(x[-1], collapse = " "))
    names(help_per_param) <- param_names


    ## manually specify types here
    logical_params <- c("outputInputInVCFFormat", "readAware", "regenerateInput", "keepInterimFiles", "keepTempDir", "generateInputOnly", "useSoftClippedBases", "regenerateInputWithDefaultValues", "plotHapSumDuringIterations", "plotAfterImputation", "save_sampleReadsInfo")
    integer_params <- c("K", "chrStart", "chrEnd", "regionStart", "regionEnd", "buffer", "niterations", "nCores", "Jmax", "pseudoHaploidModel", "switchModelIteration", "outputBlockSize", "inputBundleBlockSize", "reference_phred", "reference_iterations", "gridWindowSize")
    double_params <- c("nGen", "downsampleToCov", "downsampleFraction", "maxDifferenceBetweenReads", "maxEmissionMatrixDifference", "alphaMatThreshold", "emissionThreshold", "iSizeUpperLimit", "bqFilter", "expRate", "maxRate", "minRate", "downsampleSamples", "initial_min_hapProb", "initial_max_hapProb")
    character_params <- c("chr", "posfile", "outputdir", "tempdir", "bamlist", "cramlist", "reference", "genfile", "method", "shuffleHaplotypeIterations", "splitReadIterations", "originalRegionName", "environment", "restartIterations", "refillIterations", "downsampleSamplesKeepList", "subsetSNPsfile", "reference_haplotype_file", "reference_legend_file", "reference_sample_file", "reference_populations", "vcf_output_name", "output_filename")
    ## deprecated
    integer_params <- c(integer_params, "diploidModel")
    logical_params <- c(logical_params, "outputHaplotypeProbabilities")
    
    logical_params <- check_and_add_param(other_logical_params, logical_params, integer_params, double_params, character_params)
    integer_params <- check_and_add_param(other_integer_params, integer_params, logical_params, double_params, character_params)
    double_params <- check_and_add_param(other_double_params, double_params, logical_params, integer_params, character_params)
    character_params <- check_and_add_param(other_character_params, character_params, logical_params, integer_params, double_params)    
    
    param_type <- array(NA, length(param_names))
    param_type[match(logical_params, param_names)] <- "logical"
    param_type[match(integer_params, param_names)] <- "integer"
    param_type[match(double_params, param_names)] <- "double"
    param_type[match(character_params, param_names)] <- "character"
    if (sum(is.na(param_type)) > 0)
        stop(
            paste0(
                "Unassigned parameter types:",
                paste0(param_names[is.na(param_type)], collapse = ",")
            )
        )

    names(param_type) <- param_names

    ## further, note that some are in fact vectors, or can be vectors
    ## so down below, need to split intelligently

    ## hmm, not sure how to do this otherwise
    ## just manually specify numeric
    ## either nothing after term, in which case no default, or get default underneath
    defaults <- sapply(param_names, function(param) {
        ## if it has an equals, it has a default
        x <- paste0("    ", param, " = ")
        i <- which(
            substr(a, 1, nchar(param) + 7) == x
        )
        if (length(i) == 0) {
            return(NULL)
        } else {
            return(strsplit(strsplit(a[i], x)[[1]][2], ",")[[1]])
        }
    })

    cat(
        "#!/usr/bin/env Rscript\n",
        "\n",
        "if (!suppressPackageStartupMessages(require(\"optparse\")))\n",
        "    install.packages(\"optparse\", repos=\"http://cran.rstudio.com/\")\n",
        "\n",
        "option_list <- list(\n",
        sep = "", file = cli_output_file, append = FALSE
    )
    for(param in param_names) {
        default <- defaults[[param]]
        ori_default_length <- 1
        if (length(default) > 1) {
            ori_default_length <- length(default)
            default <- paste0(default, collapse = ",")
        }
        cat(
            "    make_option(\n",
            "        \"--", param, "\",\n",
            "        type = \"", param_type[param], "\",\n",
            "        help = \"", help_per_param[param],
            sep = "", file = cli_output_file, append = TRUE
        )
        if (is.null(default) == FALSE) {
            default_string <- default
            if (default == "\"\"") {
                default_string <- "\\\"\\\""
            }
            if (default == "\"diploid\"") {
                default_string <- "diploid"
            }
            if (default == "\"server\"") {
                default_string <- "server"
            }
            if (ori_default_length > 1) {
                default_string <- default
                default <- paste0("\"", default, "\"")
            }
            ## determine what it is - numeric or otherwise
            cat(
                " [default ", default_string, "] \"",
                ",\n",
                "        default = ", default,
                sep = "", file = cli_output_file, append = TRUE
            )
        } else {
            cat("\"", sep = "", file = cli_output_file, append = TRUE)
        }
        if (param != param_names[length(param_names)])
            cat("\n", "    ), ", "\n", sep = "", file = cli_output_file, append = TRUE)
    }
    cat("\n    )\n)\n", sep = "", file = cli_output_file, append = TRUE)
    ## why of why is parse_args giving warnings when it seems to be working fine
    cat("opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))\n", sep = "", file = cli_output_file, append = TRUE)
    cat("suppressPackageStartupMessages(library(STITCH))\n", sep = "", file = cli_output_file, append = TRUE)
    cat("Sys.setenv(PATH = paste0(Sys.getenv(\"PATH\"), \":\", getwd()))\n", sep = "", file = cli_output_file, append = TRUE)

    ## cat("print(opt)\n", sep = "", file = cli_output_file, append = TRUE)

    cat(function_name, "(", sep = "", file = cli_output_file, append = TRUE)

    for(param in param_names) {
        ## evaluate the ones that are vectors
        ## this screws up things like NULL
        if (param %in% integer_vectors | param %in% character_vectors) {
            cat(
                "\n    ", param, " = ",
                "eval(parse(text=opt$", param, "))",
                sep = "", file = cli_output_file, append = TRUE
            )
        } else {
            ## "as.", param_type[param], "(
            cat(
                "\n    ", param, " = ",
                "opt$", param,
                sep = "", file = cli_output_file, append = TRUE
            )
        }
        if (param != param_names[length(param_names)])
            cat(",", sep = "", file = cli_output_file, append = TRUE)
    }

    cat("\n)\n", sep = "", file = cli_output_file, append = TRUE)

}



check_param <- function(new_params, should_not_be_in_set) {
    w <- match(new_params, should_not_be_in_set)
    if (sum(is.na(w) == FALSE) > 0) {
        stop(paste0("bad assignment:", new_params[is.na(w) == FALSE]))
    }
}

check_and_add_param <- function(new_params, old_params, old_params1, old_params2, old_params3) {
    if (is.null(new_params) == FALSE) {
        ## 1 - check not in any other set
        ## 2 - add then uniqueify
        check_param(new_params, old_params1)
        check_param(new_params, old_params2)
        check_param(new_params, old_params3)
        old_params <- unique(c(new_params, old_params))
    }
    return(old_params)
}


hack_for_1_2_5 <- function(x) {
    if (length(x) == 61) {
        y <- grep("initial_min_hapProb", x[59])
        if (length(y) == 1) {
            if (y == 1) {
                x[59] <- substr(x[59], 2, 1000)
            }
        }
        y <- grep("initial_max_hapProb", x[60])
        if (length(y) == 1) {
            if (y == 1) {
                x[60] <- substr(x[60], 2, 1000)
            }
        }
    }
    return(x)
}
