#' @title Sequencing To Imputation Through Constructing Haplotypes
#' @param chr What chromosome to run. Should match BAM headers
#' @param posfile Where to find file with positions to run. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>
#' @param K How many founder / mosaic haplotypes to use
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure
#' @param outputdir What output directory to use
#' @param tempdir What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/
#' @param bamlist Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai
#' @param cramlist Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into STITCH
#' @param sampleNames_file Optional, if not specified, sampleNames are taken from the SM tag in the header of the BAM / CRAM file. This argument is the path to file with sampleNames for samples. It is used directly to name samples in the order they appear in the bamlist / cramlist
#' 
#' @param reference Path to reference fasta used for making cram files. Only required if cramlist is defined
#' @param genfile Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header
#' @param method How to run imputation - either diploid, pseudoHaploid, or diploid-inbred. Please see main README for more information. All methods assume diploid samples. diploid is the most accurate but slowest, while pseudoHaploid may be advantageous for large sample sizes and K. diploid-inbred assumes all samples are inbred and invokes an internal haploid mathematical model but outputs diploid genotypes and probabilities
#' @param output_format one of bgvcf (i.e. bgziped VCF) or bgen (Layout = 2, CompressedSNPBlocks = 1)
#' @param B_bit_prob when using bgen, how many bits to use to store each double. Optiosn are 8, 16, 24 or 32
#' @param outputInputInVCFFormat Whether to output the input in vcf format
#' @param downsampleToCov What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage
#' @param downsampleFraction Downsample BAMs by choosing a fraction of reads to retain. Must be value 0<downsampleFraction<1
#' @param readAware Whether to run the algorithm is read aware mode. If false, then reads are split into new reads, one per SNP per read
#' @param chrStart When loading from BAM, some start position, before SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param chrEnd When loading from BAM, some end position, after SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop.
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param maxDifferenceBetweenReads How much of a difference to allow the reads to make in the forward backward probability calculation. For example, if P(read | state 1)=1 and P(read | state 2)=1e-6, re-scale so that their ratio is this value. This helps prevent any individual read as having too much of an influence on state changes, helping prevent against influence by false positive SNPs
#' @param maxEmissionMatrixDifference Similar to maxDifferenceBetweenReads, specifies ratio of how much larger the most probable state can be than the least probable state, but across all reads rather than for a single read. This helps to limit overflow in C++ calculations
#' @param alphaMatThreshold Minimum (maximum is 1 minus this) state switching into probabilities
#' @param emissionThreshold Emission probability bounds. emissionThreshold < P(alt read | state k) < (1-emissionThreshold)
#' @param iSizeUpperLimit Do not use reads with an insert size of more than this value
#' @param bqFilter Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used
#' @param niterations Number of EM iterations.
#' @param shuffleHaplotypeIterations Iterations on which to perform heuristic attempt to shuffle founder haplotypes for better fit. To disable set to NA.
#' @param splitReadIterations Iterations to try and split reads which may span recombination breakpoints for a better fit
#' @param nCores How many cores to use
#' @param expRate Expected recombination rate in cM/Mb
#' @param maxRate Maximum recomb rate cM/Mb
#' @param minRate Minimum recomb rate cM/Mb
#' @param Jmax Maximum number of SNPs on a read
#' @param regenerateInput Whether to regenerate input files. If this is FALSE, please using the same regionStart, regionEnd, buffer and posfile as you used to generate the input. Setting any of those to different values can cause the previous input data to be improperly interpreted. Please also see originalRegionName and regenerateInputWithDefaultValues
#' @param originalRegionName If regenerateInput is FALSE (i.e. using existing data), this is the name of the original region name (chr.regionStart.regionEnd). This is necessary to load past variables
#' @param keepInterimFiles Whether to keep interim parameter estimates
#' @param keepTempDir Whether to keep files in temporary directory
#' @param switchModelIteration Whether to switch from pseudoHaploid to diploid and at what iteration (NA for no switching)
#' @param generateInputOnly Whether to just generate input data then quit
#' @param restartIterations In pseudoHaploid method, which iterations to look for collapsed haplotype prnobabilities to resolve
#' @param refillIterations When to try and refill some of the less frequently used haplotypes
#' @param downsampleSamples What fraction of samples to retain. Useful for checking effect of N on imputation. Not meant for general use
#' @param downsampleSamplesKeepList When downsampling samples, specify a numeric list of samples to keep
#' @param subsetSNPsfile If input data has already been made for a region, then subset down to a new set of SNPs, as given by this file. Not meant for general use
#' @param useSoftClippedBases Whether to use (TRUE) or not use (FALSE) bases in soft clipped portions of reads
#' @param outputBlockSize How many samples to write out to disk at the same time when making temporary VCFs that are later pasted together at the end to make the final VCF. Smaller means lower RAM footprint, larger means faster write.
#' @param outputSNPBlockSize How many SNPs to write to disk at one time to reduce RAM usage when making VCFs
#' @param inputBundleBlockSize If NA, disable bundling of input files. If not NA, bundle together input files in sets of <= inputBundleBlockSize together
#' @param reference_haplotype_file Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)
#' @param reference_legend_file Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)
#' @param reference_sample_file Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)
#' @param reference_populations Vector with character populations to include from reference_sample_file e.g. CHB, CHS
#' @param reference_phred Phred scaled likelihood or an error of reference haplotype. Higher means more confidence in reference haplotype genotypes, lower means less confidence
#' @param reference_iterations When using reference haplotypes, how many iterations to use to train the starting data
#' @param reference_shuffleHaplotypeIterations When using reference haplotypes, how much shuffling to do to lead to better global fit
#' @param output_filename Override the default bgzip-VCF / bgen output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple STITCH runs are processing on the same region then they may over-write each others inputs and outputs
#' @param initial_min_hapProb Initial lower bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param initial_max_hapProb Initial upper bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param regenerateInputWithDefaultValues If regenerateInput is FALSE and the original input data was made using regionStart, regionEnd and buffer as default values, set this equal to TRUE
#' @param plotHapSumDuringIterations Boolean TRUE/FALSE about whether to make a plot that shows the relative number of individuals using each ancestral haplotype in each iteration
#' @param plot_shuffle_haplotype_attempts Boolean TRUE/FALSE about whether to make a plot that tries to show the selection of ancestral haplotypes to check for shuffling / flipping
#' @param plotAfterImputation Boolean TRUE/FALSE about whether to make plots after imputation has run (can be set to FALSE if this throws errors on systems without x11)
#' @param save_sampleReadsInfo Experimental. Boolean TRUE/FALSE about whether to save additional information about the reads that were extracted
#' @param gridWindowSize Whether to work on a grid where reads are binned into windows of this size (1 based, i.e. first bin is bases 1-gridWindowSize). This is particularly appropriate for very low coverage data (e.g. less than 0.2X) and can substantially speed up analyses
#' @param shuffle_bin_nSNPs Parameter that controls how to detect ancestral haplotypes that are shuffled during EM for possible re-setting. If set (not NULL), then break per-SNP (or per-grid) every this many SNPs / grids, and compare each to detect whether haplotypes either 1) are more likely to stay where they are or 2) switch from one haplotype to another. Note that only one of shuffle_bin_nSNPs or shuffle_bin_radius should be non-NULL
#' @param shuffle_bin_radius Parameter that controls how to detect ancestral haplotypes that are shuffled during EM for possible re-setting. If set (not NULL), then recombination rate is calculated around pairs of SNPs in window of twice this value, and those that exceed what should be the maximum (defined by nGen and maxRate) are checked for whether they are shuffled
#' @param keepSampleReadsInRAM Whether to (generally) keep sampleReads in RAM or store them in the temporary directory. STITCH is substantially faster if this is FALSE at the expense of RAM
#' @param useTempdirWhileWriting Whether to use temporary directory while writing output file (TRUE), or to keep result in RAM (FALSE). Using temporary directory is slower but uses less RAM
#' @param output_haplotype_dosages Whether to output ancestral haplotype dosages, i.e. the expected number of ancestral haplotypes carried by that sample at that locus
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
STITCH <- function(
    chr,
    nGen,
    posfile,
    K,
    outputdir,
    tempdir = NA,
    bamlist = "",
    cramlist = "",
    sampleNames_file = "",
    reference = "",
    genfile = "",
    method = "diploid",
    output_format = "bgvcf",
    B_bit_prob = 16,
    outputInputInVCFFormat = FALSE,
    downsampleToCov = 50,
    downsampleFraction = 1,
    readAware = TRUE,
    chrStart = NA,
    chrEnd = NA,
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    maxDifferenceBetweenReads = 1000,
    maxEmissionMatrixDifference = 1e10,
    alphaMatThreshold = 1e-4,
    emissionThreshold = 1e-4,
    iSizeUpperLimit = as.integer(600),
    bqFilter = as.integer(17),
    niterations = 40,
    shuffleHaplotypeIterations = c(4, 8, 12, 16),
    splitReadIterations = 25,
    nCores = 1,
    expRate = 0.5,
    maxRate = 100,
    minRate = 0.1,
    Jmax = 1000,
    regenerateInput = TRUE,
    originalRegionName = NA,
    keepInterimFiles = FALSE,
    keepTempDir = FALSE,
    outputHaplotypeProbabilities = FALSE,
    switchModelIteration = NA,
    generateInputOnly = FALSE,
    restartIterations = NA,
    refillIterations = c(6, 10, 14, 18),
    downsampleSamples = 1,
    downsampleSamplesKeepList = NA,
    subsetSNPsfile = NA,
    useSoftClippedBases = FALSE,
    outputBlockSize = 1000,
    outputSNPBlockSize = 10000,
    inputBundleBlockSize = NA,
    reference_haplotype_file = "",
    reference_legend_file = "",
    reference_sample_file = "",
    reference_populations = NA,
    reference_phred = 20,
    reference_iterations = 40,
    reference_shuffleHaplotypeIterations = c(4, 8, 12, 16),
    output_filename = NULL,
    initial_min_hapProb = 0.4,
    initial_max_hapProb = 0.6,
    regenerateInputWithDefaultValues = FALSE,
    plotHapSumDuringIterations = FALSE,
    plot_shuffle_haplotype_attempts = FALSE,
    plotAfterImputation = TRUE,
    save_sampleReadsInfo = FALSE,
    gridWindowSize = NA,
    shuffle_bin_nSNPs = NULL,
    shuffle_bin_radius = 5000,
    keepSampleReadsInRAM = FALSE,
    useTempdirWhileWriting = FALSE,
    output_haplotype_dosages = FALSE        
) {

    ## capture command line
    x <- as.list(environment())
    command_line <- paste0(
        "STITCH(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

    minimizeSwitchingIterations = NA

    ## #' @param pseudoHaploidModel How to model read probabilities in pseudo diploid model (shouldn't be changed)    
    pseudoHaploidModel <- 9 ## remove this from being settable
    environment <- "server" ## I'm assuming this is true anyway. more or less deprecated
    
    ##    K_subset = NA,
    ##    K_random = NA
    
    phasefile <- "" ## do not export for now
    ##  #' @param phasefile Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header

    outputHaplotypeProbabilities <- FALSE ## do not export for now
    ## #' @param outputHaplotypeProbabilities Whether to output haplotype probabilities in files

    ## if doing e.g. STITCH again
    ## add this to the VCF pieces
    ## do not put in tempdir as might be too small for outputs
    vcf.piece_unique <- paste0(toupper(letters[sample(10)]), collapse = "")

    ##
    ## specify output description
    ##
    regionName <- chr
    options(scipen = 999) ## Dangit rounding
    if(is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE)
        regionName <- paste0(chr, ".", regionStart,".", regionEnd)


    ##
    ## define some other variables - remove eventually
    ##
    outputInputGen <- FALSE # whether to output input in gen format
    startIterations <- 0 # how many start iterations to run

    windowSNPs <- 5000 # if doing the window start, how many SNPs to loop over
    if(method=="diploid")
        restartIterations <- NA


    ##
    ## external dependency checks
    ##
    check_program_dependency("rsync")
    check_program_dependency("bgzip")



    ##
    ##
    ## validate parameters
    ##
    ##
    validate_method(method)
    validate_chr(chr)
    validate_nGen(nGen)
    validate_nCores(nCores)
    validate_posfile(posfile)
    validate_genfile(genfile)
    validate_K(K)
    validate_K_subset(method, K, K_subset)
    validate_outputdir(outputdir)
    validate_tempdir(tempdir)
    validate_outputBlockSize(outputBlockSize)
    validate_downsampleFraction(downsampleFraction)
    validate_downsampleSamples(downsampleSamples)
    validate_plotHapSumDuringIterations(plotHapSumDuringIterations)
    validate_gridWindowSize(gridWindowSize)

    validate_output_format(output_format)    
    validate_output_filename(output_filename, output_format)
    validate_B_bit_prob(B_bit_prob, output_format)
    validate_output_haplotype_dosages(output_haplotype_dosages, output_format)    

    ## more involved checks
    validate_regionStart_regionEnd_and_buffer(regionStart, regionEnd, buffer)
    validate_bamlist_and_cramlist_for_input_generation(regenerateInput, originalRegionName, bamlist, cramlist, regionStart, regionEnd, buffer, regenerateInputWithDefaultValues, reference)
    validate_inputBundleBlockSize(inputBundleBlockSize, readAware)
    validate_reference_files(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations)
    validate_refillIterations(refillIterations, niterations)
    validate_shuffleHaplotypeIterations(shuffleHaplotypeIterations, niterations)
    validate_shuffleHaplotypeIterations(reference_shuffleHaplotypeIterations, reference_iterations)
    validate_hapProb(initial_min_hapProb, initial_max_hapProb)


    ##
    ##
    ## make necessary directory
    ##
    ##
    out <- initialize_directories(
        tempdir = tempdir,
        keepTempDir = keepTempDir,
        outputdir = outputdir
    )
    outputdir <- out$outputdir ## path.expand
    tempdir <- out$tempdir
    inputdir <- out$inputdir


    ##
    ## print out date
    ##
    print_message("Program start")
    if(generateInputOnly==FALSE) {
        date <- date()        
        file <- file.path(outputdir, "RData", paste0("start.", regionName ,".RData"))
        save(date, file = file)
    }


    ##
    ## load the positions, genotypes and phases
    ##
    out <- get_and_validate_pos_gen_and_phase(
        posfile = posfile,
        genfile = genfile,
        phasefile = phasefile,
        chr = chr,
        verbose = TRUE
    )
    pos <- out$pos
    gen <- out$gen
    phase <- out$phase
    nSNPs <- out$nSNPs
    L <- out$L




    ##
    ## shrink (if regionStart and regionEnd are NA)
    ##
    out <- shrink_region(
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        L = L,
        pos = pos,
        gen = gen,
        phase = phase
    )
    pos <- out$pos
    gen <- out$gen
    phase <- out$phase
    L <- out$L
    nSNPs <- out$nSNPs
    inRegionL <- out$inRegionL
    start_and_end_minus_buffer <- out$start_and_end_minus_buffer



    ##
    ## determine chrStart and chrEnd
    ##
    out <- initialize_chrStart_and_chrEnd(
        chrStart = chrStart,
        chrEnd = chrEnd,
        L = L,
        iSizeUpperLimit = iSizeUpperLimit
    )
    chrStart <- out$chrStart
    chrEnd <- out$chrEnd




    ##
    ##
    ## get sample names
    ##
    ##
    out <- get_sample_names(
        bamlist = bamlist,
        cramlist = cramlist,
        nCores = nCores,
        outputdir = outputdir,
        regionName = regionName,
        originalRegionName = originalRegionName,
        sampleNames_file = sampleNames_file
    )
    N <- out$N
    sampleNames <- out$sampleNames
    bam_files <- out$bam_files
    cram_files <- out$cram_files


    ##
    ## get matches to high coverage gen and phase
    ##
    out <- match_gen_and_phase_to_samples(
        sampleNames = sampleNames,
        gen = gen,
        phase = phase
    )
    highCovInLow <- out$highCovInLow
    samples_with_phase <- out$samples_with_phase




    ##
    ## if inputBundleBlockSize is not NA
    ## get bundling matrix for input files
    ##
    bundling_info <- get_bundling_position_information(
        N = N,
        nCores = nCores,
        blockSize = inputBundleBlockSize
    )

    ##
    ## either generate the data, or load it from before
    ##
    generate_or_refactor_input(regenerateInput = regenerateInput, bundling_info = bundling_info, L = L, pos = pos, nSNPs = nSNPs, bam_files = bam_files, cram_files = cram_files, reference = reference, iSizeUpperLimit = iSizeUpperLimit, bqFilter = bqFilter, chr = chr, outputdir = outputdir, N = N, downsampleToCov = downsampleToCov, sampleNames = sampleNames, inputdir = inputdir, useSoftClippedBases = useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, generateInputOnly = generateInputOnly, environment = environment, nCores = nCores, save_sampleReadsInfo = save_sampleReadsInfo)
    ## if only generating input data, we are done
    if (generateInputOnly)
        return(NULL)





    ##
    ## if necessary, shrink BAMs, but only if regenerateInput = FALSE
    ##
    shrinkReads(N = N, nCores = nCores, originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL, environment = environment, regenerateInput = regenerateInput, inputBundleBlockSize = inputBundleBlockSize)



    ##
    ## if we want, re-intersect with new set of positions
    ##
    if(is.na(subsetSNPsfile)==FALSE & subsetSNPsfile!="NA") {
        print_message(paste0("Subsetting SNPs from file ", subsetSNPsfile))
        ## load in positions
        out <- subsetSNPsFunction(N=N,subsetSNPsfile=subsetSNPsfile,regionName=regionName,tempdir=tempdir,L=L,environment=environment,nCores=nCores,outputdir=outputdir)
        print_message(paste0("Back to main ", subsetSNPsfile))
        keep <- out$keep
        T <- as.integer(sum(keep))
        L <- L[keep]
        pos <- pos[keep,]
        gen <- gen[keep,]
        print_message(paste0("Done subsetting section ",subsetSNPsfile))
        ## save(nSNPs , L, keep, pos, gen, file = paste0(outputdir, "/debug/files.RData"))
    }



    ##
    ## downsample the number of samples
    ##
    if(downsampleSamples < 1) {
        out <- downsample_the_samples(
            downsampleSamplesKeepList = downsampleSamplesKeepList,
            downsampleSamples = downsampleSamples,
            highCovInLow = highCovInLow,
            nCores = nCores,
            blockSize = blockSize,
            N = N,
            inputBundleBlockSize = inputBundleBlockSize,
            bundling_info = bundling_info,
            sampleNames = sampleNames,
            regionName = regionName,
            tempdir = tempdir,
            outputdir = outputdir,
            environment = environment
        )
        N <- out$N
        sampleNames <- out$sampleNames
        bundling_info <- out$bundling_info
        highCovInLow <- out$highCovInLow
    }





    ##
    ## downsample to percent
    ##
    if(downsampleFraction < 1) {
        downsampleToFraction(N=N,nCores=nCores,downsampleFraction=downsampleFraction,regionName=regionName,tempdir=tempdir, bundling_info = bundling_info)
    }

    ##
    ## at this point (could be earlier), keep sampleReads in RAM if so desired
    ##
    if (keepSampleReadsInRAM) {
        allSampleReads <- load_all_sampleReads_into_memory(
            N = N,
            nCores = nCores,
            tempdir = tempdir,
            regionName = regionName,
            bundling_info = bundling_info
        )
    } else {
        allSampleReads <- NULL
    }


    ##
    ## build alleleCount
    ##
    alleleCount <- buildAlleleCount(
        nSNPs = nSNPs,
        N = N,
        nCores = nCores,
        regionName = regionName,
        tempdir = tempdir,
        bundling_info = bundling_info,
        allSampleReads = allSampleReads
    )


    ##
    ## toggle read aware
    ##
    if (readAware == FALSE) {
        allSampleReads <- split_reads_completely(
            N = N,
            nCores = nCores,
            tempdir = tempdir,
            regionName = regionName,
            bundling_info = bundling_info,
            allSampleReads = allSampleReads
        )
    }

    ##
    ## if we're so inclined, output input in VCF format, then exit
    ##
    if (outputInputInVCFFormat == TRUE) {
        out <- outputInputInVCFFunction(outputdir = outputdir, pos = pos, nSNPs = nSNPs, tempdir = tempdir, N = N, nCores = nCores, regionName = regionName, environment = environment, sampleNames = sampleNames, outputBlockSize = outputBlockSize, bundling_info = bundling_info, output_filename = output_filename, vcf.piece_unique = vcf.piece_unique, output_format = output_format, allSampleReads = allSampleReads)
        return(NULL)
    }


    ## set up grid
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    snps_in_grid_1_based <- out$snps_in_grid_1_based

    ## determine output regions
    blocks_for_output <- determine_snp_and_grid_blocks_for_output(
        grid = grid,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        outputSNPBlockSize = outputSNPBlockSize
    )
    

    ##
    ## initialize variables
    ##
    out <- initialize_parameters(reference_haplotype_file = reference_haplotype_file, reference_legend_file = reference_legend_file, reference_sample_file = reference_sample_file, reference_populations = reference_populations, reference_phred = reference_phred, reference_iterations = reference_iterations, nSNPs = nSNPs, K = K, L = L, pos = pos, inputBundleBlockSize = inputBundleBlockSize, nCores = nCores, regionName = regionName, alleleCount = alleleCount, startIterations = startIterations, windowSNPs = windowSNPs, expRate = expRate, nGen = nGen, tempdir = tempdir, outputdir = outputdir, pseudoHaploidModel = pseudoHaploidModel, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, minRate = minRate, maxRate = maxRate, regionStart = regionStart, regionEnd = regionEnd, buffer = buffer, niterations = niterations, grid = grid, grid_distances = grid_distances, nGrids = nGrids, reference_shuffleHaplotypeIterations = reference_shuffleHaplotypeIterations, L_grid = L_grid, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, shuffle_bin_radius = shuffle_bin_radius, snps_in_grid_1_based = snps_in_grid_1_based)
    sigmaCurrent <- out$sigmaCurrent
    eHapsCurrent <- out$eHapsCurrent
    alphaMatCurrent <- out$alphaMatCurrent
    ##    hapSumCurrent <- out$hapSumCurrent ## not needed?
    priorCurrent <- out$priorCurrent
    reference_panel_SNPs <- out$reference_panel_SNPs




    ##
    ## snap reads to grid (or not)
    ##
    if (is.na(gridWindowSize) == FALSE) {
        allSampleReads <- snap_reads_to_grid(
            N = N,
            nCores = nCores,
            regionName = regionName,
            tempdir = tempdir,
            bundling_info = bundling_info,
            grid = grid,
            downsampleToCov = downsampleToCov,
            sampleNames = sampleNames,
            allSampleReads = allSampleReads
        )
    }

    ## comes after snap_reads_to_grid as makes srp
    if (method == "pseudoHaploid") {
        out <- initialize_readProbs(
            N = N,
            nCores = nCores,
            pseudoHaploidModel = pseudoHaploidModel,
            tempdir = tempdir,
            bundling_info = bundling_info,
            regionName = regionName,
            initial_min_hapProb = initial_min_hapProb,
            initial_max_hapProb = initial_max_hapProb,
            allSampleReads = allSampleReads
        )
    }


    

    ##
    ## run EM algorithm here
    ##
    print_message("Begin EM")
    date <- date()
    save(date, file = file.path(outputdir, "RData", paste0("startEM.", regionName, ".RData")))
    print_message(paste0("Number of samples: ", N))
    print_message(paste0("Number of SNPs: ", nSNPs))
    if (nGrids != nSNPs) {
        print_message(paste0("Number of grids: ", nGrids))
    }

    eHapsCurrent_t <- t(eHapsCurrent) ## do this once, eventually remove entirely
    alphaMatCurrent_t <- t(alphaMatCurrent)
    
    
    for(iteration in 1:niterations) {
        ##
        ## switch methods if appropriate
        ##
        if(is.na(switchModelIteration)==FALSE) {
            if (iteration == switchModelIteration) {
                method <- "diploid"
                print_message(paste0("Switching from pseudoHaploid to diploid - iteration ",iteration))
            }
        }
        ##
        ## fork out and get results
        ##
        results <- completeSampleIteration(N=N,tempdir=tempdir,chr=chr,K=K,K_subset = K_subset, K_random = K_random, nSNPs = nSNPs,nGrids = nGrids,priorCurrent=priorCurrent,eHapsCurrent_t=eHapsCurrent_t,alphaMatCurrent_t=alphaMatCurrent_t,sigmaCurrent=sigmaCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,maxEmissionMatrixDifference = maxEmissionMatrixDifference,Jmax=Jmax,highCovInLow=highCovInLow,iteration=iteration,method=method,expRate=expRate,minRate=minRate,maxRate=maxRate,niterations=niterations,splitReadIterations=splitReadIterations,shuffleHaplotypeIterations=shuffleHaplotypeIterations,nCores=nCores,L=L,nGen=nGen,emissionThreshold=emissionThreshold,alphaMatThreshold=alphaMatThreshold,gen=gen,outputdir=outputdir,environment=environment,pseudoHaploidModel=pseudoHaploidModel,outputHaplotypeProbabilities=outputHaplotypeProbabilities,switchModelIteration=switchModelIteration,regionName=regionName,restartIterations=restartIterations,refillIterations=refillIterations,outputBlockSize=outputBlockSize, bundling_info = bundling_info, alleleCount = alleleCount, phase = phase, samples_with_phase = samples_with_phase, vcf.piece_unique = vcf.piece_unique, grid = grid, grid_distances = grid_distances, L_grid = L_grid, B_bit_prob = B_bit_prob, nSNPsInRegion = nSNPsInRegion, start_and_end_minus_buffer = start_and_end_minus_buffer, shuffle_bin_nSNPs = shuffle_bin_nSNPs, shuffle_bin_radius = shuffle_bin_radius, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, blocks_for_output = blocks_for_output, allSampleReads = allSampleReads, snps_in_grid_1_based = snps_in_grid_1_based, minimizeSwitchingIterations = minimizeSwitchingIterations, useTempdirWhileWriting = useTempdirWhileWriting)
        ##
        if (iteration == niterations) {        
            allAlphaBetaBlocks <- results$allAlphaBetaBlocks
            ## this is a "result" for this iteration
            hapSumCurrent_t <- results$hapSum_t
        } else {
            ##
            ## save previous results to disk
            ##
            eHapsFuture_t <- results$eHapsFuture_t
            alphaMatFuture_t <- results$alphaMatFuture_t
            sigmaFuture <- results$sigmaFuture
            priorFuture <- results$priorFuture
            hapSumCurrent_t <- results$hapSum_t
            allAlphaBetaBlocks <- results$allAlphaBetaBlocks
            rm(results)
            ## save everything if needed
            if (keepInterimFiles) {
                save(
                    hapSumCurrent_t, eHapsCurrent_t,alphaMatCurrent_t,sigmaCurrent,priorCurrent,eHapsFuture_t,alphaMatFuture_t ,sigmaFuture,priorFuture,
                    file = file.path(outputdir, "RData", paste0("interim.",regionName,".iteration",iteration,".RData"))
                )
            }
            ## plot interim plots to understand performance better
            if (plotHapSumDuringIterations == TRUE) {
                interim_plotter(
                    outputdir = outputdir,
                    regionName = regionName,
                    iteration = iteration,
                    L_grid = L_grid,
                    hapSumCurrent_t = hapSumCurrent_t,
                    alphaMatCurrent_t = alphaMatCurrent_t,
                    sigmaCurrent = sigmaCurrent,
                    N = N
                )                 
            }
            ## perform switchover here
            eHapsCurrent_t <- eHapsFuture_t
            alphaMatCurrent_t <- alphaMatFuture_t
            sigmaCurrent <- sigmaFuture
            priorCurrent <- priorFuture
            rm(eHapsFuture_t, alphaMatFuture_t, sigmaFuture, priorFuture)
        }
        ## GC a bit if big enough data
        if ((nSNPs > 3000) | (N > 3000)) {
            gc(reset = TRUE); gc(reset = TRUE)
        }
    }

    ##
    ## build final output
    ##
    out <- make_and_write_output_file(
        output_filename = output_filename,
        outputdir = outputdir,
        regionName = regionName,
        output_format = output_format,
        blocks_for_output = blocks_for_output,
        allAlphaBetaBlocks = allAlphaBetaBlocks,
        reference_panel_SNPs = reference_panel_SNPs,
        priorCurrent = priorCurrent,
        sigmaCurrent = sigmaCurrent,
        alphaMatCurrent_t = alphaMatCurrent_t,
        eHapsCurrent_t = eHapsCurrent_t,
        N = N,
        method = method,
        sampleNames = sampleNames,
        nSNPs = nSNPs,
        nCores = nCores,
        B_bit_prob = B_bit_prob,
        bundling_info = bundling_info,
        tempdir = tempdir,
        grid = grid,
        nGrids = nGrids,
        alleleCount = alleleCount,
        pos = pos,
        K = K,
        highCovInLow = highCovInLow,
        start_and_end_minus_buffer = start_and_end_minus_buffer,
        allSampleReads = allSampleReads,
        niterations = niterations,
        maxEmissionMatrixDifference = maxEmissionMatrixDifference,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = Jmax,
        useTempdirWhileWriting = useTempdirWhileWriting,
        output_haplotype_dosages = output_haplotype_dosages
    )
    
    gen_imp <- out$gen_imp
    estimatedAlleleFrequency <- out$estimatedAlleleFrequency
    info <- out$info
    hwe <- out$hwe
    hweCount <- out$hweCount
    passQC <- info > 0.4 & hwe > 1e-6 # one interpretation of passQC

    ##
    ## now, if there was a buffer, remove the unwanted SNPs from what's saved to disk
    ## not sure how useful it is to save without buffer?
    ##
    if (is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE) {
        ##
        ## first, save (with buffer) to disk
        ##
        save(
            eHapsCurrent_t, alphaMatCurrent_t,
            sigmaCurrent, priorCurrent,
            hapSumCurrent_t,
            alleleCount,
            estimatedAlleleFrequency,
            pos, gen, L,
            nCores, N, nSNPs, chr, K, niterations, nGen,
            gen_imp,
            grid, grid_distances, L_grid, nGrids, snps_in_grid_1_based,
            inRegionL, start_and_end_minus_buffer,
            reference_panel_SNPs,
            info, hwe, passQC, hweCount,
            file = file.path(
                outputdir, "RData", paste0("EM.all.", regionName, ".withBuffer.RData")
            )
        )
        ##
        ## next, remove buffer
        ##
        out <- remove_buffer_from_variables(L = L,  regionStart = regionStart, regionEnd = regionEnd, pos = pos, gen = gen, phase = phase, alleleCount =  alleleCount, highCovInLow = highCovInLow, gen_imp = gen_imp, alphaMatCurrent_t = alphaMatCurrent_t, sigmaCurrent = sigmaCurrent, eHapsCurrent_t = eHapsCurrent_t, hapSumCurrent_t = hapSumCurrent_t, grid = grid, grid_distances = grid_distances, L_grid = L_grid, nGrids = nGrids, gridWindowSize = gridWindowSize, hwe = hwe, hweCount = hweCount, info = info, passQC = passQC, estimatedAlleleFrequency = estimatedAlleleFrequency)
        pos <- out$pos
        gen <- out$gen
        phase <- out$phase
        gen_imp <- out$gen_imp
        alleleCount <- out$alleleCount
        L <- out$L
        nSNPs <- out$nSNPs
        alphaMatCurrent_t <- out$alphaMatCurrent_t
        sigmaCurrent <- out$sigmaCurrent
        eHapsCurrent_t <- out$eHapsCurrent_t
        hapSumCurrent_t <- out$hapSumCurrent_t
        priorCurrent <- out$priorCurrent
        grid <- out$grid
        grid_distances <- out$grid_distances
        L_grid <- out$L_grid
        nGrids <- out$nGrids
        snps_in_grid_1_based <- out$snps_in_grid_1_based        
        reference_panel_SNPs <- out$reference_panel_SNPs
        hwe <- out$hwe
        hweCount <- out$hweCount
        info <- out$info
        estimatedAlleleFrequency <- out$estimatedAlleleFrequency
        passQC <- out$passQC
    }



    ##
    ## save important files from EM output
    ##
    print_message("Save RData objects to disk")
    save(
        eHapsCurrent_t, alphaMatCurrent_t,
        sigmaCurrent, priorCurrent,
        hapSumCurrent_t,
        alleleCount,
        estimatedAlleleFrequency,
        pos, gen, L,
        nCores, N, nSNPs, chr, K, niterations, nGen,
        gen_imp,
        grid, grid_distances, L_grid, nGrids, snps_in_grid_1_based,
        inRegionL, start_and_end_minus_buffer,
        reference_panel_SNPs,
        info, hwe, passQC, hweCount,
        file = file.path(
            outputdir, "RData", paste0("EM.all.", regionName, ".RData")
        )
    )



    if (plotAfterImputation) {
        ##
        ## do some plotting
        ##
        print_message("Make metrics plot")
        plotMetricsForPostImputationQC(
            iSample = NULL, highCovList = NULL, gen = gen, gen_imp = gen_imp,
            alleleCount = alleleCount, chr = chr, L = L,
            estimatedAlleleFrequency = estimatedAlleleFrequency, info = info,
            outputdir = outputdir, hwe = hwe, regionName = regionName
        )
        ##
        ## plot hapProbs sum - should tell what states are being used
        ##
        print_message("Make HapSum plot")
        outname <- file.path(outputdir, "plots", paste0("hapProbs.run.",regionName,".jpg"))
        plotHapSumCurrent_t(
            outname = outname,
            L_grid = L_grid,
            K = K,
            hapSumCurrent_t = hapSumCurrent_t,
            nGrids = nGrids,
            N = N
        )
        ##
        ## plot estimated AF against real 0.1X pileups (get r2 as well)
        ##
        if(sum(passQC)>1) {
            print_message("Make estimated against real")
            plotEstimatedAgainstReal(
                outputdir = outputdir,alleleCount=alleleCount,
                estimatedAlleleFrequency=estimatedAlleleFrequency,
                which=passQC,chr=chr,regionName=regionName
            )
        }
    }




    ##
    ## clean up and terminate
    ##
    print_message("Clean up and end")
    warnings()
    ## remove things from tempdir
    if (keepTempDir == FALSE) {
        unlink(tempdir, recursive = TRUE)
    }

    date <- date()
    print_message("Program done")
    save(date, file = file.path(outputdir, "RData", paste0("end.",regionName,".RData")))

    return(NULL)

}




#' @export
initialize_directories <- function(
    tempdir,
    keepTempDir,
    outputdir,
    STITCH_again = FALSE,
    keep_generated_input = FALSE
) {
    ## expand to remove tilda's, weird problems sometimes
    outputdir <- path.expand(outputdir)
    ## not ideal - overriding input parameter
    ## also not ideal - same name as function
    if (is.na(tempdir)) {
        tempdir <- tempdir()
    }
    tempdir <- path.expand(tempdir)
    tempdir <- file.path(
        tempdir,
        paste(
            toupper(letters[sample(26, 10, replace = TRUE)]), collapse=""
        )
    )
    if (keepTempDir == TRUE) {
        print_message(paste0("tempdir = ", tempdir))
    }
    dirs_to_create <- c(
        tempdir,
        outputdir,
        file.path(outputdir, "RData"),
        file.path(outputdir, "debug"),
        file.path(outputdir, "plots")
    )
    if (STITCH_again) {
        if (keep_generated_input) {
            inputdir <- file.path(outputdir, "input_again")
        } else {
            inputdir <- tempfile(tmpdir = tempdir)
        }
    } else {
        inputdir <- file.path(outputdir, "input")
    }
    dirs_to_create <- c(dirs_to_create, inputdir)
    for(dir in dirs_to_create) {
        dir.create(dir, showWarnings = FALSE)
        ## throw error if there was a failure to create the folder
        ## indicates some file system problems, e.g. slow NFS
        if (dir.exists(dir) == FALSE) {
            stop(paste0(
                "Unable to make the required directory ", dir, " while running. ",
                "You can try re-starting STITCH, but if the problem consists, please contact your system administrator."
            ))
        }
    }
    return(
        list(
            outputdir = outputdir,
            tempdir = tempdir,
            inputdir = inputdir
        )
    )
}

#' @export
initialize_chrStart_and_chrEnd <- function(chrStart, chrEnd, L, iSizeUpperLimit) {
    if ((as.integer(is.na(chrStart)) + as.integer(is.na(chrEnd))) == 1)
        stop("Please set either both or neither of chrStart and chrEnd to NA")
    if (is.na(chrStart) & is.na(chrEnd)) {
        chrStart <- min(L) - iSizeUpperLimit
        chrEnd <- max(L) + iSizeUpperLimit
    }
    return(
        list(
            chrStart = chrStart,
            chrEnd = chrEnd
        )
    )
}


## throw an error if dependencies are not installed
check_program_dependency <- function(program) {
    check <- Sys.which(program)
    if (check == "")
        stop(paste0(
            "The program ", program, " is not available in the PATH. ",
            "STITCH requires ", program, " to function. ",
            "Please make ", program, " available from the PATH"
        ))
    return(NULL)
}


#' @export
remove_buffer_from_variables <- function(
    L,
    regionStart,
    regionEnd,
    pos = NULL,
    gen = NULL,
    phase = NULL,
    alleleCount = NULL,
    highCovInLow = NULL,
    gen_imp = NULL,
    alphaMatCurrent_t = NULL,
    sigmaCurrent = NULL,
    eHapsCurrent_t = NULL,
    hapSumCurrent_t = NULL,
    info = NULL,
    hwe = NULL,
    estimatedAlleleFrequency = NULL,
    phase_entropy = NULL,
    strandedness = NULL,
    eHapsUpdate_numer = NULL,
    eHapsUpdate_denom = NULL,
    grid = NULL,
    grid_distances = NULL,
    L_grid = NULL,
    nGrids = NULL,
    gridWindowSize = NULL,
    reference_panel_SNPs = NULL,
    hweCount = NULL,
    passQC = NULL,
    verbose = TRUE
) {
    if (verbose)
        print_message("Remove buffer region from output")
    ## determine where the region is
    inRegion2 <- regionStart <= L & L <= regionEnd
    inRegion2L <- which(inRegion2)
    ## pos, gen, L
    pos <- pos[inRegion2, ]
    gen <- gen[inRegion2, ]
    phase <- phase[inRegion2, , , drop = FALSE]
    if ( length(highCovInLow) > 0 )
        gen_imp <- gen_imp[inRegion2, ]
    alleleCount <- alleleCount[inRegion2, ]
    L <- L[inRegion2]
    nSNPs <- as.integer(nrow(pos))
    reference_panel_SNPs <- reference_panel_SNPs[inRegion2]
    ## hmm, there is no perfect answer here
    ## for now, save everything that intersected
    to_keep <- unique(grid[inRegion2L]) + 1
    hapSumCurrent_t <- hapSumCurrent_t[, to_keep, drop = FALSE]
    ## now for these there is one fewer point    
    to_keep <- to_keep[-length(to_keep)]
    alphaMatCurrent_t <- alphaMatCurrent_t[, to_keep, drop = FALSE]
    sigmaCurrent <- sigmaCurrent[to_keep]
    ##
    ##
    eHapsCurrent_t <- eHapsCurrent_t[, inRegion2]
    priorCurrent <- hapSumCurrent_t[, 1] / sum(hapSumCurrent_t[, 1])
    ##
    hwe <- hwe[inRegion2]
    hweCount <- hweCount[inRegion2, , drop = FALSE]    
    info <- info[inRegion2]
    passQC <- passQC[inRegion2]
    estimatedAlleleFrequency <- estimatedAlleleFrequency[inRegion2]
    ## phasing stuff
    phase_entropy <- phase_entropy[inRegion2]
    strandedness <- strandedness[inRegion2, , drop = FALSE]
    eHapsUpdate_numer <- eHapsUpdate_numer[inRegion2, , drop = FALSE]
    eHapsUpdate_denom <- eHapsUpdate_denom[inRegion2, , drop = FALSE]
    ## grid stuff
    out <- assign_positions_to_grid(
        L = L,
        gridWindowSize = gridWindowSize
    )
    grid <- out$grid
    grid_distances <- out$grid_distances
    L_grid <- out$L_grid
    nGrids <- out$nGrids
    snps_in_grid_1_based <- out$snps_in_grid_1_based
    return(
        list(
            inRegion2 = inRegion2,
            pos = pos,
            gen = gen,
            phase = phase,
            gen_imp = gen_imp,
            alleleCount = alleleCount,
            L = L,
            nSNPs = nSNPs,
            alphaMatCurrent_t = alphaMatCurrent_t,
            sigmaCurrent = sigmaCurrent,
            eHapsCurrent_t = eHapsCurrent_t,
            hapSumCurrent_t = hapSumCurrent_t,
            priorCurrent = priorCurrent,
            info = info,
            hwe = hwe,
            hweCount = hweCount,
            passQC = passQC,
            estimatedAlleleFrequency = estimatedAlleleFrequency,
            phase_entropy = phase_entropy,
            strandedness = strandedness,
            eHapsUpdate_numer = eHapsUpdate_numer,
            eHapsUpdate_denom = eHapsUpdate_denom,
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based,
            reference_panel_SNPs = reference_panel_SNPs
        )
    )
}


#' @export
match_gen_and_phase_to_samples <- function(
    sampleNames,
    gen = NULL,
    phase = NULL
) {
    highCovInLow <- NULL
    if (is.null(gen) == FALSE) {
        highCovInLow <- match(colnames(gen), sampleNames)
        x <- is.na(highCovInLow)
        if (sum(x) > 0)
            stop(paste0(
                "The following gen file samples could not ",
                "be matched to a sample specified from bamlist:",
                paste(colnames(gen)[x], collapse = ", ")
            ))
    }
    ## phase
    samples_with_phase <- NULL
    if (is.null(phase) == FALSE) {
        samples_with_phase <- match(dimnames(phase)[[2]], sampleNames)
        x <- is.na(samples_with_phase)
        if (sum(x) > 0)
            stop(paste0(
                "The following phase file samples could not ",
                "be matched to a sample specified from bamlist:",
                paste(dimnames(phase)[[2]], collapse = ", ")
            ))
    }
    return(
        list(
            highCovInLow = highCovInLow,
            samples_with_phase = samples_with_phase
        )
    )
}






#' @export
get_sample_names <- function(
    bamlist = "",
    cramlist = "",
    nCores = 1,
    outputdir = NULL,
    regionName = NULL,
    originalRegionName = NULL,
    save = TRUE,
    verbose = TRUE,
    sampleNames_file = ""
) {

    bam_files <- NULL
    cram_files <- NULL
    if (bamlist != "" | cramlist != "") {
        if (cramlist != "") {
            files <- as.character(readLines(cramlist))
            cram_files <- files
            file_type <- "CRAM"
        } else {
            files <- as.character(readLines(bamlist))
            bam_files <- files
            file_type <- "BAM"
        }
        if (sampleNames_file == "") {
            sampleNames <- get_sample_names_from_bam_or_cram_files(
                files,
                nCores,
                file_type,
                verbose
            )
        } else {
            if (file.exists(sampleNames_file) == FALSE) {
                stop(paste0("Cannot find sampleNames_file:", sampleNames_file))
            }
            sampleNames <- read.table(sampleNames_file, stringsAsFactors = FALSE)[, 1]
            if (
            (length(sampleNames) != length(bam_files)) &
            (length(sampleNames) != length(cram_files))
            ) {
                stop(paste0("There are ", length(sampleNames), " sample names according to the sampleNames_file but there are ", max(c(length(bam_files), length(cram_files))), " according to the bam / cram list"))
            }
        }

        if(length(unique(sampleNames))!=length(sampleNames)) {
            stop("There are repeat sample names")
        }

        ## save sample names
        if (save) {
            save(
                sampleNames,
                file = file.path(
                    outputdir, "RData", paste0("sampleNames.", regionName, ".RData"))
            )
        }

    } else {
        ## otherwise, we won't be re-doing vinput anyway!
        ## load up sampleNames, fileHeaders
        file <- file.path(
            outputdir, "RData", paste0("sampleNames.", originalRegionName, ".RData")
        )
        if (file.exists(file) == TRUE) {
            load(file)
        } else {
            stop(paste0("Cannot find file:", file))
        }
    }
    N <- as.integer(length(sampleNames)) # number of samples

    return(
        list(
            sampleNames = sampleNames,
            N = N,
            bam_files = bam_files,
            cram_files = cram_files
        )
    )

}



## check regionStart, regionEnd, shrink things if appropriate
## so first, we start with posfile with "pos"
## then we define a sub-region from shrinking the "pos" to regionStart - buffer <= regionEnd + buffer
## then we define another sub-region that includes vs excludes the buffer
## AFTER SHRINKING
## defined by regionStart, regionEnd and buffer: pos, L, inRegionL (1-based from larger pos, what was kept)
## defined by regionStart, regionEnd: start_and_end_minus_buffer (1-based what to keep from current SNPs to remove the buffer)

#' @export
shrink_region <- function(
    regionStart,
    regionEnd,
    buffer,
    L,
    pos = NULL,
    gen = NULL,
    phase = NULL,
    validate = TRUE
) {
    if ( is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE ) {
        ## need to shrink pos, L, gen
        ## then, need to fix input
        inRegion <- ((regionStart-buffer) <= L) & (L <= (regionEnd + buffer))
        w <- which(inRegion)
        inRegionL <- c(head(w, 1), tail(w, 1)) ## 1-based, first SNP in buffer, last SNP outside buffer
        L <- L[inRegion]
        pos <- pos[inRegion, ]
        gen <- gen[inRegion, , drop = FALSE]
        phase <- phase[inRegion, , , drop = FALSE]
        nSNPs <- as.integer(nrow(pos))
        ## now, compare shrunk L to buffer region
        inRegion <- L >= (regionStart) & L <= (regionEnd)
        w <- which(inRegion)
        start_and_end_minus_buffer <- c(head(w, 1), tail(w, 1)) ## first SNP
        if (validate) {
            validate_region_to_impute_when_using_regionStart(L, regionStart, regionEnd, buffer)
        }
    } else {
        nSNPs <- length(L) 
        inRegionL <- c(1, nSNPs)
        start_and_end_minus_buffer <- c(1, nSNPs)
        nSNPsInRegionMinusBuffer <- nSNPs
    }
    ## print off number of SNPs in regions
    ## mostly useful to remember how to define
    ## nSNPsInRegionMinusBuffer <- start_and_end_minus_buffer[2] - start_and_end_minus_buffer[1] + 1    
    return(
        list(
            pos = pos,
            gen = gen,
            phase = phase,
            L = L,
            nSNPs = nSNPs,
            inRegionL = inRegionL,
            start_and_end_minus_buffer = start_and_end_minus_buffer
        )
    )
}






initialize_parameters <- function(
    reference_haplotype_file,
    reference_legend_file,
    reference_sample_file,
    reference_populations,
    reference_phred,
    reference_iterations,
    nSNPs,
    K,
    L,
    pos,
    inputBundleBlockSize,
    nCores,
    regionName,
    alleleCount,
    startIterations,
    windowSNPs,
    expRate,
    nGen,
    tempdir,
    outputdir,
    pseudoHaploidModel,
    emissionThreshold,
    alphaMatThreshold,
    minRate,
    maxRate,
    regionStart,
    regionEnd,
    buffer,
    niterations,
    grid,
    grid_distances,
    nGrids,
    reference_shuffleHaplotypeIterations,
    L_grid,
    plot_shuffle_haplotype_attempts,
    shuffle_bin_radius,
    snps_in_grid_1_based
) {

    print_message("Begin parameter initialization")
    
    ## default values
    eHapsCurrent <- matrix(runif(nSNPs * K), ncol = K)
    if (startIterations > 0) # but, make the rest of them uniform - get info later
        eHapsCurrent[(2 * windowSNPs + 1):nSNPs, ] <- 0.25
        ## initialize using rate assuming 0.5 cM/Mb, nSNPs = 100
    alphaMatCurrent <- matrix(1 / K, nrow = (nGrids - 1), ncol = K)
    priorCurrent <- rep(1 / K, K)
    hapSumCurrent <- array(1 / K, c(nGrids, K))
    reference_panel_SNPs <- array(FALSE, nSNPs)
    ##
    if (is.null(grid_distances)) {
        dl <- diff(L)
    } else {
        dl <- grid_distances
    }
    sigmaCurrent <- exp(-nGen * expRate / 100 / 1000000 * dl)

    
    if (reference_haplotype_file == "") {
        ##
        to_out <- list(
            sigmaCurrent = sigmaCurrent,
            eHapsCurrent = eHapsCurrent,
            alphaMatCurrent = alphaMatCurrent,
            hapSumCurrent = hapSumCurrent,
            priorCurrent = priorCurrent,
            reference_panel_SNPs = reference_panel_SNPs
        )
        ##
    } else {
        ## 
        to_out <- get_and_initialize_from_reference(
            sigmaCurrent = sigmaCurrent, eHapsCurrent = eHapsCurrent, alphaMatCurrent = alphaMatCurrent, hapSumCurrent = hapSumCurrent, priorCurrent = priorCurrent, ## first bit is default
            reference_haplotype_file = reference_haplotype_file, reference_legend_file = reference_legend_file, reference_sample_file = reference_sample_file, reference_populations = reference_populations, reference_phred = reference_phred, reference_iterations = reference_iterations, nSNPs = nSNPs, K = K, L = L, pos = pos, inputBundleBlockSize = inputBundleBlockSize, nCores = nCores, regionName = regionName, alleleCount = alleleCount, startIterations = startIterations, windowSNPs = windowSNPs, expRate = expRate, nGen = nGen, tempdir = tempdir, outputdir = outputdir, pseudoHaploidModel = pseudoHaploidModel, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, minRate = minRate, maxRate = maxRate, regionStart = regionStart, regionEnd = regionEnd, buffer = buffer, niterations = niterations, grid = grid, grid_distances = grid_distances, nGrids = nGrids, reference_shuffleHaplotypeIterations = reference_shuffleHaplotypeIterations, L_grid = L_grid, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, shuffle_bin_radius = shuffle_bin_radius, snps_in_grid_1_based = snps_in_grid_1_based)
    }

    print_message("Done parameter initialization")

    return(to_out)
}














## ugh this needs a clean
downsample_the_samples <- function(
    downsampleSamplesKeepList,
    downsampleSamples,
    highCovInLow,
    nCores,
    blockSize,
    N,
    inputBundleBlockSize,
    bundling_info,
    sampleNames,
    regionName,
    tempdir,
    outputdir,
    environment
) {
    
    print_message(paste0("Begin downsample samples to ",100 * downsampleSamples , "% of the original numbers"))
    keep <- array(FALSE,N)
    keelList <- NULL
    ## keep - either those on a list, or the high coverage
    if (is.na(downsampleSamplesKeepList) == FALSE) {
        if (file.exists(downsampleSamplesKeepList)) {
            keepList <- read.table(downsampleSamplesKeepList)
            keep[unlist(keepList)]=TRUE
        }
    }
    if (length(highCovInLow) > 0) {
        keep[highCovInLow] <- TRUE
    }
    ## sample remaining samples to keep
    which <- sample(
        1:sum(keep == FALSE),
        N * downsampleSamples - sum(keep)
    )
    keep[keep == FALSE][which] <- TRUE

    new_N <- sum(keep)
    keepL <- which(keep)
    new_sampleNames <- sampleNames[keep]
    
    if (length(highCovInLow > 0)) {
        new_highCovInLow <- match(sampleNames[highCovInLow], new_sampleNames)
    }
    
    
    new_bundling_info <- get_bundling_position_information(
        N = new_N,
        nCores = nCores,
        blockSize = inputBundleBlockSize
    )
    
    sampleRanges <- getSampleRange(new_N, nCores)
    
    f <- function(sampleRange, keepL, tempdir, regionName, bundling_info, new_bundling_info) {
        bundledSampleReads <- NULL
        for(iSample in sampleRange[1]:sampleRange[2]) {
            
            oldSample <- keepL[iSample]
            newSample <- iSample
            
            out <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = oldSample,
                bundling_info = bundling_info,
                bundledSampleReads = bundledSampleReads
            )
            sampleReads <- out$sampleReads
            bundledSampleReads <- out$bundledSampleReads

            ## hack - if inputBundleBlockSize = NA, then would overwrite
            ## so save to slightly different name, then push
            regionNameLocal <- regionName
            if (length(bundling_info) == 0)
                regionNameLocal <- paste0(regionName, ".TEMP.")
            
            save(
                sampleReads,
                file = file_sampleReads(tempdir, newSample, regionNameLocal),
                compress = FALSE
            )
            
            if (length(new_bundling_info) > 0) {
                last_in_bundle <- new_bundling_info$matrix[newSample, "last_in_bundle"]
                if (last_in_bundle == 1) {
                    bundle_inputs_after_generation(
                        bundling_info = new_bundling_info,
                        iBam = newSample,
                        dir = tempdir,
                        regionName = regionName
                    )
                }
            }
            
        }
    }
    
    
    if(environment=="server") {
        out2=mclapply(sampleRanges, mc.cores=nCores,FUN=f, keepL= keepL, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info, new_bundling_info = new_bundling_info)
    }
    if(environment=="cluster") {
        cl = makeCluster(nCores, type = "FORK")
        out2 = parLapply(cl, sampleRanges, fun=f, keepL= keepL, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info, new_bundling_info = new_bundling_info)
        stopCluster(cl)
    }
    check_mclapply_OK(out2, "There has been an error downsampling the samples. Please see error message above")

    ## continuation of hack from above
    ## do not want to over-write files
    if (length(bundling_info) == 0) {
        ## remove original ones, push new ones
        for(iSample in 1:N) {
            unlink(file_sampleReads(tempdir, iSample, regionName))
        }
        regionNameLocal <- paste0(regionName, ".TEMP.")
        for(iSample in 1:new_N) {
            system(
                paste0(
                    "mv ",
                    "'", file_sampleReads(tempdir, iSample, regionNameLocal), "'",
                    " ",
                    "'", file_sampleReads(tempdir, iSample, regionName), "'"
                )
            )
        }
    }
    

    save(
        sampleNames,
        keepL,
        file = file.path(
            outputdir, "RData", paste0("keepL.", regionName, ".RData")
        )
    )

  print_message("End downsample samples")
    
    return(
        list(
            N = new_N,
            sampleNames = new_sampleNames,
            bundling_info = new_bundling_info,
            highCovInLow = new_highCovInLow
        )
    )
}




## this generates the input data
## three options
## regenerateInput == TRUE, is.na(inputBundleBlockSize) == TRUE
##    make per-sample BAMs
## regenerateInput == TRUE, is.na(inputBundleBlockSize) == FALSE
##    make bundled inuput
## regenerateIput == FALSE, is.na(inputBundleBlockSize) == FALSE
##    check or rebundle input

#' @export
generate_or_refactor_input <- function(
    regenerateInput,
    bundling_info,
    L,
    pos,
    nSNPs,
    bam_files,
    cram_files,
    reference,
    iSizeUpperLimit,
    bqFilter,
    chr,
    outputdir,
    N,
    downsampleToCov,
    sampleNames,
    inputdir,
    useSoftClippedBases,
    regionName,
    tempdir,
    chrStart,
    chrEnd,
    generateInputOnly,
    environment,
    nCores,
    save_sampleReadsInfo
) {
    if (regenerateInput == TRUE) {
        ## handles both bundled and unbundled inputs
        out <- generate_input(bundling_info = bundling_info, L = L, pos = pos, nSNPs = nSNPs, bam_files = bam_files, cram_files = cram_files, reference = reference, iSizeUpperLimit = iSizeUpperLimit, bqFilter = bqFilter, chr = chr, N = N, downsampleToCov = downsampleToCov, sampleNames = sampleNames, inputdir = inputdir, useSoftClippedBases = useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, environment = environment, nCores = nCores, save_sampleReadsInfo = save_sampleReadsInfo)
        if(generateInputOnly==TRUE) {
            save(
                pos,
                file = file.path(
                    outputdir, "RData",
                    paste0("pos.", regionName ,".RData")
                )
            )
            return(NULL)
        }
    } else { # regenerateInput = FALSE
        if (length(bundling_info) > 0) {
            rebundle_input(
                inputdir = inputdir,
                regionName = regionName,
                bundling_info = bundling_info,
                N = N,
                tempdir = tempdir,
                nCores = nCores,
                outputdir = outputdir
            )
        }
    }
}











## so were definitely generating input here
## decide whether it is on a server or a cluster
## that's it!
generate_input <- function(
  bundling_info,
  L,
  pos,
  nSNPs,
  bam_files,
  cram_files,
  reference,
  iSizeUpperLimit,
  bqFilter,
  chr,
  N,
  downsampleToCov,
  sampleNames,
  inputdir,
  useSoftClippedBases,
  regionName,
  tempdir,
  chrStart,
  chrEnd,
  environment,
  nCores,
  save_sampleReadsInfo
) {
    print_message("Generate inputs")
    chrLength <- NA
    if (is.na(chrStart)) {
        print_message("Get chromosome length")
        chrLength <- get_chromosome_length(iBam = 1, bam_files, cram_files, chr)
    }
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    if(environment == "server") {
        out <- mclapply(1:length(sampleRanges), mc.cores=nCores, FUN=loadBamAndConvert_across_a_range,sampleRanges = sampleRanges,bundling_info=bundling_info,L=L,pos=pos, nSNPs = nSNPs,bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
    }
    if(environment=="cluster") {
        cl = makeCluster(length(sampleRanges), type = "FORK")
        out = parLapply(cl, 1:length(sampleRanges), mc.cores=nCores, fun=loadBamAndConvert_across_a_range,sampleRange=sampleRange,bundling_info=bundling_info, L=L,pos=pos, nSNPs = nSNPs, bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
        stopCluster(cl)
    }
    if (length(unlist(out)) > 0)  {
        print_message(out[[which(sapply(out, length) > 0)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }
    print_message("Done generating inputs")
    return(NULL)
}


# loadBams over a range
# if appropriate, dump to disk
loadBamAndConvert_across_a_range <- function(
  iCore,
  sampleRanges,
  bundling_info,
  L,
  pos,
  nSNPs,
  bam_files,
  cram_files,
  reference,
  iSizeUpperLimit,
  bqFilter,
  chr,
  N,
  downsampleToCov,
  sampleNames,
  inputdir,
  useSoftClippedBases,
  regionName,
  tempdir,
  chrStart,
  chrEnd,
  chrLength,
  save_sampleReadsInfo
) {
    ##print_message(paste0("Convert inputs for iCore=", iCore))
    sampleRange <- sampleRanges[[iCore]]
    for(iBam in sampleRange[1]:sampleRange[2]) {
        loadBamAndConvert(iBam = iBam, L=L,pos=pos, nSNPs = nSNPs,bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
        ## if a bundle iteration, scoop up samples, dump to disk
        if (length(bundling_info) > 0) {
            last_in_bundle <- bundling_info$matrix[iBam, "last_in_bundle"]
            if (last_in_bundle == 1) {
                bundle_inputs_after_generation(
                    bundling_info = bundling_info,
                    iBam = iBam,
                    dir = inputdir,
                    regionName = regionName
                )
            }
        }
    }
    ##print_message(paste0("Done converting inputs for iCore=", iCore))
    return(NULL)
}


















get_RG_lines_from_SeqLib <- function(file) {
    header_string <- get_header_using_SeqLib(file)
    x <- strsplit(header_string, "\n")[[1]]
    return(x[substr(x, 1, 3) == "@RG"])
}

convert_sam_header_rg_tags_to_sample_name <- function(header, file) {
    rg_spots <- lapply(header, strsplit, split = "\t")
    if (length(rg_spots) == 0)
        stop(paste0("There is no RG tag with sample name in file:", file))
    sm <- sapply(rg_spots, function(x) {
        rg <- x[[1]]
        sm <- substr(rg[substr(rg, 1, 3) == "SM:"], 4, 1000)
        return(sm)
    })
    if (sum(sapply(sm, length)) == 0)
        stop(paste0("The RG tags do not contain SM entries for file:", file))
    if (length(unique(sm)) > 1)
        stop(paste0("There is more than one sample name in the header for:", file))
    return(sm[1])
}

get_sample_name_from_bam_file_using_external_samtools <- function(file) {
    header <- system(paste0("samtools view -H ", file, " | grep ^@RG"), intern = TRUE)
    return(convert_sam_header_rg_tags_to_sample_name(header = header, file = file))
}

get_sample_name_from_bam_file_using_SeqLib <- function(file) {
    header <- get_RG_lines_from_SeqLib(file)
    return(convert_sam_header_rg_tags_to_sample_name(header = header, file = file))
}



## for a set of bam or cram files
## use samtools to get the names of the samples
get_sample_names_from_bam_or_cram_files <- function(
    files,
    nCores,
    file_type = "BAM",
    verbose = TRUE
) {

    if (verbose)
        print_message(paste0("Get ", file_type, " sample names"))

    for(file in files)
        if (file.exists(file) == FALSE)
            stop(paste0("Cannot find ", file_type, " file:", file))

    sampleNames <- mclapply(
        files,
        mc.cores = nCores,
        get_sample_name_from_bam_file_using_SeqLib
    )

    sampleNames <- as.character(unlist(sampleNames))

    if (verbose)
        print_message(paste0(
            "Done getting ", file_type, " sample names"
        ))

    return(sampleNames)
}













## for a list of bam files or cram files
## for a bam, either return the bam name,
## or if a cram file, first decompress to bam format
## in the temporary directory
get_bam_name_and_maybe_convert_cram <- function(
    iBam,
    bam_files,
    cram_files,
    reference,
    tempdir,
    chr,
    chrStart,
    chrEnd
) {
    if (length(bam_files) > 0) {
        bamName <- bam_files[iBam]
        if (file.exists(bamName) == FALSE)
            stop(paste0("Cannot find file:", bamName))
    } else if (length(cram_files) > 0) {
        cramName <- cram_files[iBam]
        if (file.exists(cramName) == FALSE)
            stop(paste0("Cannot find file:", cramName))
        bamName <- cramName
    } else {
        stop("Both bam_files and cram_files are empty")
    }
    return(bamName)
}



get_index_for_bamName <- function(bamName) {
    idx1 <- paste(substr(bamName,1,nchar(bamName)-3),"bai",sep="")
    idx2 <- paste(bamName,".bai",sep="")
    if(file.exists(idx1)) {
        return(idx1)
    } else if (file.exists(idx2)) {
        return(idx2)
    } else {
        stop(
            paste0(
                "there is no index for bam file:",
                bamName
            )
        )
    }
}


deal_with_soft_clipped_bases <- function(splitCigarRead, useSoftClippedBases, posRead, seqRead, qualRead) {
    for(iRead in 1:length(seqRead)) {
        out <- cpp_deal_with_soft_clipped_bases(
            splitCigarRead = splitCigarRead[[iRead]],
            useSoftClippedBases = useSoftClippedBases,
            posRead = posRead[iRead],
            seqRead = seqRead[iRead],
            qualRead = qualRead[iRead]
        )
        splitCigarRead[[iRead]] <- out$splitCigarRead
        posRead[iRead] <- out$posRead
        seqRead[iRead] <- out$seqRead
        qualRead[iRead] <- out$qualRead
    }
    return(
        list(
            splitCigarRead = splitCigarRead ,
            posRead = posRead,
            seqRead = seqRead,
            qualRead = qualRead
        )
    )
}






## previous sampleReads comes straight from C++
## and is not linked by reads
## here, match by reads
merge_reads_from_sampleReadsRaw <- function(
    sampleReadsRaw,
    qname,
    strand,
    readStart,
    readEnd,
    iSizeUpperLimit,
    save_sampleReadsInfo = FALSE,
    qname_all,
    readStart_all,
    readEnd_all
) {
    ## wif: 1-based which read it came from
    wif <- sapply(sampleReadsRaw, function(x) x[[5]]) + 1
    qnameUnique <- unique(qname[wif])
    qnameInteger <- match(qname[wif], qnameUnique)
    ord <- order(qnameInteger) - 1 ## 0-based for C++
    qname_ord <- qname[ord + 1]
    qnameInteger_ord <- c(qnameInteger[ord + 1], - 2)
    readStart_ord <- readStart[ord + 1]
    readEnd_ord <- readEnd[ord + 1]
    ## qnameUnique is unique read names of used reads
    ## qnameInteger is integer version of that
    ## ord is 0-based ordered version of qnameInteger
    ## qnameInteger_ord is the ordered version of qnameInteger
    out <- cpp_read_reassign(
        ord = ord,
        qnameInteger_ord = qnameInteger_ord,
        qname = qname,
        strand = strand,
        sampleReadsRaw = sampleReadsRaw,
        readStart_ord = readStart_ord,
        readEnd_ord = readEnd_ord,
        readStart = readStart,
        readEnd = readEnd,
        iSizeUpperLimit = iSizeUpperLimit,
        save_sampleReadsInfo = save_sampleReadsInfo 
    )
    sampleReads <- out$sampleReads
    sampleReadsInfo <- out$sampleReadsInfo
    ##
    if (save_sampleReadsInfo) {
        ## 
        readMin <- array(-1, nrow(sampleReadsInfo))
        readMax <- array(-1, nrow(sampleReadsInfo))
        match_vec <- match(qname_all, sampleReadsInfo[, "qname"]) - 1 ## make 0-based
        match_vec[is.na(match_vec)] <- -1
        ##
        ## pass by reference
        get_min_from_position(
            match_vec = match_vec,
            readStart_all = readStart_all,
            readEnd_all = readEnd_all,
            readMin = readMin,
            readMax = readMax
        )
        ##
        sampleReadsInfo[, "minReadStart"] <- readMin
        sampleReadsInfo[, "maxReadEnd"] <- readMax
    }
    ## 
    ## which reads are saved (i.e. do not violate iSizeUpperLimit)    
    ##    save_read <-  head(out$save_read, length(qnameUnique))
    ## add information on "left-most" edge of read, "right-most" edge of read
    ## are read start, read end always ordered small to large?
    ##    sampleReadsInfo <- data.frame(
    ##        qname = qnameUnique[save_read],
    ##       strand = strand[match(qnameUnique, qname)][save_read] ## problem if varies
    ## )
    return(
        list(
            sampleReads = sampleReads,
            sampleReadsInfo = sampleReadsInfo
        )
    )
}


## if there are no reads, give the same some fake-ish reads
## anyway
get_fake_sampleReads <- function(name, mess) {
    print_message(
        paste0(
            "WARNING - sample ", name, " ",
            mess, ". ",
            "It is being given random reads. ",
            "Consider removing from analysis"
        )
    )
    return(make_fake_sampleReads())
}

make_fake_sampleReads <- function() {
    ## nU,d,phiU,pRU
    ## 1st - nU = 0-based number in read
    ## 2nd - d = "representative" position, 0-based index
    ## 3rd - phiU = for each SNP in read, probability of alt allele
    ## 4th - pRU = for each SNP in read, 0-based index in L
    sampleReads <- lapply(
        1:3,
        function(i)
        list(
            0,
            as.integer(i),
            matrix(as.integer(i),ncol=1,nrow=1),
            matrix(as.integer(i),ncol=1,nrow=1)
        )
    )
return(sampleReads)
}


## necessary to load chromosome in chunks
## assume the same length for all samples
get_chromosome_length <- function(iBam, bam_files, cram_files, chr) {
    if (length(bam_files) > 0) {
        file <- bam_files[iBam]
    } else if (length(cram_files) > 0) {
        file <- cram_files[iBam]
    }
    header <- system(paste0("samtools view -H ", file, " | grep ^@SQ"), intern = TRUE)
    seq_spots <- lapply(header, strsplit, split = "\t")
    if (length(seq_spots) == 0)
        stop(paste0("There is no @SQ tag (with sample name) for:", file))
    seq <- t(sapply(seq_spots, function(x) {
        y <- x[[1]]
        chrName <- substr(y[substr(y, 1, 3) == "SN:"], 4, 1000)
        chrLength <- substr(y[substr(y, 1, 3) == "LN:"], 4, 1000)
        return(c(chrName, chrLength))
    }))
    seq <- seq[seq[, 1] == chr, , drop = FALSE]
    if (nrow(seq) == 0)
        stop(paste0("Could not find chromosome length for file:", file))
    return(seq[1, 2])
}











## takes a lot of inputs
## converts a BAM file for one subject
## into sampleReads, saved into
## file_sampleReads(inputdir, iBam, regionName)
loadBamAndConvert <- function(
    iBam,
    L,
    pos,
    nSNPs,
    bam_files = NULL,
    cram_files = NULL,
    reference = "",
    iSizeUpperLimit = 600,
    bqFilter = 17,
    chr,
    N,
    downsampleToCov = 100,
    sampleNames,
    inputdir,
    useSoftClippedBases = FALSE,
    regionName,
    tempdir,
    chrStart,
    chrEnd,
    chrLength = NA,
    save_sampleReadsInfo = FALSE
) {

    sampleReadsInfo <- NULL ## unless otherwise created

    if ((iBam %% 100) == 0) {
        if (length(bam_files) > 0) {
            what <- "BAM"
        } else {
            what <- "CRAM"
        }
        print_message(
            paste0(
                "Load and convert ", what, " ",
                iBam, " of ", N
            )
        )
    }

    file_name <- get_bam_name_and_maybe_convert_cram(
        iBam,
        bam_files,
        cram_files,
        reference,
        tempdir,
        chr,
        chrStart,
        chrEnd
    )

    ref <- as.character(pos[, "REF"])
    alt <- as.character(pos[, "ALT"])
    out <- get_sampleReadsRaw_from_SeqLib(
        useSoftClippedBases = useSoftClippedBases,
        bqFilter = bqFilter,
        iSizeUpperLimit = iSizeUpperLimit,
        ref = ref,
        alt = alt,
        nSNPs = nSNPs,
        L = L,
        region = paste0(chr, ":", chrStart, "-", chrEnd),
        file_name = file_name,
        reference = reference,
        save_sampleReadsInfo = save_sampleReadsInfo
    )
    sampleReadsRaw <- out$sampleReadsRaw
    qname <- out$qname
    strand <- out$strand
    readStart <- out$readStart
    readEnd <- out$readEnd
    ## relevant if save_sampleReadsInfo = TRUE
    qname_all <- out$qname_all
    readStart_all <- out$readStart_all
    readEnd_all <- out$readEnd_all

    if (length(sampleReadsRaw) == 0) {
        sampleReads <- get_fake_sampleReads(
            name = sampleNames[iBam],
            mess = "has no informative reads"
        )
    } else {
        out <- merge_reads_from_sampleReadsRaw(
            sampleReadsRaw = sampleReadsRaw,
            qname = qname,
            strand = strand,
            readStart = readStart,
            readEnd = readEnd,
            iSizeUpperLimit = iSizeUpperLimit,
            save_sampleReadsInfo = save_sampleReadsInfo,
            qname_all = qname_all,
            readStart_all = readStart_all,
            readEnd_all = readEnd_all
        )
        sampleReads <- out$sampleReads
        sampleReadsInfo <- out$sampleReadsInfo
        out <- add_central_SNP_to_sampleReads(
            sampleReads = sampleReads,
            sampleReadsInfo = sampleReadsInfo,
            L = L
        )
        sampleReads <- out$sampleReads
        sampleReadsInfo <- out$sampleReadsInfo
    }

    out <- downsample(
        sampleReads = sampleReads,
        iBam = iBam,
        downsampleToCov = downsampleToCov,
        sampleNames = sampleNames,
        sampleReadsInfo = sampleReadsInfo
    )
    sampleReads <- out$sampleReads
    sampleReadsInfo <- out$sampleReadsInfo

    ## these are saved once - can be compressed!
    save(
        sampleReads,
        file = file_sampleReads(inputdir, iBam, regionName),
        compress = FALSE
    )

    if (save_sampleReadsInfo)
        save(
            sampleReadsInfo,
            file = file_sampleReadsInfo(inputdir, iBam, regionName),
            compress = FALSE
        )

    return(NULL)
}


##
## take a single, representative SNP from each - use this instead for position, etc.
##
add_central_SNP_to_sampleReads <- function(
    sampleReads,
    sampleReadsInfo,
    L
){
    sampleReads <- lapply(
        sampleReads,
        function(x) {
        if(length(x[[4]]) > 2) {
            ## take middle spot
            y <- L[x[[4]] + 1]
            x[[2]] <- as.integer(
                x[[4]][which.min(abs(y - mean(y)))]
            )
        } else if (length(x[[4]]) == 2) {
            ## take at random
            x[[2]] <- as.integer(sample(x[[4]], 1))
        } else {
            ## otherwise, make them the same
            x[[2]] <- as.integer(x[[4]])
        }
        return(x)
    }
    )
    t1 <- unlist(lapply(sampleReads,function(x) x[[2]]))
    sampleReads <- sampleReads[order(t1)]
    sampleReadsInfo <- sampleReadsInfo[order(t1), ]
    return(
        list(
            sampleReads = sampleReads,
            sampleReadsInfo = sampleReadsInfo
        )
    )
}





#' @export
shrinkReads <- function(
    N,
    nCores,
    originalRegionName,
    regionName,
    bundling_info,
    tempdir,
    inputdir,
    inRegionL,
    environment,
    regenerateInput,
    inputBundleBlockSize
) {

    if((regenerateInput == TRUE) | (regionName == originalRegionName)) {

        print_message("Copying files onto tempdir")
        if (is.na(inputBundleBlockSize)) {
            file_with_files_to_transfer <- file.path(tempdir, "files_to_transfer.txt")
            command1 <- paste0(
                'cd ', shQuote(inputdir), ' && find . -name "',
                'sample.*.input.', regionName, '.RData',
                '" > ', shQuote(file_with_files_to_transfer)
            )
            system(command1)
        } else {
            file_with_files_to_transfer <- file.path(tempdir, "files_to_transfer.txt")
            command1 <- paste0(
                'cd ', shQuote(inputdir), ' && find . -name "',
                'bundledSamples.*-*.', regionName, '.RData',
                '" > ', shQuote(file_with_files_to_transfer)
            )
            system(command1)
        }
        command2 <- paste0(
            "rsync -a --files-from=",
            shQuote(file_with_files_to_transfer),  " ",
            shQuote(inputdir), " ",
            shQuote(tempdir)
        )
        system(command2)
        print_message("Done copying files onto tempdir")
    } else {

        print_message("Begin shrink reads")
        sampleRanges <- getSampleRange(N = N, nCores = nCores)
        if(environment == "server") {
            out <- mclapply(sampleRanges,FUN=shrinkReads_on_range, mc.cores=nCores, originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL)
        }
        if(environment == "cluster") {
            cl <- makeCluster(nCores, type = "FORK")
            out <- parLapply(cl, sampleRanges,fun=shrinkReads_on_range,originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL)
            stopCluster(cl)
        }
        error_check <- sapply(out, class) == "try-error"
        if (sum(error_check) > 0) {
            print_message(out[[which(error_check)[1]]])
            stop("There has been an error generating the input. Please see error message above")
        }
        print_message("End shrink readsk")

    }

    return(NULL)
}





shrinkReads_on_range <- function(
  sampleRange,
  originalRegionName,
  regionName,
  bundling_info,
  tempdir,
  inputdir,
  inRegionL
) {
    ## here, we have a range to operate on
    ## loop over samples. get read. if no bundling, replace into tempdir
    ## otherwise, write every once in a while
    bundledSampleReads <- NULL
    for(iSample in sampleRange[1]:sampleRange[2]) {
        out <- get_sampleReads_from_dir_for_sample(
            dir = inputdir,
            regionName = originalRegionName,
            iSample = iSample,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads

        sampleReads <- shrinkReads_on_sample_reads(sampleReads = sampleReads, inRegionL = inRegionL)
        save(
            sampleReads,
            file = file_sampleReads(tempdir, iSample, regionName),
            compress = FALSE
        )
        if (length(bundling_info) > 0) {
            last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
            if (last_in_bundle == 1) {
                bundle_inputs_after_generation(
                    bundling_info = bundling_info,
                    iBam = iSample,
                    dir = tempdir,
                    regionName = regionName
                )
            }
        }
    }
    return(NULL)
}





# pass it a sampleReads
# shrink it
shrinkReads_on_sample_reads <- function(
  sampleReads,
  inRegionL
) {
    ##
    ## figure out which intersect - keep
    ##
    ## snps intersected
    x=lapply(sampleReads,function(x) x[[4]]) # slow - next 3 fast
    snpPerRead=unlist(lapply(x,length))
    wrif=rep(1:length(sampleReads),snpPerRead)
    snpi=unlist(x)+1 # SNPs intersected - 1 based
    ## keep those with SNPs completely within range
    w=  snpi>=inRegionL[1] & snpi<=inRegionL[2] # inRegionL is 1based
    readsToKeep=unique(wrif[w])
    sampleReads=sampleReads[readsToKeep]
    ##
    ## now - remove those with SNPs outside
    ##
    x=lapply(sampleReads,function(x) x[[4]])
    snpPerRead=unlist(lapply(x,length))
    wrif=rep(1:length(sampleReads),snpPerRead)
    snpi=unlist(x)+1 # SNPs intersected
    ## keep those with SNPs completely within range
    w=  snpi<inRegionL[1] | snpi>inRegionL[2]
    readsToRemove=unique(wrif[w])
    if(length(readsToRemove) > 0)
        sampleReads <- sampleReads[ - readsToRemove]
    ##
    ## now - change start and end
    ##
    ##a=0 # 0-based - what we keep
    ##b=diff(inRegionL)-1
    sampleReads <- lapply(
        sampleReads,
        function(x, inRegionL) {
        x[[2]] <- x[[2]] - (inRegionL[1] - 1)
        x[[4]] <- x[[4]] - (inRegionL[1] - 1)
        return(x)
    },
    inRegionL = inRegionL)
    return(sampleReads)
}




getCentral=function(LL)
{
    if(length(LL)==1) return(1)
    if(length(LL)==2) return(sample(2,1))
    if(length(LL)>=3) return(which.min(sapply(1:length(LL),function(i) sum(abs(LL[i]-LL[-i])))))
}




#subsetSNPsfile="/well/myers/rwdavies/converge/imputation.chr20.list.txt"
subsetSNPsFunction=function(N,subsetSNPsfile,regionName,tempdir,L,environment,nCores,outputdir)
{
  # load positions
  subsetList=scan(subsetSNPsfile)
  keep=is.na(match(L,subsetList))==FALSE
  # for those which are kept, do re-mapping
  keepL=array(NA,length(keep))
  keepL[keep]=(1:sum(keep))-1
  # loop through files
  f=function(iSample,tempdir,regionName,L,keep,outputdir)
  {
    load(file = file_sampleReads(tempdir, iSample, regionName))
    # for loop through on sampleReads? fast enough?
    for(i in 1:length(sampleReads))
    {
      c=sampleReads[[i]][[3]]
      d=sampleReads[[i]][[4]] # 0-based
      e=keep[d+1]
      if(sum(e==FALSE)==0) # keep everything!
      {
        # re-map only
        sampleReads[[i]][[2]]=keepL[sampleReads[[i]][[2]]+1]
        sampleReads[[i]][[4]]=keepL[sampleReads[[i]][[4]]+1]
      } else {
        if(sum(e==TRUE)==0) # otherwise, remove!
        {
          sampleReads[[i]]=NA
        } else {
          # keepL is 0-based
          d2=keepL[matrix(d[e],ncol=1)+1] # 0-based
          c2=matrix(c[e],ncol=1)
          # get central position
          sampleReads[[i]]=list(sum(e)-1,d2[getCentral(L[d2+1])],c2,d2)
        }
      }
    }
    # remove those which are now NA
    sampleReads=sampleReads[sapply(sampleReads,length)==4]
    # done!
    save(sampleReads,file=paste(tempdir,"sample.",iSample,".input.",regionName,".RData",sep=""))
  }
  if(environment=="server")
  {
    out2=mclapply(1:N,FUN=f,tempdir=tempdir,regionName=regionName,L=L,keep=keep,mc.cores=nCores,outputdir=outputdir)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, as.list(1:N), fun=f,tempdir=tempdir,regionName=regionName,L=L,keep=keep,outputdir=outputdir)
    save(out2,file=paste(outputdir,"debug/out2.RData",sep=""))
    stopCluster(cl)
  }
  #
  # done!
  #
  return(list(keep=keep))
}



get_depth_per_SNP_for_sampleReads <- function(sampleReads = NULL, pb_pos = NULL) {
    if (is.null(pb_pos))
        pb_pos <- unlist(lapply(sampleReads, function(x) x[[4]])) ## 0-based
    ## get depth for every sites
    depth_per_SNP <- increment2N(
        y = as.numeric(rep(1, length(pb_pos))),
        z = as.numeric(pb_pos),
        yT = as.integer(length(pb_pos)),
        xT = as.integer(max(pb_pos))
    )
    return(depth_per_SNP)
}


## remove reads as necessary
## to ensure that for that sample, there is no greater than downsampleToCov coverage
## this is not done in the most intelligent way, but what should be a fast way for R
## removing SNPs with >50X coverage shouldn't affect a low coverage imputation method much
downsample <- function(
    sampleReads,
    iBam,
    downsampleToCov,
    sampleNames,
    sampleReadsInfo,
    verbose = TRUE
) {
    ## get every position found in the reads
    pb_pos <- unlist(lapply(sampleReads, function(x) x[[4]])) ## 0-based
    ## get depth at every site
    depth_per_SNP <- get_depth_per_SNP_for_sampleReads(pb_pos = pb_pos)
    ## determine which are bad SNPs that require downsampling
    offending_SNPs <- which(depth_per_SNP > downsampleToCov) ## 1-based
    if (length(offending_SNPs) > 0) {
        ## for offending SNPs, figure out what reads intersect them
        a <- sapply(sampleReads, function(x) return(length(x[[3]])))
        pb_read <- unlist(lapply(1:length(a), function(x) rep(x, a[x]))) - 1 ## 0-based
        ## figure out what are the offending bases, and hence reads, that need to be considered
        which_bases <- is.na(match(pb_pos, offending_SNPs - 1)) == FALSE
        reads_per_SNP <- tapply(pb_read[which_bases], pb_pos[which_bases], I) ## hopefully not slow
        toRemove <- array(FALSE, length(sampleReads))
        ## now - loop over offending SNPs, if still offending, remove reads
        for(bad_snp in offending_SNPs) {
            ## check again
            if (depth_per_SNP[bad_snp] > downsampleToCov) {
                reads_at_SNP <- reads_per_SNP[[as.character(bad_snp - 1)]] ## 0-based
                ## remove from consideration if they've been removed already
                reads_at_SNP <- reads_at_SNP[toRemove[reads_at_SNP + 1] == FALSE]
                ## remove those than span SNPs more than once
                t <- table(reads_at_SNP)
                reads_to_remove_multi_span <- as.integer(names(t)[t > 1])
                reads_single_span <- as.integer(names(t)[t == 1])
                cov_due_to_multi_span_reads <- sum(is.na(match(reads_at_SNP, reads_to_remove_multi_span)) == FALSE)
                ## now, remove reads at random that intersect this SNP until desired coverage reached
                extra_coverage <- depth_per_SNP[bad_snp] - downsampleToCov - cov_due_to_multi_span_reads
                if (extra_coverage > 0) {
                    reads_to_remove_single_span <- sample(
                        reads_single_span, extra_coverage, replace = FALSE
                    )
                } else {
                    reads_to_remove_single_span <- NULL
                }
                reads_to_remove <- c(reads_to_remove_single_span, reads_to_remove_multi_span)
                toRemove[reads_to_remove + 1] <- TRUE
                ## modify depth_per_SNP appropriately from removing all those reads
                for(r in (reads_to_remove + 1)) {
                    w <- sampleReads[[r]][[4]] + 1 # 1-based
                    for(w2 in w)
                        depth_per_SNP[w2] <- depth_per_SNP[w2] - 1
                }
            }
        }
        ## done!
        if (verbose)
            print_message(paste0(
                "downsample sample ", sampleNames[iBam], " - ",
                sum(toRemove), " of ", length(sampleReads),
                " reads removed "
            ))
        sampleReads <- sampleReads[toRemove == FALSE]
        sampleReadsInfo <- sampleReadsInfo[toRemove == FALSE, ]
    }
    return(
        list(
            sampleReads = sampleReads,
            sampleReadsInfo = sampleReadsInfo
        )
    )
}




#' @export
convertScaledBQtoProbs <- function(x) {
    ## for a matrix with 1 column, with base qualities, get the probabilities of the ref or alternate
    ## recall < 0 is ref, >1 is alternate
    ## output is a matrix, two columns
    n <- nrow(x)
    o <- array(0, c(n, 2))
    e <- 10 ** ( - abs(x)/10)
    ## fill in
    w <- cbind(1:n, 2 - as.integer(x < 0))
    o[w] <- 1 - e
    w[, 2] <- 3 - w[, 2]
    o[w] <- e * 1 / 3
    return(o)
}


buildAlleleCount_subfunction <- function(
    i                                        ,
    nSNPs,
    tempdir,
    regionName,
    bundling_info,
    allSampleReads
) {
    alleleCount <- array(0,c(nSNPs, 3))
    bundledSampleReads <- NULL
    for(iSample in i[1]:i[2]) {
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
        ## get positions and intensities
        a=unlist(sapply(sampleReads,function(x) x[[3]]))
        b=unlist(sapply(sampleReads,function(x) x[[4]]))
        bqProbs=convertScaledBQtoProbs(matrix(a,ncol=1))
        ## y is numeric, z = integer counts (can have 0),
        c1 <- increment2N(
            y = as.numeric(bqProbs[,1]),
            z = as.numeric(b),
            yT = as.integer(nrow(bqProbs)),
            xT = as.integer(nSNPs - 1)
        )
        c2 <- increment2N(
            y = as.numeric(bqProbs[,2]),
            z = as.numeric(b),
            yT = as.integer(nrow(bqProbs)),
            xT = as.integer(nSNPs - 1)
        )
        alleleCount[,1]=alleleCount[,1]+c2 # fixed nov 6 2015 - was backward
        alleleCount[,2]=alleleCount[,2]+c1+c2
    }
    return(
        list(
            alleleCount = alleleCount
        )
    )
}



#' @export
buildAlleleCount <- function(
    nSNPs,
    N,
    nCores,
    regionName,
    tempdir,
    bundling_info,
    allSampleReads
) {

    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    print_message("Generate allele count")

    out2 <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = buildAlleleCount_subfunction,
        nSNPs = nSNPs,
        tempdir = tempdir,
        regionName = regionName,
        bundling_info = bundling_info,
        allSampleReads = allSampleReads
    )

    check_mclapply_OK(out2)

    alleleCount <- array(0, c(nSNPs, 3))
    for(i in 1:length(sampleRanges)) {
        alleleCount <- alleleCount + out2[[i]][["alleleCount"]]
    }
    alleleCount[, 3] <- alleleCount[, 1] / alleleCount[, 2]

    alleleCount[is.na(alleleCount[,3]),3] <- 0 # not sure why these would exist anymore
    print_allele_count(alleleCount, N)
    print_message("Done generating allele count")

    return(alleleCount)
}



## basically print quantile using message
print_allele_count <- function(alleleCount, N) {
    print_message("Quantiles across SNPs of per-sample depth of coverage")
    ## or, print it my way
    prob <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    x <- quantile(alleleCount[,2] / N, prob = prob)
    fancy_x <- format(round(x, 3), nsmall = 3)
    top <- ""
    bottom <- ""
    for (i in 1:length(x)) {
        c1 <- nchar(names(x)[i])
        c2 <- nchar(fancy_x[i])
        ## pad with difference + 1 for spacing
        c_width <- max(c1, c2) + 1
        top <- paste0(top, paste0(rep(" ", c_width - c1), collapse = ""), names(x)[i])
        bottom <- paste0(bottom, paste0(rep(" ", c_width - c2), collapse = ""), fancy_x[i])
    }
    ##
    print_message(top)
    print_message(bottom)
}


downsampleToFraction_a_range <- function(
    sampleRange,
    tempdir,
    regionName,
    downsampleFraction,
    bundling_info
) {
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
        bundledSampleReads <- out$bundledSampleReads
        
        ## subset
        keep <- runif(length(sampleReads))<downsampleFraction
        ## keep 1 no matter what
        if(sum(keep)==0) {
            keep[sample(length(keep),1)]=TRUE
        }
        sampleReads=sampleReads[keep]
        ## save result back to disk
        save(sampleReads,file=file_sampleReads(tempdir, iSample, regionName), compress = FALSE)

        ## bundle back together
        if (length(bundling_info) > 0) {
            last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
            if (last_in_bundle == 1) {
                bundle_inputs_after_generation(
                    bundling_info = bundling_info,
                    iBam = iSample,
                    dir = tempdir,
                    regionName = regionName
                )
            }
        }
    }
    return(NULL)
}


## build allele count matrix from input RData files

#' @export
downsampleToFraction <- function(
    N,
    nCores,
    downsampleFraction,
    regionName,
    tempdir,
    environment,
    bundling_info
) {
    print_message("Begin downsampling reads")
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    out2 <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        FUN = downsampleToFraction_a_range,
        tempdir = tempdir,
        regionName = regionName,
        downsampleFraction = downsampleFraction,
        bundling_info = bundling_info
    )
    check_mclapply_OK(out2, stop_message = "There has been an error downsampling the reads. Please see error message above")
    print_message("Done downsampling reads")
    return(NULL)
}




## for a sampleRange, like two numbers c(10, 20)
## and some outputBlockSize, like 5
## get minus 1 based range
## so gives back 9, 15, 20
## telling you to range from 9+1 -> 15, then 15+1 -> 20
getOutputBlockRange <- function(
    sampleRange,
    outputBlockSize
) {
    if (sampleRange[1] == sampleRange[2])
        return(c(sampleRange[1] - 1, sampleRange[2]))
    d <- sampleRange[2] - sampleRange[1]
    if(outputBlockSize > d)
        outputBlockSize = d
    outputBlockRange <- c(
        sampleRange[1] - 1 + which((sampleRange[1]:sampleRange[2] -sampleRange[1])%%outputBlockSize==0),
        sampleRange[length(sampleRange)]
    )
    ## subtract one from first entry
    outputBlockRange[1] <- outputBlockRange[1] - 1
    if (outputBlockRange[length(outputBlockRange)] == outputBlockRange[length(outputBlockRange) - 1])
        outputBlockRange <- outputBlockRange[-length(outputBlockRange)]
    return(outputBlockRange)
}



## basically, replicate how I mclapply
## used to be manual, not basically uses cut and reformats
getSampleRange <- function(
  N,
  nCores
) {
    ## upon closer inspection, this might be slow? oh well
    if (nCores == 1) {
        x <- rep(1, N)
    } else {
        x <- as.integer(cut(1:N, nCores))
    }
    w <- which(diff(x) > 0)
    start <- c(1, w + 1)
    end <- c(w, N)
    sampleRange <- lapply(1:nCores, function(i_core) {
        x <- c(start[i_core], end[i_core])
        if (sum(is.na(x)) > 0)
            return(NULL)
        return(x)
    })
    sampleRange <- sampleRange[sapply(sampleRange, is.null) == FALSE]
    return(sampleRange)
}



## for N, nCores, blockSize
## return a matrix and a list
## -- for the matrix,
## with N rows, one per sample, in sampleNames order
## where first column is what core it is in
## where second column is what bundle in that core it is in
## where third column is position within that bundle
## -- for the list,
## bundlingList[[iCore]][[iBundle]] gives start and end

#' @export
get_bundling_position_information <- function(
    N,
    nCores,
    blockSize
) {
    if (is.na(blockSize)){
        return(NULL)
    }
    ## one entry per core, the range
    x3 <- getSampleRange(N = N, nCores = nCores)
    out <- lapply(1:nCores, function(i_core) {
        x2 <- getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize = blockSize
        )
        ## for each of these block ranges, get the bundle
        inputBundle <- unlist(lapply(1:(length(x2) - 1), function(i_sub_block) {
            return(rep(i_sub_block, x2[i_sub_block + 1] - x2[i_sub_block]))
        }))
        ## within each core and bundle, get the position
        positionWithinBundle <- unlist(lapply(1:(length(x2) - 1), function(i_sub_block) {
            return(1:(x2[i_sub_block + 1] - x2[i_sub_block]))
        }))
        ## this is 0-based - transform to 1 based and return
        ## also - return a list
        start <- x2[-length(x2)] + 1
        end <- x2[-1]
        bundleList <- lapply(1:length(start), function(i) {
            c(start[i], end[i])
        })
        return(
            list(
                inputCore = rep(i_core, diff(x3[[i_core]]) + 1),
                inputBundle = inputBundle,
                positionWithinBundle = positionWithinBundle,
                bundleList = bundleList
            ))
    })
    ## turn into a matrix
    iCore <- unlist(lapply(out, function(x) x$inputCore))
    iBundle <- unlist(lapply(out, function(x) x$inputBundle))
    iPosition <- unlist(lapply(out, function(x) x$positionWithinBundle))
    bundlePosition = cbind(
        iCore = iCore,
        iBundle = iBundle,
        iPosition = iPosition,
        last_in_bundle = 0
    )
    ## also - add "last_in_bundle
    w <- which(diff(bundlePosition[, "iPosition"]) < 0)
    w <- unique(c(w, N))
    bundlePosition[w, "last_in_bundle"] <- 1
    bundlingList <- lapply(out, function(x) x$bundleList)
    return(
        list(
            matrix = bundlePosition,
            list = bundlingList
        )
    )
}




# go above
# if pseudoHaploid, do this, before csi
initialize_readProbs <- function(
    N,
    nCores,
    pseudoHaploidModel,
    tempdir,
    bundling_info,
    regionName,
    initial_min_hapProb,
    initial_max_hapProb,
    allSampleReads
) {
    print_message("Initialize readProbs")
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    out <- mclapply(1:length(sampleRanges), mc.cores = nCores, function(iCore) {
        sampleRange <- sampleRanges[[iCore]]
        bundledSampleReads <- NULL
        for(iSample in sampleRange[1]:sampleRange[2]) {
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
            ## now - can make some hapProbs
            out <- get_default_hapProbs(
                pseudoHaploidModel,
                sampleReads,
                initial_min_hapProb,
                initial_max_hapProb
            )
            pRgivenH1 <- out$pRgivenH1
            pRgivenH2 <- out$pRgivenH2
            srp <- out$srp
            save(
                pRgivenH1, pRgivenH2, srp, 
                file = file_sampleProbs(tempdir, iSample, regionName)
            )
            if (length(bundling_info) > 0) {
                last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
                if (last_in_bundle == 1) {
                    bundle_inputs_after_generation(
                        bundling_info = bundling_info,
                        iBam = iSample,
                        dir = tempdir,
                        regionName = regionName,
                        what = "sampleProbs"
                    )
                }
            }
        }
        return(NULL)
    })
    error_check <- sapply(out, class) == "try-error"
    if (sum(error_check) > 0) {
        print_message(out[[which(error_check)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }
    print_message(paste0("Done initializing readProbs"))
    return(NULL)
}


get_default_hapProbs <- function(
  pseudoHaploidModel,
  sampleReads,
  initial_min_hapProb = 0.4,
  initial_max_hapProb = 0.6
) {
    a1 <- initial_min_hapProb
    b1 <- initial_max_hapProb
    pRgivenH1=a1+(b1-a1)*runif(length(sampleReads))
    pRgivenH2=a1+(b1-a1)*runif(length(sampleReads))
    srp <- unlist(lapply(sampleReads,function(x) x[[2]]))
    return(list(
        pRgivenH1 = pRgivenH1,
        pRgivenH2 = pRgivenH2,
        srp = srp
    ))
}





run_forward_backwards <- function(
    tempdir = tempdir,
    regionName = "region",
    iSample = 1,
    sampleReads,
    method,
    K,
    priorCurrent,
    alphaMatCurrent_t,
    eHapsCurrent_t,
    transMatRate_t_H,
    transMatRate_t_D,    
    iteration = NA,
    niterations = NA,
    K_random = NA,
    grid_eHaps_distance = NA,
    maxDifferenceBetweenReads = 1000,
    maxEmissionMatrixDifference = 1e10,
    highCovInLow = NULL,
    Jmax = 10,
    nor = NULL,
    pRgivenH1 = NULL,
    pRgivenH2 = NULL,
    srp = NULL,
    pseudoHaploidModel = 9,
    run_fb_subset = FALSE,
    alphaBetaBlock = NULL,
    run_fb_grid_offset = 0,
    suppressOutput = as.integer(1),
    generate_fb_snp_offsets = FALSE,    
    blocks_for_output = array(0, c(1, 1)),
    i_snp_block_for_alpha_beta = 1,
    return_a_sampled_path = FALSE,
    snp_start_1_based = NA,
    snp_end_1_based = NA,
    grid = -1,
    return_genProbs = FALSE,
    return_extra = FALSE,
    update_in_place = FALSE,
    gammaUpdate_t = array(0, c(1, 1, 1)),
    jUpdate_t = array(0, c(1, 1)),
    hapSum_t = array(0, c(1, 1)),
    priorSum = array(0, 1),
    pass_in_alphaBeta = FALSE,        
    alphaHat_t = array(0, c(1, 1)),
    betaHat_t = array(0, c(1, 1)),
    output_haplotype_dosages = FALSE
) {

    Jmax_local <- get_Jmax_wrt_iteration(Jmax, iteration, niterations)

    if (!return_genProbs) {
        if ((is.na(iteration) == FALSE) && (iteration == niterations)) {
            return_genProbs <- TRUE
        } else if (iSample %in% highCovInLow) {
            return_genProbs <- TRUE
        }
    }

    if (return_genProbs) {
        if (grid[1] == -1) {
            stop("Include grid with return_genProbs")
        }
    }

    if (is.na(snp_start_1_based)) {
        snp_start_1_based <- 1
        snp_end_1_based <- length(grid)
    }
    
    ## re-declare obvious stuff
    nSNPs <- ncol(eHapsCurrent_t)
    if (length(nor) == 0) {
        if (method == "pseudoHaploid") {
            nor <- 2
        } else {
            nor <- 1
        }
    }

    ## dummy up input
    if (is.null(alphaBetaBlock)) {
        alphaBetaBlock <- lapply(1:nor, function(nor) {
            return(
                list(
                    alphaHatBlocks_t = array(0, c(1, 1)),
                    betaHatBlocks_t = array(0, c(1, 1))
                )
            )
        })
    }

    fbsoL <- as.list(1:nor)

    if(method=="pseudoHaploid") {
        for (iNor in 1:nor) {        
            if (iNor==1) {
                pRgivenH1L <- pRgivenH1
                pRgivenH2L <- pRgivenH2
            }
            if (iNor==2) {
                pRgivenH1L <- pRgivenH2
                pRgivenH2L <- pRgivenH1
            }
            fbsoL[[iNor]] <- forwardBackwardHaploid(
                sampleReads = sampleReads,
                nReads = as.integer(length(sampleReads)),
                pi = priorCurrent,
                transMatRate_t_H = transMatRate_t_H,
                alphaMat_t = alphaMatCurrent_t,
                eHaps_t = eHapsCurrent_t,
                maxDifferenceBetweenReads = as.double(maxDifferenceBetweenReads),
                maxEmissionMatrixDifference = as.double(maxEmissionMatrixDifference),
                Jmax = Jmax_local,
                suppressOutput = suppressOutput,
                model = as.integer(pseudoHaploidModel),                
                pRgivenH1 = pRgivenH1L,
                pRgivenH2 = pRgivenH2L,
                run_pseudo_haploid = TRUE,
                alphaStart = alphaBetaBlock[[iNor]]$alphaHatBlocks_t[, i_snp_block_for_alpha_beta],
                betaEnd = alphaBetaBlock[[iNor]]$betaHatBlocks_t[, i_snp_block_for_alpha_beta],
                run_fb_subset = run_fb_subset,                
                run_fb_grid_offset = run_fb_grid_offset,
                blocks_for_output = blocks_for_output,
                generate_fb_snp_offsets = generate_fb_snp_offsets,
                return_extra = return_extra,
                update_in_place = update_in_place,
                gammaUpdate_t = gammaUpdate_t,
                jUpdate_t = jUpdate_t,
                hapSum_t = hapSum_t,
                priorSum = priorSum,
                pass_in_alphaBeta = pass_in_alphaBeta,
                alphaHat_t = alphaHat_t,
                betaHat_t = betaHat_t,
                snp_start_1_based = snp_start_1_based,
                snp_end_1_based = snp_end_1_based,
                grid = grid,
                output_haplotype_dosages = output_haplotype_dosages
            )
            fbsoL[[iNor]]$gammaK_t <- fbsoL[[iNor]]$gamma_t
        }
    }
    if(method=="diploid-inbred") {
        ## is diploid, but fully inbred, so use statistical haploid model
        iNor <- 1
        fbsoL[[iNor]] <- forwardBackwardHaploid(
            sampleReads = sampleReads,
            nReads = as.integer(length(sampleReads)),
            Jmax = Jmax_local,
            pi = priorCurrent,
            pRgivenH1 = as.double(0),
            pRgivenH2 = as.double(0),
            eHaps_t = eHapsCurrent_t,
            alphaMat_t = alphaMatCurrent_t,
            transMatRate_t_H= transMatRate_t_H,
            maxDifferenceBetweenReads = as.double(maxDifferenceBetweenReads),
            maxEmissionMatrixDifference = as.double(maxEmissionMatrixDifference),
            suppressOutput = suppressOutput,
            model = -1, ## irrelevant for haploid
            run_pseudo_haploid = FALSE,
            run_fb_subset = run_fb_subset,
            alphaStart = alphaBetaBlock[[iNor]]$alphaHatBlocks_t[, i_snp_block_for_alpha_beta],
            betaEnd = alphaBetaBlock[[iNor]]$betaHatBlocks_t[, i_snp_block_for_alpha_beta],
            run_fb_grid_offset = run_fb_grid_offset,
            blocks_for_output = blocks_for_output,
            generate_fb_snp_offsets = generate_fb_snp_offsets,
            return_extra = return_extra,
            update_in_place = update_in_place,
            gammaUpdate_t = gammaUpdate_t,
            jUpdate_t = jUpdate_t,
            hapSum_t = hapSum_t,
            priorSum = priorSum,
            pass_in_alphaBeta = pass_in_alphaBeta,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t,
            snp_start_1_based = snp_start_1_based,
            snp_end_1_based = snp_end_1_based,
            grid = grid,
            output_haplotype_dosages = output_haplotype_dosages
        )
        fbsoL[[iNor]]$gammaK_t <- fbsoL[[iNor]]$gamma_t            
    }
    
    if(method=="diploid") {

        fbsoL[[1]] <- forwardBackwardDiploid(
            sampleReads = sampleReads,
            nReads = as.integer(length(sampleReads)),
            pi = priorCurrent,
            eHaps_t = eHapsCurrent_t,
            alphaMat_t = alphaMatCurrent_t,
            transMatRate_t_D = transMatRate_t_D,
            Jmax = Jmax_local,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            suppressOutput = suppressOutput,
            run_fb_subset = run_fb_subset,
            alphaStart = alphaBetaBlock[[1]]$alphaHatBlocks_t[, i_snp_block_for_alpha_beta],
            betaEnd = alphaBetaBlock[[1]]$betaHatBlocks_t[, i_snp_block_for_alpha_beta],
            run_fb_grid_offset = run_fb_grid_offset,
            blocks_for_output = blocks_for_output,
            generate_fb_snp_offsets = generate_fb_snp_offsets,
            return_a_sampled_path = return_a_sampled_path,
            return_genProbs = return_genProbs,
            snp_start_1_based = snp_start_1_based,
            snp_end_1_based = snp_end_1_based,
            grid = grid,
            return_extra = return_extra,
            update_in_place = update_in_place,
            gammaUpdate_t = gammaUpdate_t,
            jUpdate_t = jUpdate_t,
            hapSum_t = hapSum_t,
            priorSum = priorSum,
            pass_in_alphaBeta = pass_in_alphaBeta,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t,
            output_haplotype_dosages = output_haplotype_dosages
        )
    }

    return(
        list(
            fbsoL = fbsoL
        )
    )

}




within_EM_per_sample_heuristics <- function(
    sampleReads,
    iSample,
    iiSample,
    K,
    fbsoL,
    method,
    regionName,
    L,
    L_grid,
    nGrids,
    grid,
    nor,
    K_subset,
    nbreaks,
    break_results,
    fromMat,
    srp,
    pRgivenH1,
    pRgivenH2,
    nSNPs,
    eHapsCurrent_t,
    alphaMatCurrent_t,
    sigmaCurrent,
    maxDifferenceBetweenReads,
    maxEmissionMatrixDifference,
    outputdir,
    iteration,
    restartIterations,
    highCovInLow,
    nsplit,
    best_K_for_sample,
    restartMatrixList, 
    readsTotal,
    readsSplit,
    pseudoHaploidModel,
    tempdir,
    samples_with_phase,
    i_core,
    sampleRange,
    fbd_store,
    plot_shuffle_haplotype_attempts,
    grid_distances,
    return_a_sampled_path,
    sampledPathList
) {
    ## 
    ##
    if (nbreaks > 0) {
        ## for each sample, get the changes between every 100
        for(iNor in 1:nor) {
            for(iBreak in 1:nbreaks) {
                from <- break_results[iBreak, "left_grid_break_0_based"]
                to <- break_results[iBreak, "right_grid_break_0_based"]
                hp1 <- fbsoL[[iNor]]$gammaK_t[, from + 1]
                hp2 <- fbsoL[[iNor]]$gammaK_t[, to + 1]
                fromMat[iBreak, , ] <-
                    fromMat[iBreak, , ] + hp1 %*% t(hp2)
            }
        }
        ## also, if icore is 1, then save up until 20 samples, or if < 20, plot once all done
        if (plot_shuffle_haplotype_attempts) {
            ## ignore pseudoHaploid effectively, just take first
            if (i_core == 1) {
                M <- min(20, sampleRange[2]) ## how many to plot
                if (iSample == 1) {
                    fbd_store <- as.list(1:M)
                }
                if (iSample <= M) {
                    fbd_store[[iSample]] <- fbsoL[[1]]
                }
                if (iSample == M) {
                    save(fbd_store, file = file_fbdStore(tempdir, regionName, iteration))
                    fbd_store <- NULL
                }
            }
        }
    }
    ##
    ## now - consider whether to re-start - ie fix gammas
    ##
    if(iteration %in% restartIterations) {
        out <- restartSample(
            sampleReads = sampleReads, srp = srp, pRgivenH1 = pRgivenH1, pRgivenH2 = pRgivenH2, fbsoL=fbsoL, nSNPs = nSNPs,eHapsCurrent_t=eHapsCurrent_t,sigmaCurrent=sigmaCurrent,alphaMatCurrent_t=alphaMatCurrent_t,Jmax=Jmax,priorCurrent=priorCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,maxEmissionMatrixDifference = maxEmissionMatrixDifference,outputdir=outputdir,iteration = iteration,restartIterations = restartIterations, method = method
        )
        restartMatrixList[[iSample]] <- out$m
        fbsoL <- out$fbsoL
    }
    ##
    ## update pseudo diploid probabilities
    ##
    if (method == "pseudoHaploid") {
        r <- pseudoHaploidUpdate(pRgivenH1=pRgivenH1,pRgivenH2=pRgivenH2,fbsoL=fbsoL,pseudoHaploidModel=pseudoHaploidModel,srp=srp, K = K)
        pRgivenH1 <- r$pRgivenH1
        pRgivenH2 <- r$pRgivenH2
        save(
            srp, pRgivenH1, pRgivenH2,
            file = file.path(tempdir, paste0("sample.",iSample,".readProbs.",regionName,".RData"))
        )
    }
    ##
    ## check high coverage samples if applicable
    ##
    if (is.na(match(iSample, highCovInLow)) ==FALSE) {
        gp_t <- calculate_gp_t_from_fbsoL(
            eHapsCurrent_t = eHapsCurrent_t,
            grid = grid,
            method = method,
            fbsoL = fbsoL
        )
        dosage <- gp_t[2, ] + 2 * gp_t[3, ]
        save(dosage, file=file_dosages(tempdir, iSample, regionName), compress = FALSE)
    } # end of check on high coverage sample
    ##
    ## check phasing samples if applicable
    ##
    if (method == "pseudoHaploid" &&
        is.na(match(iSample, samples_with_phase))==FALSE) {
        haplotypes <- array(0, c(nSNPs, 2))
        for(i_hap in 1:2)
            haplotypes[, i_hap] <- colSums(fbsoL[[i_hap]]$gamma_t * eHapsCurrent_t)
        save(haplotypes, file=file_haplotypes(tempdir, iSample, regionName))
    }
    ##
    ## do read splitting if the correct iteration
    ##
    if (nsplit > 0 & method != "pseudoHaploid") {
        ## BUGwerwer - this has never been on
        ## including never for pseudoHaploid
        ## fix this in the future
        if(method=="pseudoHaploid") {
            gammaK_t <- fbsoL[[1]]$gammaK_t + fbsoL[[2]]$gammaK_t
        } else {
            gammaK_t <- fbsoL[[1]]$gammaK_t
        }
        out <- findRecombinedReadsPerSample(
            gammaK_t = fbsoL[[1]]$gammaK_t,
            eHapsCurrent_t = eHapsCurrent_t,
            K = K,
            L = L,
            iSample = iSample,
            verbose = FALSE,
            sampleReads = sampleReads,
            tempdir = tempdir,
            regionName = regionName,
            grid = grid
        )
        readsTotal[iSample] <- out$readsTotal
        readsSplit[iSample] <- out$readsSplit
    }
    if (return_a_sampled_path) {
        ## add to list?
        sampledPathList[[iiSample]] <- lapply(fbsoL, function(x) x$sampled_path)
    }
    return(
        list(
            readsTotal = readsTotal,
            readsSplit = readsSplit,
            restartMatrixList = restartMatrixList,
            fromMat = fromMat,
            fbd_store = fbd_store,
            sampledPathList = sampledPathList
        )
    )
}




## so given variables below
## calculate output
calculate_gp_t_from_fbsoL <- function(
    eHapsCurrent_t,
    grid,
    method,
    fbsoL,
    snp_start_1_based = NA,
    snp_end_1_based = NA,
    grid_offset_0_based = 0
) {
    ## another check
    if (method != "diploid") {
        for(inor in 1:length(fbsoL)) {
            for(w in c(1, ncol(fbsoL[[inor]]$gamma_t))) {
                r <-  sum(fbsoL[[inor]]$gamma_t[, w])
                if ((r < 0.9) | (1.1 < r)) {
                    stop("Calculated gamma should sum to exactly 1 while there are observed entries below 0.9 or above 1.1. Please report this")
                }
            }
        }
    }
    if (is.na(snp_start_1_based)) {
        snp_start_1_based <- 1
        snp_end_1_based <- length(grid)
    }
    nSNPs <- snp_end_1_based - snp_start_1_based + 1
    K <- nrow(eHapsCurrent_t)
    ## if pseudo-haploid, get probabilities
    ## disable outputting for now
    if (method == "pseudoHaploid" && 1 == 0) {
        read_proportions <- estimate_read_proportions(
            sampleReads = sampleReads,
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,
            nSNPs = nSNPs
        )
    } else {
        read_proportions <- NULL
    }
    ## get allele counts for HWE, info calcs, AF
    ## get most likely genotype - add to count
    if (method == "diploid") {
        ## should have previously been created
        if (!("genProbs_t" %in% names(fbsoL[[1]]))) {
            stop("need genProbs_t in fbsoL")
        }
        gp_t <- fbsoL[[1]][["genProbs_t"]]
    } else if (method == "diploid-inbred") {
        ## this is a diploid organisms, so output genotypes appropriately
        ## there are only two posteriors here!
        gp_t <- array(0, c(3, nSNPs))        
        w2 <- snp_start_1_based:snp_end_1_based
        w <- grid[w2] + 1 - grid_offset_0_based
        gp_t[1, ] <- colSums(
            fbsoL[[1]]$gamma_t[, w, drop = FALSE] * (1 - eHapsCurrent_t[, w2, drop = FALSE])
        )
        gp_t[3, ] <- 1 - gp_t[1, ]        
    } else if (method == "pseudoHaploid") {
        gp_t <- array(0, c(3, nSNPs))
        w2 <- snp_start_1_based:snp_end_1_based
        w <- grid[w2] + 1 - grid_offset_0_based
        g10 <- colSums(fbsoL[[1]]$gamma_t[, w, drop = FALSE] * (1-eHapsCurrent_t[, w2, drop = FALSE]))
        g20 <- colSums(fbsoL[[2]]$gamma_t[, w, drop = FALSE] * (1-eHapsCurrent_t[, w2, drop = FALSE]))
        gp_t[1, ] <- g10 * g20
        gp_t[2, ] <- g10 * (1 - g20) + (1 - g10) * g20
        gp_t[3, ] <- (1 - g10) * (1 - g20)
    }
    r <- sum(gp_t[, 1])
    if ((r < 0.9) | (1.1 < r)) {
        stop("Calculated genotype probabilities contain entries below 0.9 or above 1.1 while these should be exactly 1. Please report this")
    }
    return(gp_t)
}





subset_of_complete_iteration <- function(sampleRange,tempdir,chr,K,K_subset, K_random, nSNPs, nGrids, priorCurrent,eHapsCurrent_t,alphaMatCurrent_t,sigmaCurrent,maxDifferenceBetweenReads,maxEmissionMatrixDifference, Jmax,highCovInLow,iteration,method,nsplit,expRate,minRate,maxRate,gen,outputdir,pseudoHaploidModel,outputHaplotypeProbabilities,switchModelIteration,regionName,restartIterations,refillIterations,outputBlockSize, bundling_info, transMatRate_t_H, transMatRate_t_D, sampleRanges, N, niterations, L, samples_with_phase, nbreaks, break_results, vcf.piece_unique, grid, B_bit_prob, start_and_end_minus_buffer, plot_shuffle_haplotype_attempts, blocks_for_output, allSampleReads, L_grid, grid_distances, minimizeSwitchingIterations, useTempdirWhileWriting) {

    ## know what core we are in for writing
    i_core <- match(sampleRange[1], sapply(sampleRanges, function(x) x[[1]]))
    who_to_run <- sampleRange[1]:sampleRange[2]    
    ## initialize bundling variables
    bundledSampleReads <- NULL
    bundledSampleProbs <- NULL
    
    ##
    ## initialize output matrices here
    ##
    priorSum <- array(0, K)
    alphaMatSum_t <- array(0, c(K, nGrids - 1))
    gammaSum_t <- array(0, c(K, nSNPs, 2))
    hapSum_t <- array(0,c(K, nGrids))

    ## decide if to re-generate
    pass_in_alphaBeta <- TRUE
    if (pass_in_alphaBeta) {
        if (method == "diploid") {
            alphaHat_t <- array(0, c(K * K, nGrids))
            betaHat_t <- array(0, c(K * K, nGrids))
        } else if ((method == "diploid-inbred") | (method == "pseudoHaploid")) {
            alphaHat_t <- array(0, c(K, nGrids))
            betaHat_t <- array(0, c(K, nGrids))
        } else {
            stop("bad method")
        }
    } else {
        alphaHat_t <- array(0, c(1, 1))
        betaHat_t <- array(0, c(1, 1))
    }

    ##
    ## other things sometimes needed
    ##
    if (nsplit == 1) {
        readsSplit <- array(0, N)
        readsTotal <- array(0, N)
    } else {
        readsSplit <- NULL
        readsTotal <- NULL
    }
    
    if (nbreaks > 0) {
        fromMat <- array(0, c(nbreaks, K, K))
    } else {
        fromMat <- NULL
    }
    if(iteration %in% restartIterations) {    
        restartMatrixList <- as.list(1:N)
    } else {
        restartMatrixList <- NULL
    }

    ##
    if (iteration == niterations) {
        allAlphaBetaBlocks <- as.list(1:length(who_to_run))
    } else {
        allAlphaBetaBlocks <- NULL
    }

    ## run either once (default) or more than once
    nor <- 1 # number of runs - 1 for normal diploid method, or haploid
    ## run both haplotypes
    if (method == "pseudoHaploid") {
        nor <- 2
    }
    best_K_for_sample <- 1:K ##
    pRgivenH1 <- NULL
    pRgivenH2 <- NULL
    srp <- NULL
    bundledSampleProbs <- NULL
    fbd_store <- NULL ## for plotting shuffle haplotype itetself

    ## on this iteration, sample paths
    if (iteration %in% minimizeSwitchingIterations) {
        return_a_sampled_path <- TRUE
        sampledPathList <- as.list(1:length(who_to_run))
    } else {
        return_a_sampled_path <- FALSE
        sampledPathList <- NULL
    }

    if (iteration == niterations) {
        generate_fb_snp_offsets <- TRUE
    } else {
        generate_fb_snp_offsets <- FALSE
    }

    for (iiSample in 1:length(who_to_run)) {
        
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

        out <- run_forward_backwards(
            tempdir = tempdir,
            regionName = regionName,
            iSample = iSample,
            sampleReads = sampleReads,
            method = method,
            K = K,
            priorCurrent = priorCurrent,
            alphaMatCurrent_t = alphaMatCurrent_t,
            eHapsCurrent_t = eHapsCurrent_t,
            transMatRate_t_H = transMatRate_t_H,
            transMatRate_t_D = transMatRate_t_D,            
            iteration = iteration,
            niterations = niterations,
            K_random = K_random,
            grid_eHaps_distance = grid_eHaps_distance,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            highCovInLow = highCovInLow,
            Jmax = Jmax,
            nor = nor,
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,
            srp = srp,
            pseudoHaploidModel = pseudoHaploidModel,
            blocks_for_output = blocks_for_output,
            generate_fb_snp_offsets = generate_fb_snp_offsets,
            return_a_sampled_path = return_a_sampled_path,
            grid = grid,
            update_in_place = TRUE,
            gammaUpdate_t = gammaSum_t, ## do switch here
            jUpdate_t = alphaMatSum_t,
            hapSum_t = hapSum_t,
            priorSum = priorSum,
            pass_in_alphaBeta = pass_in_alphaBeta,
            alphaHat_t = alphaHat_t,
            betaHat_t = betaHat_t
        )
        fbsoL <- out$fbsoL

        ## kind of like 
        if (iteration == niterations) {
            if (useTempdirWhileWriting) {
                alphaBetaBlocks <- lapply(fbsoL, function(x) x[["alphaBetaBlocks"]])
                save(alphaBetaBlocks, file = file_alphaBetaBlocks(tempdir, iSample, regionName), compress = FALSE)
                if (length(bundling_info) > 0) {
                    ## bundle here
                    last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
                    if (last_in_bundle == 1) {
                        bundle_inputs_after_generation(
                            bundling_info = bundling_info,
                            iBam = iSample,
                            dir = tempdir,
                            regionName = regionName,
                            what = "alphaBetaBlocks"
                        )
                    }
                }
            } else {
                allAlphaBetaBlocks[[iiSample]] <- lapply(fbsoL, function(x) x[["alphaBetaBlocks"]])
            }
        }
        
        out <- within_EM_per_sample_heuristics(
            sampleReads = sampleReads,
            iSample = iSample,
            iiSample = iiSample,
            K = K,
            fbsoL = fbsoL,
            method = method,
            regionName = regionName,
            L = L,
            L_grid = L_grid,
            nGrids = nGrids,
            grid = grid,
            nor = nor,
            K_subset = K_subset,
            nbreaks = nbreaks,
            break_results = break_results,
            fromMat = fromMat,
            srp = srp,
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,
            nSNPs = nSNPs,
            eHapsCurrent_t = eHapsCurrent_t,
            alphaMatCurrent_t = alphaMatCurrent_t,
            sigmaCurrent = sigmaCurrent,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            outputdir = outputdir,
            iteration = iteration,
            restartIterations = restartIterations,
            highCovInLow = highCovInLow,
            nsplit = nsplit,
            best_K_for_sample = best_K_for_sample,
            restartMatrixList = restartMatrixList,
            readsTotal = readsTotal,
            readsSplit = readsSplit,
            pseudoHaploidModel = pseudoHaploidModel,
            tempdir = tempdir,
            samples_with_phase = samples_with_phase,
            i_core = i_core,
            sampleRange = sampleRange,
            fbd_store = fbd_store,
            plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts,
            grid_distances = grid_distances,
            return_a_sampled_path = return_a_sampled_path,
            sampledPathList = sampledPathList
        )
        readsTotal <- out$readsTotal
        readsSplit <- out$readsSplit
        restartMatrixList <- out$restartMatrixList
        fromMat <- out$fromMat
        fbd_store <- out$fbd_store
        sampledPathList <- out$sampledPathList

    }

    if (method == "pseudoHaploid") {
        hapSum_t <- hapSum_t / 2
    }

    return(
        list(
            alphaMatSum_t = alphaMatSum_t,
            gammaSum_t = gammaSum_t,
            priorSum = priorSum,
            hapSum_t = hapSum_t,
            fromMat = fromMat,
            readsSplit = readsSplit,
            readsTotal = readsTotal,
            restartMatrixList = restartMatrixList,
            allAlphaBetaBlocks = allAlphaBetaBlocks,
            sampledPathList = sampledPathList
        )
    )
}




get_transMatRate <- function(method, sigmaCurrent) {
    if (method == "diploid") {
        x <- sigmaCurrent
        transMatRate_t <- rbind(x ** 2, x * (1 - x), (1 - x) ** 2)
    } else if ((method == "pseudoHaploid") | (method == "diploid-inbred")) {
        transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)
    } else {
        stop("bad method")
    }
    return(transMatRate_t)
}

check_mclapply_OK <- function(out, stop_message = "An error occured during STITCH. The first such error is above") {
    te <- sapply(out, class) == "try-error"
    if (sum(te) > 0) {
        print_message(out[[which(te)[1]]]) # print first error
        stop(stop_message)
    }
    return(NULL)
}


## just designed to start slower for larger regions to impute
## informally, has shown to improve imputation performance by making larger starting blocks
get_Jmax_wrt_iteration <- function(Jmax, iteration, niterations) {
    if (is.na(niterations) == FALSE & is.na(iteration) == FALSE) {
        if (niterations > 10) {
            if(iteration <= 10) {
                Jmax <- as.integer(iteration - 1)
            }
        }
    }
    return(Jmax)
}


completeSampleIteration <- function(N,tempdir,chr,K,K_subset, K_random, nSNPs, nGrids,priorCurrent,eHapsCurrent_t,alphaMatCurrent_t,sigmaCurrent,maxDifferenceBetweenReads,maxEmissionMatrixDifference, Jmax,aSW=NA,bSW=NA,highCovInLow,iteration,method,expRate,minRate,maxRate,niterations,splitReadIterations,shuffleHaplotypeIterations,nCores,L,nGen,alphaMatThreshold,emissionThreshold,gen,outputdir,environment,pseudoHaploidModel,outputHaplotypeProbabilities,switchModelIteration,regionName,restartIterations,refillIterations,outputBlockSize,bundling_info, alleleCount, phase, samples_with_phase, vcf.piece_unique, grid, grid_distances, L_grid, B_bit_prob, nSNPsInRegion, start_and_end_minus_buffer, shuffle_bin_nSNPs, shuffle_bin_radius, plot_shuffle_haplotype_attempts, blocks_for_output, allSampleReads, snps_in_grid_1_based, minimizeSwitchingIterations, useTempdirWhileWriting
) {

    print_message(paste0("Start of iteration ", iteration))
    if (is.na(match(iteration,restartIterations))==FALSE)
        print_message("Restart read probabilities")

    ##
    ## set up switching iteration (on shuffleHaplotypeIteration)
    ##
    out <- get_nbreaks(
        iteration = iteration,
        tempdir = tempdir,
        regionName = regionName,
        shuffleHaplotypeIterations = shuffleHaplotypeIterations,
        nGrids = nGrids,
        shuffle_bin_nSNPs = shuffle_bin_nSNPs,
        shuffle_bin_radius = shuffle_bin_radius
    )
    nbreaks <- out$nbreaks
    break_results <- out$break_results
    
    ## set up splitting iteration
    nsplit <- 0
    if(is.na(match(iteration,splitReadIterations))==FALSE) {
        nsplit <- 1
    }
    ##
    ## transition matrix amount
    ##
    transMatRate_t_H <- get_transMatRate(method = "diploid-inbred", sigmaCurrent = sigmaCurrent)    
    transMatRate_t_D <- get_transMatRate(method = "diploid", sigmaCurrent = sigmaCurrent)

    sampleRanges <- getSampleRange(N, nCores)    

    single_iteration_results <- mclapply(
        sampleRanges,
        mc.cores=nCores,
        FUN = subset_of_complete_iteration,
        tempdir=tempdir,chr=chr,K=K, K_subset = K_subset, K_random = K_random, nSNPs = nSNPs,nGrids = nGrids, priorCurrent=priorCurrent,eHapsCurrent_t=eHapsCurrent_t,alphaMatCurrent_t=alphaMatCurrent_t,sigmaCurrent=sigmaCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,maxEmissionMatrixDifference = maxEmissionMatrixDifference,Jmax=Jmax,highCovInLow=highCovInLow,iteration=iteration,method=method,nsplit=nsplit,expRate=expRate,minRate=minRate,maxRate=maxRate,gen=gen,outputdir=outputdir,pseudoHaploidModel=pseudoHaploidModel,outputHaplotypeProbabilities=outputHaplotypeProbabilities,switchModelIteration=switchModelIteration,regionName=regionName,restartIterations=restartIterations,refillIterations=refillIterations,outputBlockSize=outputBlockSize, bundling_info = bundling_info, transMatRate_t_H = transMatRate_t_H, transMatRate_t_D = transMatRate_t_D, sampleRanges = sampleRanges, N = N, niterations = niterations, L = L, samples_with_phase = samples_with_phase, nbreaks = nbreaks, break_results = break_results, vcf.piece_unique = vcf.piece_unique, grid = grid, B_bit_prob = B_bit_prob, start_and_end_minus_buffer = start_and_end_minus_buffer, plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts, blocks_for_output = blocks_for_output, allSampleReads = allSampleReads, L_grid = L_grid, grid_distances = grid_distances, minimizeSwitchingIterations = minimizeSwitchingIterations, useTempdirWhileWriting = useTempdirWhileWriting
    )

    check_mclapply_OK(single_iteration_results)

    if (iteration == niterations) {
        ## if this is the final iteration, just save forward backwards, possibly to disk
        ## then exit
        return(reorder_alphaBetaBlocks(single_iteration_results, sampleRanges, N, K, nGrids, useTempdirWhileWriting))
    } else {
        allAlphaBetaBlocks <- NULL
    }

    
    ##
    ## update results
    ##
    out <- calculate_updates(out2 = single_iteration_results, sampleRanges = sampleRanges, K = K, nSNPs = nSNPs, N = N, nGen = nGen, expRate = expRate, minRate = minRate, maxRate = maxRate, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, L = L, grid_distances = grid_distances, alleleCount = alleleCount)
    sigmaSum <- out$sigmaSum
    sigmaSum_unnormalized <- out$sigmaSum_unnormalized
    priorSum <- out$priorSum
    alphaMatSum_t <- out$alphaMatSum_t
    gammaSum_t <- out$gammaSum_t
    hapSum_t <- out$hapSum_t
    rm(out)

    ##
    ## get other updates
    ##
    out <- calculate_misc_updates(N = N, nsplit = nsplit, nbreaks = nbreaks, sampleRanges = sampleRanges, out2 = single_iteration_results, fromMat = fromMat, K = K, iteration = iteration, restartIterations = restartIterations)
    restartMatrixList <- out$restartMatrixList
    readsTotal <- out$readsTotal
    readsSplit <- out$readsSplit
    fromMat <- out$fromMat
    rm(out)


    ##
    ## set up switching iteration (on shuffleHaplotypeIteration)
    ##
    if (iteration %in% (shuffleHaplotypeIterations - 1)) {
        if (is.null(shuffle_bin_radius) == FALSE) {
            define_and_save_breaks_to_consider(
                tempdir = tempdir,
                regionName = regionName,
                sigmaSum_unnormalized = sigmaSum_unnormalized,
                L_grid = L_grid,
                grid_distances = grid_distances,
                nGrids = nGrids,
                nGen = nGen,
                minRate = minRate,
                maxRate = maxRate,
                iteration = iteration,                
                shuffle_bin_radius = shuffle_bin_radius,
                plot_shuffle_haplotype_attempts = plot_shuffle_haplotype_attempts
            )
        }
    }

    ##
    ##
    ## do haplotype shuffling
    ##
    if(nbreaks > 0) {
        out <- getBetterSwitchesSimple(
            fromMat = fromMat,
            nbreaks = nbreaks,
            break_results = break_results,
            K = K,
            eHapsFuture_t = gammaSum_t,
            alphaMatFuture_t = alphaMatSum_t,
            grid = grid,
            iteration = iteration,
            snps_in_grid_1_based = snps_in_grid_1_based
        )
        gammaSum_t <- out$eHapsFuture_t
        alphaMatSum_t <- out$alphaMatFuture_t
        whichIsBest <- out$whichIsBest
        rm(out)
        if (plot_shuffle_haplotype_attempts) {
            load(file = file_fbdStore(tempdir, regionName, iteration)) ## to load fbdStore
            plot_attempt_to_find_shuffles(
                grid_distances = grid_distances,
                L_grid = L_grid,
                fbd_store = fbd_store,
                tempdir = tempdir,
                outputdir = outputdir,
                regionName = regionName,
                iteration = iteration,
                whichIsBest = whichIsBest,
                is_reference = FALSE ## obvious!
            )
            ## for now, re-save everything to replot, etc
            ## print("REMOVE ME")
            ## save(fromMat, nbreaks, break_results, K, gammaSum_t, alphaMatSum_t, grid, snps_in_grid_1_based, grid_distances, L_grid, fbd_store, tempdir, outputdir, regionName, iteration, file = file_fbdStore(tempdir, regionName, iteration))
        }
    }
    ##
    ## look at refilling - round 2
    ##
    if (is.na(match(iteration, refillIterations)) == FALSE) {

        ## OK - try to make this more complicated
        ## previous iteration - generate paths
        ## this generation - figure out common paths from those

        
        print_message(paste0("Iteration - ",iteration," - refill infrequently used haplotypes"))
        out <- refillSimple(
            hapSum_t = hapSum_t,
            nGrids = nGrids,
            K = K,
            gammaSum_t = gammaSum_t,
            N = N,
            L_grid = L_grid,
            grid = grid
        )
        gammaSum_t<- out$gammaSum_t
        ever_changed <- out$ever_changed
        noise <- 0 ## consider adding in noise so things aren't exact
        sEC <- sum(ever_changed)
        if ((sEC > 2) & (noise > 0)) {
            gammaSum_t[, ever_changed] <-
                (1 - noise) * gammaSum_t[, ever_changed] +
                noise * array(runif(sEC * K), c(K, sEC))
            priorSum <- noise * rep(1 / K, K) + (1 - noise) * priorSum # restart these as well
            ## this is arguably the worst one
            alphaMatSum_t <-
                noise * matrix(1 / K, ncol = (nGrids - 1), nrow = K) +
                (1 - noise) * alphaMatSum_t
        }
        rm(out)
    }
    ##
    ## look at splitting
    ##
    if(nsplit > 0) {
        print_message(paste0(
            "Split reads, average N=",round(mean(readsSplit))," (",round(100*mean(readsSplit)/mean(readsTotal),3), " %)"
        ))
        save(readsTotal,readsSplit,file = file.path(outputdir, "RData", paste0("splitReads.",regionName,".iteration.",iteration,".RData")))
    }
    ##
    ## finally, if we have high coverage samples, run them
    ##
    if(length(highCovInLow)>0) {
        currentR2 <- getR2DuringEMAlgorithm(highCovInLow,gen=gen,tempdir=tempdir,regionName=regionName, alleleCount = alleleCount)
        ## output!
        write.table(
            matrix(c(iteration,round(currentR2,4),date()),nrow=1),
            file = file.path(outputdir, "RData", paste0("interim.r2",chr,".txt")),
            row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=iteration>1
        )
        if (iteration == 1)
            print_message("Printing out per-sample and mean estimates of correlation between provided genotypes from genfile and sample dosages, oriented so that the major allele has dosage 0 and the minor allele has dosage 1")
        y <- round(currentR2, 3)
        print_message(paste0(
            "iteration=", iteration, ", sample(r2)=",
            paste(y[-length(y)], collapse = ", "),
            " - mean=", y[length(y)]
        ))
    }
    ##
    ## also, if there is phasing data, calculate PSE
    ##
    if (is.null(phase) == FALSE & method == "pseudoHaploid") {
        pse <- get_phase_switch_error(
            samples_with_phase = samples_with_phase,
            tempdir = tempdir,
            regionName = regionName,
            phase = phase
        )
        if (iteration == 1)
            print_message("Printing out per-sample phase switch error")
        print_message(paste0(
            "iteration = ", iteration, ", sample(pse %)=",
            paste(round(pse, 2) * 100, collapse = ", "), ", "
        ))
    }
    ##
    ## also, if its a restart iteration, print out some stats
    ##
    ##if (iteration %in% restartIterations) {
    ##    restartMatrix <- cbind(rep(1:N,sapply(restartMatrixList,nrow)),matrix(unlist(lapply(restartMatrixList,t)),ncol=5,byrow=TRUE))
    ##    save(restartMatrix,file=paste0(outputdir,"RData/restartMatrix.iteration.",iteration,".",regionName,".RData"))
    ##}
    return(
        list(
            eHapsFuture_t = gammaSum_t,
            sigmaFuture = sigmaSum,
            priorFuture = priorSum,
            alphaMatFuture_t = alphaMatSum_t,
            hapSum_t = hapSum_t,
            allAlphaBetaBlocks = allAlphaBetaBlocks
        )
    )

}



calculate_misc_updates <- function(
    N,
    nsplit,
    nbreaks,
    sampleRanges,
    out2,
    fromMat,
    K,
    iteration,
    restartIterations
) {
    if(nbreaks > 0) {    
        fromMat <- array(0, c(nbreaks, K, K))
    } else {
        fromMat <- NULL
    }
    ## note - already defined fromMat before the function f above
    if(nsplit > 0) {    
        readsTotal <- array(0, N)
        readsSplit <- array(0, N)
    } else {
        readsTotal <- NULL
        readsSplit <- NULL
    }
    if(iteration %in% restartIterations) {    
        restartMatrixList <- as.list(1:N)
    } else {
        restartMatrixList <- NULL
    }
    for(i in 1:length(sampleRanges)) {
        if(nbreaks > 0) {
            fromMat <- fromMat + out2[[i]][["fromMat"]]
        }
        if(nsplit > 0) {
            readsTotal <- readsTotal + out2[[i]][["readsTotal"]]
            readsSplit <- readsSplit + out2[[i]][["readsSplit"]]
        }
        if(iteration %in% restartIterations) {            
            restartMatrixList[sampleRanges[[i]][1]:sampleRanges[[i]][2]]=out2[[i]]$restartMatrixList[sampleRanges[[i]][1]:sampleRanges[[i]][2]]
        }
    }
    return(
        list(
            restartMatrixList = restartMatrixList,
            readsTotal = readsTotal,
            readsSplit = readsSplit,
            fromMat = fromMat
        )
    )
}


calculate_updates <- function(
    out2, sampleRanges, K, nSNPs, N,
    nGen, expRate, minRate, maxRate,
    emissionThreshold, alphaMatThreshold, L,
    grid_distances, alleleCount = NA
) {

    nGrids <- ncol(out2[[1]]$alphaMatSum_t) + 1
    priorSum <- array(0, K)
    alphaMatSum_t <- array(0, c(K, nGrids - 1))
    gammaSumBoth_t <- array(0, c(K, nSNPs, 2))
    hapSum_t <- array(0, c(K, nGrids))
    ## sum and sum

    for(i in 1:length(sampleRanges)) {
        priorSum <- priorSum + out2[[i]]$priorSum
        alphaMatSum_t <- alphaMatSum_t + out2[[i]]$alphaMatSum_t
        gammaSumBoth_t <- gammaSumBoth_t + out2[[i]]$gammaSum_t
        hapSum_t <- hapSum_t + out2[[i]]$hapSum_t
    }

    ## normalize prior
    priorSum <- priorSum / sum(priorSum)
    ## note - it is possible in simulations not to have sampled the first SNP
    ## in which case there is no posterior information, so set to random
    if(is.na(priorSum[1]))
        priorSum <- rep(1 / K, K)

    minPriorSum <- (1 / K) / 100
    priorSum[priorSum < minPriorSum] <- minPriorSum
    priorSum <- priorSum / sum(priorSum)


    ## normalize sigma
    sigmaSum <- colSums(alphaMatSum_t) / N / 2
    sigmaSum <- exp(-sigmaSum)

    which <- is.na(sigmaSum)
    if (is.null(grid_distances)) {
        dl <- diff(L)
    } else {
        dl <- grid_distances
    }
    if (sum(which) > 0) {
        sigmaSum[which] <- exp(-nGen * expRate / 100 / 1000000 * dl)[which]
    }
    sigmaSum_unnormalized <- sigmaSum
    
    ## morgans per SNP assuming T=100, 0.5 cM/Mb
    x1 <- exp(-nGen * minRate * dl/100/1000000) # lower
    x2 <- exp(-nGen * maxRate * dl/100/1000000) # upper
    ## recombination average rate
    ## we estiamte the compound parameter T * sigma_t
    ## so we need to use nGen to bound apppropriately
    sigmaSum[sigmaSum > x1] <- x1[sigmaSum > x1]
    sigmaSum[sigmaSum < x2] <- x2[sigmaSum < x2]

    ## needed if only 1 SNP. stupid R and not able to drop selective dimensions
    gammaSum_t <- array(0, c(K, nSNPs))    
    gammaSum_t[, ] <- gammaSumBoth_t[, , 1] / gammaSumBoth_t[, , 2] 
    
    ## if missing, use alleleCount to fill in
    if (is.na(alleleCount[1])) {
        gammaSum_t[is.na(gammaSum_t)] <- 0.5 ## blank out entire SNP if no reads anywhere
    } else {
        for(k in 1:K) {
            w <- is.na(gammaSum_t[k, ])
            gammaSum_t[k, w] <- alleleCount[w, 3]
        }
    }
    
    gammaSum_t[gammaSum_t > (1 - emissionThreshold)] <- (1 - emissionThreshold)
    gammaSum_t[gammaSum_t < emissionThreshold] <- emissionThreshold

    ## honestly what was I even doing before
    ## this looks right, possibly inefficient
    alphaMatSum_t <- alphaMatSum_t / rep(colSums(alphaMatSum_t), each = K)
    
    ## now reset columns below threshold to value to make to sum to 1
    alphaMatSum_t[alphaMatSum_t < alphaMatThreshold] <- 0
    ## if this happens, reset the averages of those columns
    ## to keep each column having sum 1, with minimum entry the thresold value
    how_many_cols_below_0 <- colSums(alphaMatSum_t == 0)
    if (sum(how_many_cols_below_0) > 0) {
        ## for columns with an entry below 0
        ## each 0 entry becomes threshold
        ## then rest re-scaled so whole thing has sum 1
        which_cols <- how_many_cols_below_0 > 0
        y <- alphaMatThreshold * how_many_cols_below_0[which_cols]
        x <- colSums(alphaMatSum_t[, which_cols, drop = FALSE])
        for(k in 1:K) {
            alphaMatSum_t[k, which_cols] <- (1 - y) * alphaMatSum_t[k, which_cols] / x
        }
        alphaMatSum_t[alphaMatSum_t == 0] <- alphaMatThreshold
    }

    ## now - cols, meaning SNPs, want to
    t1 <- colSums(is.na(alphaMatSum_t)) > 0
    if (sum(t1) > 0) {
        warning("alphaMat update contains NA")
        alphaMatSum_t[, t1] <- 1 / K
    }

    return(
        list(
            sigmaSum = sigmaSum,
            sigmaSum_unnormalized = sigmaSum_unnormalized,            
            priorSum = priorSum,
            alphaMatSum_t = alphaMatSum_t,
            gammaSum_t = gammaSum_t,
            hapSum_t = hapSum_t
        )
    )
}



## for a sample, first, find regions where the two haplotypes are nearly the same
## next, loop over these regions. run diploid forward backward algorithm
## once that's done, re-do probabilities of the two haplotypes
restartSample <- function(sampleReads, srp, pRgivenH1, pRgivenH2,fbsoL, nSNPs,eHapsCurrent_t,sigmaCurrent,alphaMatCurrent_t,Jmax,priorCurrent,maxDifferenceBetweenReads,maxEmissionMatrixDifference,outputdir,iteration,restartIterations, method) {
    ##if(is.na(match(iteration, restartIterations)) == TRUE)
    ##    return(list(m = NULL, fbsoL = fbsoL))
    if (method== "diploid")
        stop("Attempting to restart iterations for diploid mode - this should not happen")
    ##
    ## first, find the regions we want to run this on
    ##
    m <- pseudoHaploidFindRegion(
        srp = srp,
        pRgivenH1 = pRgivenH1,
        pRgivenH2 = pRgivenH2
    )
    m=cbind(matrix(m,ncol=3),NA,NA)
    if(sum(m[,3]==1)==0) return(list(m=m,fbsoL=fbsoL))
    ##
    ## otherwise - continue
    ##
    xx=(1:nrow(m))[m[,3]==1]
    for(iR in 1:length(xx)) {
        start=m[xx[iR],1]
        end=m[xx[iR],2]
        ##
        ## now - over this region - re-run diploid
        ## get a viterbi path back
        ##
        ## buffer another 50 reads to either side
        ##
        start=max(1,start-50)
        end=min(length(sampleReads),end+50)
        TR=range(unlist(lapply(sampleReads[start:end],function(x) x[[4]])))+1 # 1-based SNPs to use
        ##
        ## subset reads
        ##
        sampleReads=sampleReads[start:end]
        sampleReads=lapply(sampleReads,function(x) {
            ## remember - these are 0-based
            x[[2]]=x[[2]]-(TR[1]-1)
            x[[4]]=x[[4]]-(TR[1]-1)
            return(x)
        })
        ##
        ## forward backward diploid here
        ##
        transMatRate_t_D <- get_transMatRate(
            "diploid",
            sigmaCurrent[TR[1]:(TR[2]-1)]
        )
        fbd <- forwardBackwardDiploid(
            sampleReads = sampleReads,
            nReads = as.integer(length(sampleReads)),
            pi = priorCurrent,
            eHaps_t = eHapsCurrent_t[, TR[1]:TR[2]],
            alphaMat_t = alphaMatCurrent_t[, TR[1]:(TR[2]-1)],
            transMatRate_t_D = transMatRate_t_D,
            Jmax = Jmax,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,            
            suppressOutput = as.integer(1),
            blocks_for_output = array(0, c(1, 1)),
            return_gamma = TRUE, ## unsure
            update_in_place = FALSE,
            gammaUpdate_t = array(0, c(1, 1, 1)),
            jUpdate_t = array(0, c(1, 1)),
            hapSum_t = array(0, c(1, 1)),
            priorSum = array(0, 1),
            pass_in_alphaBeta = FALSE,        
            alphaHat_t = array(0, c(1, 1)),
            betaHat_t = array(0, c(1, 1))
        )
        fbd$gamma <- t(fbd$gamma_t)
        fbd$jUpdate <- t(fbd$jUpdate_t)
        ##
        ## compare to haploid run, try to properly integrate
        ##
        a1=fbsoL[[1]]$gamma[c(TR[1],TR[2]),]
        a2=fbsoL[[2]]$gamma[c(TR[1],TR[2]),]
        b1=fbd$gamma[c(1,TR[2]-TR[1]+1),]
        ##
        ## simplest - argmax. symmetry shouldn't be a problem as argmax has non random tie breaks
        ##
        z=apply(fbd$gamma,1,which.max)
        k1=z%%K
        k1[k1==0]=K
        k2=(z-k1)/K+1
        ##
        ## try to determine right way to switch together
        ##
        startD=c(k1[1],k2[1])
        endD=c(k1[length(k1)],k2[length(k2)])
        ##
        ## now - either determine whether c(hap1 ,h
        ##
        ## first option being bigger - means k1 matches hap1, k2 matches hap2
        ## otherwhise, opposite
        wS=c(a1[1,startD[1]]+a2[1,startD[2]],a2[1,startD[1]]+a1[1,startD[2]])
        wSX=which.max(wS)
        ## now look at the end
        wE=c(a1[2,endD[1]]+a2[2,endD[2]],a2[2,endD[1]]+a1[2,endD[2]])
        wEX=which.max(wE)
        ##
        ## if max or wE and wS are the same, we're OK, no switch required
        ##
        ## no matter what, push updates to gamma, based on new states
        ## wS is all that matters
        ##
        for(iHap in 1:2) {
            if(wSX==1 & iHap==1) kL=k1
            if(wSX==1 & iHap==2) kL=k2
            if(wSX==2 & iHap==1) kL=k2
            if(wSX==2 & iHap==2) kL=k1
            p=0.1 # "error" split over other states
            q=array(runif((1+diff(TR))*K),c(1+diff(TR),K)) # the new gammas for that region
            q=p * q/rowSums(q)
            q[cbind(1:nrow(q),kL)]=(1-p) + q[cbind(1:nrow(q),kL)]
            ## changes to gamma for 1 are based on k1
            fbsoL[[iHap]]$gamma[TR[1]:TR[2],]=q
        }
        if(wSX!=wEX & TR[2]!=nSNPs) # if they are different, and there's something left to do
        {
            a1=fbsoL[[1]]$gamma[(TR[2]+1):nSNPs,]
            a2=fbsoL[[2]]$gamma[(TR[2]+1):nSNPs,]
            fbsoL[[2]]$gamma[(TR[2]+1):nSNPs,]=a2
            fbsoL[[1]]$gamma[(TR[2]+1):nSNPs,]=a1
        }
        ##
        ## out of interest, record whether a flip has been done
        ##
        m[xx[iR],4:5]=c(wSX,wEX)
    } # end of loop over region
    return(list(m=m,fbsoL=fbsoL))
}










pseudoHaploidFindRegion=function(srp,pRgivenH1,pRgivenH2)
{
  # in this function, we define underperforming regions
  #srp,pRgivenH1,pRgivenH2
  # define an underperforming region as a set of >100 SNPs with av <0.0001 difference
    n=100 # region to check over
    n2=100 # when a region is smaller than this, remove it
    x=abs(pRgivenH1-pRgivenH2)
    y=cumsum(x)
    p=y[-c(1:(n-1))]-y[-(length(y)-(n-2):0)]
    w=p<1
    # turn into a set of regions
    z=array(FALSE,length(pRgivenH1))
    z[c(array(FALSE,n-1),w)]=TRUE
    z[c(w,array(FALSE,n-1))]=TRUE
    #z=array(FALSE,length(pRgivenH1))
    #z[1000:1200]=TRUE
    #z[1250:1600]=TRUE
    #z[1670:1800]=TRUE
    #z[3000:3500]=TRUE
    start=c(1,(1:length(z))[diff(z)!=0]+1)
    end=c((1:length(z))[diff(z)!=0],length(pRgivenH1))
    state=z[start]
    # recap
    # have length(pRgivenH1) reads
    # those are now in regions, which start at "start" and end at "end"
    # inside those regions, they all have the same state as "state"
    #
    # if there is one state - return
    #
    if(length(state)==1)
      return(c(start,end,state))
    #
    # keep removing the smallest state until none are less than n2 in length
    #
    m=cbind(start,end,state)
    while(sum((m[,2]-m[,1])<n2)>0)
    {
      q=m[,2]-m[,1]
      l=which.min(q) # remove this sample
      # if there are only two entries, remove smallest
      if(nrow(m)==2)
      {
        mX=m
        m=matrix(m[-l,],ncol=3)
        if(l==2) m[1,2]=mX[2,2]
        if(l==1) m[1,1]=mX[1,1]
      } else {
        if(l==1)
        {
          m=m[-1,]
          m[1,1]=1
        }
        if(l==nrow(m))
        {
          a=m[nrow(m),2]
          m=m[-nrow(m),]
          m[nrow(m),2]=a
        }
        if(l>1 & l<nrow(m))
        {
          # remove two entries here
          mX=m
          m=matrix(m[-c(l,l+1),],ncol=3)
          m[l-1,2]=mX[l+1,2]
        }
      }
    }
    return(m)
}










# derived from "A Note on Exact Tests of Hardy-Weinberg Equilibrium"
calculate_hwe_p <- function(x) {
  #x <- c(100,50,10)

    if (x[3]>x[1])
        x <- x[c(3,2,1)]

    nAA <- x[1] # individuals
    nAB <- x[2] # individuals
    nBB <- x[3] # individuals
    nA <- nAA * 2 + nAB # number of copies of A allele
    nB <- nBB * 2 + nAB # number of copies of B allele
    n <- nAA + nAB + nBB # number of (diploid) individuals

    min_het <- 0
    max_het <- nAB + 2 * min(nAA, nBB)

    ## maximum-ish value
    mid <- floor((nA * nB) / ((nA + nB)))
    if (is.na(mid))
        mid <- 0
    ## make odd if it needs to be
    if ((nA%%2) != (mid%%2))
        mid <- mid + 1

    ## determine a mid point
    probs <- array(0, max_het + 1) # note - this is 0-based
    probs[mid + 1] <- 1

    ## use recurrence relation - going down
    n_het <- mid
    n_hom_alt <- (nBB * 2 + nAB - n_het) / 2
    n_hom_ref <- (n - n_het - n_hom_alt)
    if ((mid - 2) >= min_het) {
        for (het in seq(mid - 2, min_het, -2)) {
            probs[het + 1] <- probs[het + 3] *
                n_het * (n_het-1) / (4 * (n_hom_ref + 1) * (n_hom_alt + 1))
            n_het <- n_het - 2
            n_hom_ref <- n_hom_ref + 1
            n_hom_alt <- n_hom_alt + 1
        }
    }


    ## use recurrence relationship - going up
    n_het <- mid
    n_hom_alt <- (nBB * 2 + nAB - n_het) / 2
    n_hom_ref <- (n - n_het - n_hom_alt)
    if ((mid + 2) <= max_het) {
        for (het in seq(mid + 2, max_het, 2)) {
            probs[het + 1] <- probs[het - 1] *
                (4 * n_hom_ref * n_hom_alt) / ( (n_het + 2) * (n_het + 1))
            n_het <- n_het + 2
            n_hom_ref <- n_hom_ref - 1
            n_hom_alt <- n_hom_alt - 1
        }
    }


    all_probs <- probs / sum(probs)
    p2 <- sum(all_probs[all_probs <= all_probs[nAB + 1]])

    p <- min(1, p2)

    return(p)
}




generate_hwe_on_counts <- function(
    hweCount,
    nSNPs,
    nCores
) {
    snpRanges <- getSampleRange(nSNPs, nCores) ## should work on this as well
    out <- mclapply(snpRanges, mc.cores = nCores, function(i) {
        hwe_out <- array(NA, i[2] - i[1] + 1)
        for(j in i[1]:i[2]) {
            hwe_out[j - i[1] + 1] <- rcpp_calculate_hwe_p(hweCount[j, ])
        }
        return(hwe_out)
    })
    out <- unlist(out)
    return(out)
}



## function to get phred scaled likelihood and genotype
## from input RData file
get_pl_and_rd <- function(
    sampleReads,
    nSNPs
) {
    ## get genotype likelihoods
    x <- unlist(lapply(sampleReads,function(x) x[[4]])) + 1
    bq <- unlist(lapply(sampleReads,function(x) x[[3]]))
    s <- as.integer(bq > 0)
    bq <- abs(bq)
    e <- 10^(-bq/10)
    p <- cbind(1 - e, e * (1 / 3))
    pA <- p[cbind(1:nrow(p), s + 1)]
    pB <- p[cbind(1:nrow(p), 2 - s)]
    ## hmm, just do for loop? not ideal, but OK?
    g <- array(1,c(nSNPs,3))
    rd <- array(0,c(nSNPs,2))
    for (i in 1:length(x)) {
        a <- x[i]
        ##b=y[i]
        b <- pA[i]/2
        c <- pB[i]/2
        g[a, ] <- g[a, ] * c(b + b, b + c, c + c)
        z <- which.max(c(pA[i], pB[i]))
        rd[a,z] = rd[a,z] + 1
    }
    g <- g / rowSums(g) # divide by max
    a <- g[, 1]
    w <- a < g[,2]; a[w] = g[w,2]
    w <- a < g[,3]; a[w] = g[w,3]
    g <- g / a
    ## 0/0:2,0:2:3:0,3,45
    ## turn into PLs
    ## ex: 1/1:0,1:1:3:32,3,0
    pl <- round(10 * -log10(g) )
    return(
        list(
            pl = pl,
            rd = rd
        )
    )
}








## get the r2's for all high coverage sample and the average
getR2DuringEMAlgorithm <- function(
    highCovInLow,
    gen,
    tempdir,
    regionName,
    alleleCount
) {
    n_hc <- length(highCovInLow)
    out <- lapply(1:n_hc, function(j) {
        load(file = file_dosages(tempdir, highCovInLow[j], regionName))
        return(dosage)
    })
    ## now, get r2 between generated and observed, store and output
    store <- array(0, n_hc + 1)
    to_flip <- alleleCount[, 3] > 0.5
    for(isample in 1:n_hc) {
        x <- gen[, isample]
        y <- out[[isample]] ## dosage here is 0-2 now
        ## print_message(paste0("range(x) = ", range(x, na.rm=TRUE)))
        ## print_message(paste0("range(y) = ", range(y, na.rm=TRUE)))
        x[to_flip] <- 2 - x[to_flip]
        y[to_flip] <- 2 - y[to_flip]
        store[isample] <- suppressWarnings(cor(x, y, use = "complete.obs") ** 2)
    }
    store[n_hc + 1]=mean(store[1:n_hc])
    return(store)
}


## assign this when highCovInLow assigned
get_phase_switch_error <- function(
    samples_with_phase,
    tempdir,
    regionName,
    phase
) {
    n_phase <- length(samples_with_phase)
    out <- sapply(1:samples_with_phase, function(i) {
        iBam <- samples_with_phase[i]
        load(file = file_haplotypes(tempdir, iBam, regionName))
        ## per-sample estimate of phase switch error here
        pse <- calculate_pse(
            test = haplotypes,
            truth = phase[, i, ]
        )
        return(pse)
    })
    return(out)
}


## calculate phase switch error
## for now, assume everything is an integer
## and rowSums == 1

#' @export
calculate_pse <- function(
    test,
    truth,
    seed = NULL
) {
    ## for testing
    if (is.null(seed) == FALSE)
        set.seed(seed)
    which_sites <-
        rowSums(truth == 0 | truth == 1) == 2 &
        rowSums(truth) == 1
    truth <- truth[which_sites, ]
    test <- test[which_sites, ]
    if (nrow(test) == 0)
        return(NA)
    ## round test data for now. choose hets at random
    ## for homs
    test[, 1] <- as.integer(round(test[, 1]))
    test[, 2] <- as.integer(round(test[, 2]))
    choose_at_random <- which(rowSums(test) != 1)
    if (length(choose_at_random) > 0) {
        test[choose_at_random, ] <- 0
        r <- sample(
            c(1, 2),
            length(choose_at_random),
            replace = TRUE
        )
        test[cbind(choose_at_random, r)] <- 1
    }
    ## chose best start
    if (test[1, 1] != truth[1, 1])
        test <- test[, c(2, 1)]
    ## calculate number of differences
    n_bad <- sum(diff(test[,1] - truth[,1]) != 0)
    return(c(n_bad / (nrow(test) - 1)))
}


## plot HWE, info, MAF, coloured by discrepency
plotMetricsForPostImputationQC <- function(
    iSample,
    highCovList,
    gen,
    gen_imp,
    alleleCount,
    chr,
    L,
    estimatedAlleleFrequency,
    info,
    outputdir,
    colour = TRUE,
    hwe,
    regionName
) {
    if (is.null(iSample)==FALSE) {
        m1A <- cbind(gen[,iSample] / 2, gen_imp[, iSample])
        m1A[alleleCount[,3]>0.5,] <- 1 - m1A[alleleCount[,3]>0.5,]
        dist <- abs(m1A[,1]-m1A[,2])
        dist[dist > 0.1] <- 0.1
    } else {
        dist <- rep(0, T)
    }
    colfunc <- colorRampPalette(c("blue", "red"))
    col=colfunc(11)[round(dist,2)*100+1]
    ## plot
    jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQC.",regionName,".sample",iSample,".jpg")),height=1000,width=3000,qual=100)
    par(mfrow=c(1,3))
    plot(info,log10(hwe) ,cex.lab=3 ,col=col)
    plot(info,estimatedAlleleFrequency,cex.lab=3,col=col)
    plot(estimatedAlleleFrequency,log10(hwe) ,cex.lab=3,col=col )
    dev.off()
    ## now plot along the chromosome as well
    jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQCChromosomeWide.",regionName,".sample",iSample,".jpg")),height=3000,width=3000,qual=100)
    par(mfrow=c(3,1))
    plot(L,log10(hwe) ,cex.lab=3 ,col=col)
    plot(L,info,cex.lab=3,col=col)
    plot(L,estimatedAlleleFrequency,cex.lab=3,col=col )
    dev.off()
}





### plot estimated against actual allele frequency
plotEstimatedAgainstReal <- function(outputdir,alleleCount,estimatedAlleleFrequency,which,chr,regionName) {
    m1 <- cbind(alleleCount[, 3], estimatedAlleleFrequency)
    jpeg(file.path(outputdir, "plots", paste0("r2.",regionName,".goodonly.jpg")),height=800,width=800,qual=100)
    m1 <- m1[which, ]
    corr <- suppressWarnings(
        round(cor(m1[,1], m1[,2], use="complete.obs") ** 2, 3)
    )
    main <- paste("nSnps: ",dim(m1)[1],", r2:",corr,sep="")
    plot(x=m1[,1],y=m1[,2],xlim=c(0,1),ylim=c(0,1),main=main,xlab="Real Allele Frequency",ylab="Estimated allele frequency",pch=20,cex.main=1.5,cex.lab=1.5)
    abline(0,1)
    dev.off()
}



split_function <- function(
    eHapsCurrent_t,
    sampleRead,
    K
) {
    y <- cbind(
        matrix(
            eHapsCurrent_t[, sampleRead[[4]] + 1],
            ncol = K,
            byrow = TRUE
        ),
        1 / 2
    )
    out <- sapply(
        1:(K+1),
        function(k) {
        bq <- convertScaledBQtoProbs(sampleRead[[3]])
        exp(sum(log(bq[, 2] * y[, k] + bq[, 1] * (1 - y[, k]))))
    })
    return(out)
}


## returns 1-based index of what reads are not desired
get_reads_worse_than_50_50 <- function(sampleReads, eHapsCurrent_t, K) {
    ## only sample those with more than 2 SNPs
    whichReads <- which(sapply(sampleReads,function(x) x[[1]]) > 2)
    ## pre-allocate output
    ## do not duplicate sampleReads continuously
    to_return <- array(FALSE, length(sampleReads))
    for(iRead in whichReads) {
        m <- split_function(eHapsCurrent_t, sampleReads[[iRead]], K = K)
        if (which.max(m) == (K + 1))
            to_return[iRead] <- TRUE
    }
    return(which(to_return))
}



split_a_read <- function(
    sampleReads,
    read_to_split,
    gammaK_t,
    L,
    eHapsCurrent_t,
    K,
    grid
) {

    did_split <- FALSE
    sampleRead <- sampleReads[[read_to_split]]

    ## get the likely for and after states
    ## by checking, a few reads up and downstream, what changes
    ## note - probably want to change this to per-base, not per-read!
    startRead <- max(1, read_to_split - 10)
    endRead <- min(length(sampleReads), read_to_split + 10)
    ##
    startGammaGrid <- sampleReads[[startRead]][[2]] + 1
    endGammaGrid <- sampleReads[[endRead]][[2]] + 1

    change <- apply(gammaK_t[, c(startGammaGrid, endGammaGrid)], 1, diff)
    ## from is the one that drops the most
    ## to is the one that gains the most
    from <- which.min(change)
    to <- which.max(change)

    ## try a break between all SNPs - take the best one
    y <- eHapsCurrent_t[, sampleRead[[4]] + 1]

    bqProbs <- convertScaledBQtoProbs(sampleRead[[3]])
    g1 <- bqProbs[, 2] * y[from, ] +
          bqProbs[, 1] * (1 - y[from, ])
    g2 <- bqProbs[, 2] * y[to, ] +
          bqProbs[, 1] * (1 - y[to, ])

    ## sum probabilities
    ## best split is between SNPs
    z <- sapply(1:(sampleRead[[1]]),function(a) {
        x1 <- sum(g1[1:a])
        x2 <- sum(g2[(a+1):(sampleRead[[1]] + 1)])
        sum(x1 + x2)
    })
    best_split <- which.max(z) ## between this and + 1

    u <- log(
        split_function(eHapsCurrent_t, sampleRead, K)[ K + 1]
    )

    if (max(z)>u) {
        new_read_1 <- get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, 1, best_split, L, grid
        )
        new_read_2 <- get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, best_split + 1, sampleRead[[1]] + 1, L, grid
        )
        sampleReads[[read_to_split]] <- new_read_1
        sampleReads <- append(
            sampleReads,
            list(new_read_2)
        )
        did_split <- TRUE
    }

    return(
        list(
            sampleReads = sampleReads,
            did_split = did_split
        )
    )
}





findRecombinedReadsPerSample <- function(
    gammaK_t,
    eHapsCurrent_t,
    K,
    L,
    iSample,
    sampleReads,
    tempdir,
    regionName,
    grid,
    verbose=FALSE
) {
    ## needs a full run
    ## only do for some - need at least 3 SNPs to consider
    w <- get_reads_worse_than_50_50(
        sampleReads = sampleReads,
        eHapsCurrent_t = eHapsCurrent_t,
        K = K
    )
    w <- w[w != 1 & w != length(w)]
    count <- 0
    if (length(w) > 0) {
        for (w1 in w) {
            out <- split_a_read(
                sampleReads = sampleReads,
                read_to_split = w1,
                gammaK_t = gammaK_t,
                L = L,
                eHapsCurrent_t = eHapsCurrent_t,
                K = K,
                grid = grid
            )
            sampleReads <- out$sampleReads
            count <- count + as.integer(out$did_split)
        } # end of loop on reads
        sampleReads <- sampleReads[order(unlist(lapply(sampleReads,function(x) x[[2]])))]
        save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName), compress = FALSE)
        if (verbose==TRUE)
            print_message(paste0(
                "sample ", iSample, " readsSplit ", count, " readsTotal ", length(sampleReads)
            ))
    }
    return(
        list(
            readsSplit = count,
            readsTotal = length(sampleReads)
        )
    )
}


pseudoHaploid_update_9 <- function(
    pRgivenH1,
    pRgivenH2,
    eMatHap_t1,
    eMatHap_t2,
    gamma_t1,
    gamma_t2,
    K,
    srp
) {
    x1 <- pRgivenH1 / (pRgivenH1 + pRgivenH2)
    x2 <- pRgivenH2 / (pRgivenH1 + pRgivenH2)
    y1 <- eMatHap_t1
    for(k in 1:K)
        y1[k, ] <- (y1[k, ] - x2 * pRgivenH2) / x1
    y2 <- eMatHap_t2
    for(k in 1:K)
        y2[k, ] <- (y2[k, ] - x1 * pRgivenH1) / x2
    pRgivenH1 <- colSums(gamma_t1[, srp + 1] * y1)
    pRgivenH2 <- colSums(gamma_t2[, srp + 1] * y2)
    return(
        list(
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2
        )
    )
}



pseudoHaploidUpdate <- function(
    pRgivenH1,
    pRgivenH2,
    fbsoL,
    pseudoHaploidModel,
    srp,
    K
) {
    ##
    ## once both haplotypes are available - update omega immediately
    ##
    ## for each read, get new phasing
    ## srp comes in with readProbs originally - doesnt change
    ## 4 things here
    ## 1 - pRgivenH1 - P(R | h1)
    ## 2 - pRgivenH2 - P(R | h2)
    ## note that
    ## 3 - pRandH1 - P(R and H1) = 1/2 P(R| H1)
    ## 4 - pRandH2 - P(R and H2) = 1/2 P(R| H2)
    ##
    ## MODEL 1 - divide one by other (y1/(y1+y2)) to get readProbs
    ##
    ## readProbs is y1/(y1+y2) from above
    ## then probabilities using
    ##    eMatHap(iRead,k) = eMatHap(iRead,k) * readProbs(iRead) + (1/Kd)  * (1-readProbs(iRead));
    ## then updating using multiply both by readProbs
    ## y1=(fbsoL[[1]]$gamma[srp+1,] * fbsoL[[1]]$eMatHap)
    ## y2=(fbsoL[[2]]$gamma[srp+1,] * fbsoL[[2]]$eMatHap)
    ## readProbs=y1/(y1+y2)
    ##
    ## MODEL 3, 7 - try to predict which haplotype it came from
    ##
    if(pseudoHaploidModel==7 | pseudoHaploidModel==9) {
        out <- pseudoHaploid_update_model_9(
            pRgivenH1 = pRgivenH1,
            pRgivenH2 = pRgivenH2,
            eMatHap_t1 = fbsoL[[1]]$eMatHap_t,
            eMatHap_t2 = fbsoL[[2]]$eMatHap_t,
            gamma_t1 = fbsoL[[1]]$gamma_t,
            gamma_t2 = fbsoL[[2]]$gamma_t,
            K = K,
            srp = srp
        )
        pRgivenH1 <- out$pRgivenH1
        pRgivenH2 <- out$pRgivenH2
        ## bound above below
        pRgivenH1[pRgivenH1<0.001]=0.001
        pRgivenH2[pRgivenH2<0.001]=0.001
        pRgivenH1[pRgivenH1>0.999]=0.999
        pRgivenH2[pRgivenH2>0.999]=0.999
    }
    ##
    ## MODEL 8
    ##
    if(pseudoHaploidModel==8 | pseudoHaploidModel==10) {
        x=pRgivenH1/(pRgivenH1+pRgivenH2)
        ## get eMatHap without the second part
        y1=fbsoL[[1]]$eMatHapOri # shouldn't change 1 or 2
        pRgivenH1=rowSums(fbsoL[[1]]$gamma[srp+1,] * y1)
        pRgivenH2=rowSums(fbsoL[[2]]$gamma[srp+1,] * y1)
        ## bound above below
        pRgivenH1[pRgivenH1<0.001]=0.001
        pRgivenH2[pRgivenH2<0.001]=0.001
        pRgivenH1[pRgivenH1>0.999]=0.999
        pRgivenH2[pRgivenH2>0.999]=0.999
    }
    return(list(pRgivenH1=pRgivenH1,pRgivenH2=pRgivenH2))
}




estimate_read_proportions <- function(
    sampleReads,
    pRgivenH1,
    pRgivenH2,
    nSNPs
) {
    output <- array(0, c(nSNPs, 4))
    colnames(output) <- c("ER1", "EA1", "ER2", "EA2")
    norm <- pRgivenH1 + pRgivenH2
    p1_norm<- pRgivenH1 / norm
    p2_norm <- pRgivenH2 / norm
    for(iRead in 1:length(sampleReads)) {
        sampleRead <- sampleReads[[iRead]]
        bq <- sampleRead[[3]]
        u <- sampleRead[[4]]
        pr <- convertScaledBQtoProbs(bq)
        output[u + 1, ] <- output[u + 1, ] +
            c(
                p1_norm[iRead] * pr[, 1],
                p1_norm[iRead] * pr[, 2],
                p2_norm[iRead] * pr[, 1],
                p2_norm[iRead] * pr[, 2]
            )
    }
    return(output)
}


print_message <- function(x, include_mem = FALSE) {
    if (include_mem) {
        mem <- system("ps auxww | grep 'scripts/profile.R' | grep slave | grep -v 'grep' | awk -v OFS='\t' '$1=$1' | cut -f6", intern = TRUE)
        if (length(mem) > 0) {
            mem <- paste0(paste0(round(as.integer(mem) / 2 ** 20, 3), ollapse = ", "), " - ")
        } else {
            mem <- ""
        }
    } else {
        mem <- ""
    }
    message(
        paste0(
            "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", mem, x
        )
    )
}


get_reads_at_SNP <- function(sampleReads, nSNPs) {
    srp <- sapply(sampleReads, function(x) x[[2]])
    list <- tapply(1:length(srp), srp, function(x) x - 1, simplify = FALSE)
    ## add in missing components
    reads_at_SNPs <- as.list(rep(-1, nSNPs))
    reads_at_SNPs[match(names(list), 0:(nSNPs - 1))] <- list
    return(reads_at_SNPs)
}


snap_sampleReads_to_grid <- function(sampleReads, grid) {
    sampleReads <- lapply(sampleReads, function(sampleRead) {
        sampleRead[[2]] <- grid[sampleRead[[2]] + 1]
        return(sampleRead)
    })
    return(sampleReads)
}


## assign each of the positions to a grid
## input is numbers, e.g. 3, 5, 10, 15
## and a windowSize, like 5
## output is 1-based on grid coordinates, like
## 1-5 -> 0, 6-10 -> 1, etc
## remove holes, sigma will be made able to handle it with bounding
assign_positions_to_grid <- function(
    L,
    gridWindowSize
) {
    if (is.na(gridWindowSize) == FALSE) {
        grid <- ceiling(L / gridWindowSize)
        grid <- grid - min(grid)
        ## for L_grid, get first mid-point
        L_grid_start <- gridWindowSize * (ceiling(L[1] / gridWindowSize) - 0.5)
        grid_distances <- diff(unique(grid)) * gridWindowSize
        L_grid <- L_grid_start + c(0, cumsum(grid_distances))
        grid <- match(grid, unique(grid)) - 1
        nGrids <- length(grid_distances) + 1
        if ((length(L_grid) - length(grid_distances)) != 1) {
            stop("An error has been made assigning SNP positions to grid. Please report this")
        }
        ## this allows one to go from grid to start and end of SNPs in that grid
        ## so e.g. the fifth grid with 0-based index 4
        ## has 1-based SNPs starting from
        ## snps_in_grid_1_based[4 + 1, "grid_starts"]
        ## to
        ## snps_in_grid_1_based[4 + 1, "grid_ends"]        
        snps_in_grid_1_based <- cbind(
            snps_start = match(unique(grid), grid),
            snps_end = length(grid) - match(unique(grid), grid[length(grid):1]) + 1
        )
    } else {
        grid <- 0:(length(L) - 1)
        grid_distances <- diff(L)
        L_grid <- L
        nGrids <- length(L)
        snps_in_grid_1_based <- cbind(
            snps_start = 1:nGrids,
            snps_end = 1:nGrids
        )
    }
    return(
        list(
            grid = grid,
            grid_distances = grid_distances,
            L_grid = L_grid,
            nGrids = nGrids,
            snps_in_grid_1_based = snps_in_grid_1_based
        )
    )
}



snap_reads_to_grid <- function(
    N,
    nCores,
    regionName,
    tempdir,
    bundling_info,
    grid,
    downsampleToCov,
    sampleNames,
    allSampleReads = NULL,
    whatKindOfReads = "sampleReads",
    verbose = TRUE
) {

    if (verbose) {
        print_message("Snap reads to grid")
    }
    sampleRanges <- getSampleRange(N = N, nCores = nCores)

    out <- mclapply(
        sampleRanges,
        mc.cores = nCores,
        tempdir = tempdir,
        regionName = regionName,
        grid = grid,
        bundling_info = bundling_info,
        whatKindOfReads = whatKindOfReads,
        downsampleToCov = downsampleToCov,
        verbose = verbose,
        sampleNames = sampleNames,
        allSampleReads = allSampleReads,
        N = N,
        FUN = snap_reads_to_grid_subfunction
    )
    check_mclapply_OK(out, "There has been an error snapping reads to grid. Please see above")

    ## 
    allSampleReadsStats <- array(NA, c(N, 2))
    if (is.null(allSampleReads) == FALSE) {
        allSampleReads <- as.list(1:N)
    }
    
    for(i_core in 1:length(out)) {
        sampleRange <- sampleRanges[[i_core]]
        ## only get these stats for sampleReads, not the reference ones
        if (whatKindOfReads == "sampleReads") {
            y <- out[[i_core]]$allSampleReadsStats            
            for(iSample in sampleRange[1]:sampleRange[2]) {        
                allSampleReadsStats[iSample, ] <- y[[iSample]]
            }
        }
        ## re-capture
        if (is.null(allSampleReads) == FALSE) {
            x <- out[[i_core]]$allSampleReadsOutput
            for(iSample in sampleRange[1]:sampleRange[2]) {        
                allSampleReads[[iSample]] <- x[[iSample]]
            }
        }
    }

    if (verbose & (whatKindOfReads == "sampleReads")) {
        x1 <- mean(allSampleReadsStats[, 1])
        x2 <- mean(allSampleReadsStats[, 2])        
        y <- nrow(allSampleReadsStats)
        print_message(paste0(
            "Warning - When assigning reads to grids, ",
            sum(allSampleReadsStats[, 1] > 0), " out of ", y, " ",
            "samples had reads removed, with on average (across all samples) ",
            round(x1, 3), " / ", round(x2, 3), " ",
            "(", round(x1 / x2 * 100, 3), " %) of reads removed"
        ))
    }

    if (verbose) {
        print_message("Done snap reads to grid")
    }
    
    return(allSampleReads)

}


snap_reads_to_grid_subfunction <- function(
    sampleRange,
    tempdir,
    regionName,
    grid,
    bundling_info,
    downsampleToCov,
    sampleNames,
    allSampleReads,
    N,
    whatKindOfReads = "sampleReads",
    verbose = TRUE
) {
    
    bundledSampleReads <- NULL
    allSampleReadsOutput <- as.list(1:N)
    allSampleReadsStats <- as.list(1:N)
    
    for(iSample in sampleRange[1]:sampleRange[2]) {
        
        out <- get_sampleReads_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = iSample,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads,
            what = whatKindOfReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads

        ## snap to grid
        sampleReads <- snap_sampleReads_to_grid(
            sampleReads = sampleReads,
            grid = grid
        )

        if (whatKindOfReads == "sampleReads") {
            ## OK, need new downsample ARGH
            ## basically like the old one
            ## downsample again once on the grid
            ## no point imputing with 100 reads in a grid window
            
            out <- downsample_snapped_sampleReads(
                sampleReads = sampleReads,
                iBam = iSample,
                downsampleToCov = downsampleToCov,
                sampleNames = sampleNames,
                verbose = verbose
            )
            sampleReads <- out$sampleReads
            allSampleReadsStats[[iSample]] <- out$remove_stats
        }

        ## then, either add to allSampleReads
        ## or, save and possibly bundle
        if (is.null(allSampleReads) == FALSE) {
            if (whatKindOfReads == "sampleReads") {
                allSampleReadsOutput[[iSample]] <- sampleReads
            }
        } else {
            if (whatKindOfReads == "sampleReads") {
                save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName), compress = FALSE)
            } else if (whatKindOfReads == "referenceSampleReads") {
                save(sampleReads, file = file_referenceSampleReads(tempdir, iSample, regionName), compress = FALSE)
            }
            if (length(bundling_info) > 0) {
                last_in_bundle <- bundling_info$matrix[iSample, "last_in_bundle"]
                if (last_in_bundle == 1) {
                    bundle_inputs_after_generation(
                        bundling_info = bundling_info,
                        iBam = iSample,
                        dir = tempdir,
                        regionName = regionName,
                        what = whatKindOfReads
                    )
                }
            }
        }
    }
    return(
        list(
            allSampleReadsOutput = allSampleReadsOutput,
            allSampleReadsStats = allSampleReadsStats
        )
    )
}






## remove reads as necessary
## do this per-grid point
## so e.g. <= 50 reads per site
downsample_snapped_sampleReads <- function(
    sampleReads,
    iBam,
    downsampleToCov,
    sampleNames,
    verbose = TRUE,
    print_warning_to_screen = FALSE
) {
    ## get every position found in the reads
    readSNPs_pos <- unlist(lapply(sampleReads, function(x) x[[2]])) ## 0-based
    readSNPs_per_grid <- table(readSNPs_pos)
    offending_grids <- as.integer(names(which(readSNPs_per_grid > downsampleToCov)))
    remove_stats <- c(0, length(sampleReads))
    if (length(offending_grids) > 0) {
        toRemove <- array(FALSE, length(sampleReads))
        which_reads <- is.na(match(readSNPs_pos, offending_grids)) == FALSE
        readList <- 1:length(sampleReads)
        reads_per_grid <- tapply(readList[which_reads], readSNPs_pos[which_reads], I) ## hopefully not slow
        for(bad_grid in offending_grids) {
            reads_in_grid_window <- reads_per_grid[[as.character(bad_grid)]]
            keep <- sample(reads_in_grid_window, downsampleToCov)
            remove <- setdiff(reads_in_grid_window, keep)
            toRemove[remove] <- TRUE
        }
        ## done!
        if (verbose & print_warning_to_screen) {
            print_message(paste0(
                "Warning - In gridding procedure ", 
                "downsample sample ", sampleNames[iBam], ": ",
                sum(toRemove), " of ", length(sampleReads),
                " reads removed "
            ))
        }
        remove_stats <- c(sum(toRemove), length(sampleReads))
        sampleReads <- sampleReads[toRemove == FALSE]
    }
    return(
        list(
            sampleReads = sampleReads,
            remove_stats = remove_stats
        )
    )
}



## so want to determine what are the grid points and the SNPs in each output "block"
determine_snp_and_grid_blocks_for_output <- function(
    grid,
    start_and_end_minus_buffer,
    outputSNPBlockSize = 1000
) {
    to_out <- array(NA, c(ceiling(length(grid) / outputSNPBlockSize) + 1, 4))
    ## yay! mixture of 0 and 1 based indexing
    colnames(to_out) <- c("snp_start_1_based", "snp_end_1_based", "grid_start_0_based", "grid_end_0_based")
    if (start_and_end_minus_buffer[1] == start_and_end_minus_buffer[2]) {
        ## only 1 SNP in output region!
        to_out[1, c("snp_start_1_based", "snp_end_1_based")] <- start_and_end_minus_buffer
        to_out[1, c("grid_start_0_based", "grid_end_0_based")] <- start_and_end_minus_buffer
        c2 <- 2
    } else {
        c <- 1
        c2 <- 1
        start <- start_and_end_minus_buffer[1] ## start SNP, 1-based
        prev_grid <- grid[start]        
        for(iSNP in (start_and_end_minus_buffer[1] + 1):start_and_end_minus_buffer[2]) {
            ## so iSNP is 1-based here
            c <- c + 1
            if (iSNP == (start_and_end_minus_buffer[2])) {
                to_out[c2, c("snp_start_1_based", "snp_end_1_based")] <- c(start, iSNP)
                c2 <- c2 + 1
            } else if (
                (outputSNPBlockSize < c) &
                (grid[iSNP - 1] < grid[iSNP]) &
                ((prev_grid + 1) < (grid[iSNP]))
            ) {
                to_out[c2, c("snp_start_1_based", "snp_end_1_based")] <- c(start, iSNP - 1)
                c <- 0
                c2 <- c2 + 1
                start <- iSNP
                prev_grid <- grid[start]
            }
        }
    }
    to_out <- to_out[1:(c2 - 1), , drop = FALSE]
    to_out[, "grid_start_0_based"] <- grid[to_out[, "snp_start_1_based"]]
    to_out[, "grid_end_0_based"] <- grid[to_out[, "snp_end_1_based"]]
    ## now, if last output block has 1 grid in it, merge into previous
    if ((2 <= nrow(to_out) )) {
        c2 <- nrow(to_out) + 1
        if ((to_out[(c2 - 1), "grid_end_0_based"] - to_out[(c2 - 1), "grid_start_0_based"]) == 0) {
            to_out[(c2 - 2), "snp_end_1_based"] <-  to_out[(c2 - 1), "snp_end_1_based"]
            to_out[(c2 - 2), "grid_end_0_based"] <-  to_out[(c2 - 1), "grid_end_0_based"]
            c2 <- c2 - 1
            to_out <- to_out[-nrow(to_out), , drop = FALSE]
        }
    }
    ## if last region too small (smaller than 2 SNPs or 10%), merge back in
    ## print some stats
    mess <- paste0(
        "Outputting will be done in ", nrow(to_out), " blocks with on average ",
        1 + round(mean(to_out[, "snp_end_1_based"] - to_out[, "snp_start_1_based"]), 1),
        " SNPs in them"
    )
    if (tail(grid, 1) != (length(grid) - 1)) {
        mess <- paste0(
            mess,
            " and ",
            1 + round(mean(to_out[, "grid_end_0_based"] - to_out[, "grid_start_0_based"], 1)),
            " grids"
        )
    }
    print_message(mess)
    return(to_out)
}



reorder_alphaBetaBlocks <- function(single_iteration_results, sampleRanges, N, K, nGrids, useTempdirWhileWriting) {
    if (useTempdirWhileWriting) {
        allAlphaBetaBlocks <- NULL
    } else {
        allAlphaBetaBlocks <- as.list(1:N)
        for(i_core in 1:length(single_iteration_results)) {
            sampleRange <- sampleRanges[[i_core]]
            who_to_run <- sampleRange[1]:sampleRange[2]
            abSmall <- single_iteration_results[[i_core]]$allAlphaBetaBlocks
            for(iiSample in 1:length(who_to_run)) {
                iSample <- who_to_run[[iiSample]]
                allAlphaBetaBlocks[[iSample]] <- abSmall[[iiSample]]
            }
        }
    }
    ## also want hapSum
    hapSum_t <- array(0, c(K, nGrids))
    ## sum and sum
    for(i in 1:length(sampleRanges)) {
        hapSum_t <- hapSum_t + single_iteration_results[[i]]$hapSum_t
    }
    list(
        allAlphaBetaBlocks = allAlphaBetaBlocks,
        hapSum_t = hapSum_t
    )
}
    
