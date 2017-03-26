#' @title Sequencing To Imputation Through Constructing Haplotypes
#' @param chr What chromosome to run. Should match BAM headers
#' @param posfile Where to find file with positions to run. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>
#' @param K How many founder / mosaic haplotypes to use
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure
#' @param outputdir What output directory to use
#' @param tempdir What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/
#' @param bamlist Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai
#' @param cramlist Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into STITCH
#' @param reference Path to reference fasta used for making cram files. Only required if cramlist is defined
#' @param genfile Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header
#' @param method How to run imputation - either diploid or pseudoHaploid, the former being the original method quadratic in K, the later being linear in K
#' @param outputInputInVCFFormat Whether to output the input in vcf format
#' @param downsampleToCov What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage
#' @param downsampleFraction Downsample BAMs by choosing a fraction of reads to retain. Must be value 0<downsampleFraction<1
#' @param readAware Whether to run the algorithm is read aware mode. If false, then reads are split into new reads, one per SNP per read
#' @param chrStart When loading from BAM, some start position, before SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param chrEnd When loading from BAM, some end position, after SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param maxDifferenceBetweenReads How much of a difference to allow the reads to make in the forward backward probability calculation. For example, if P(read | state 1)=1 and P(read | state 2)=1e-6, re-scale so that their ratio is this value. This helps prevent any individual read as having too much of an influence on state changes, helping prevent against influence by false positive SNPs
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
#' @param regenerateInput Whether to regenerate input files
#' @param originalRegionName If regenerateInput is FALSE (i.e. using existing data), this is the name of the original region name (chr.regionStart.regionEnd). This is necessary to load past variables
#' @param keepInterimFiles Whether to keep interim parameter estimates
#' @param keepTempDir Whether to keep files in temporary directory
#' @param environment Whether to use server or cluster multicore options
#' @param pseudoHaploidModel How to model read probabilities in pseudo diploid model (shouldn't be changed)
#' @param switchModelIteration Whether to switch from pseudoHaploid to diploid and at what iteration (NA for no switching)
#' @param generateInputOnly Whether to just generate input data then quit
#' @param restartIterations In pseudoHaploid method, which iterations to look for collapsed haplotype prnobabilities to resolve
#' @param refillIterations When to try and refill some of the less frequently used haplotypes
#' @param downsampleSamples What fraction of samples to retain. Useful for checking effect of N on imputation. Not meant for general use
#' @param downsampleSamplesKeepList When downsampling samples, specify a numeric list of samples to keep
#' @param subsetSNPsfile If input data has already been made for a region, then subset down to a new set of SNPs, as given by this file. Not meant for general use
#' @param useSoftClippedBases Whether to use (TRUE) or not use (FALSE) bases in soft clipped portions of reads
#' @param outputBlockSize How many samples to write out to disk at the same time when making temporary VCFs that are later pasted together at the end to make the final VCF. Smaller means lower RAM footprint, larger means faster write.
#' @param inputBundleBlockSize If NA, disable bundling of input files. If not NA, bundle together input files in sets of <= inputBundleBlockSize together
#' @param reference_haplotype_file Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)
#' @param reference_legend_file Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)
#' @param reference_sample_file Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)
#' @param reference_populations Vector with character populations to include from reference_sample_file e.g. CHB, CHS
#' @param reference_phred Phred scaled likelihood or an error of reference haplotype. Higher means more confidence in reference haplotype genotypes, lower means less confidence
#' @param reference_iterations When using reference haplotypes, how many iterations to use to train the starting data
#' @param vcf_output_name Override the default VCF output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple STITCH runs are processing on the same region then they may over-write each others inputs and outputs
#' @param initial_min_hapProb Initial lower bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param initial_max_hapProb Initial upper bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param regenerateInputWithDefaultValues If regenerateInput is FALSE and the original input data was made using regionStart, regionEnd and buffer as default values, set this equal to TRUE
#' @param plotHapSumDuringIterations Boolean TRUE/FALSE about whether to make a plot that shows the relative number of individuals using each ancestral haplotype in each iteration
#' @param save_sampleReadsInfo Experimental. Boolean TRUE/FALSE about whether to save additional information about the reads that were extracted
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
    reference = "",
    genfile = "",
    method = "diploid",
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
    environment = "server",
    pseudoHaploidModel = 9,
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
    inputBundleBlockSize = NA,
    reference_haplotype_file = "",
    reference_legend_file = "",
    reference_sample_file = "",
    reference_populations = NA,
    reference_phred = 20,
    reference_iterations = 10,
    vcf_output_name = NULL,
    initial_min_hapProb = 0.4,
    initial_max_hapProb = 0.6,
    regenerateInputWithDefaultValues = FALSE,
    plotHapSumDuringIterations = FALSE,
    save_sampleReadsInfo = FALSE
) {


    phasefile = "" ## do not export for now
    ##  #' @param phasefile Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header

    outputHaplotypeProbabilities <- FALSE ## do not export for now
    ## #' @param outputHaplotypeProbabilities Whether to output haplotype probabilities in files

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
    check_program_dependency("samtools")
    check_program_dependency("bgzip")



    ##
    ##
    ## validate parameters
    ##
    ##
    validate_chr(chr)
    validate_nGen(nGen)
    validate_nCores(nCores)
    validate_posfile(posfile)
    validate_K(K)
    validate_outputdir(outputdir)
    validate_tempdir(tempdir)
    validate_outputBlockSize(outputBlockSize)
    validate_downsampleFraction(downsampleFraction)
    validate_downsampleSamples(downsampleSamples)
    validate_plotHapSumDuringIterations(plotHapSumDuringIterations)

    ## more involved checks
    validate_regionStart_regionEnd_and_buffer(regionStart, regionEnd, buffer)
    validate_bamlist_and_cramlist_for_input_generation(regenerateInput, originalRegionName, bamlist, cramlist, regionStart, regionEnd, buffer, regenerateInputWithDefaultValues, reference)
    validate_vcf_output_name(vcf_output_name)
    validate_inputBundleBlockSize(inputBundleBlockSize, readAware)
    validate_reference_files(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations)
    validate_refillIterations(refillIterations, niterations)
    validate_shuffleHaplotypeIterations(shuffleHaplotypeIterations, niterations)
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
    tempdir <- out$tempdir
    inputdir <- out$inputdir


    ##
    ## print out date
    ##
    date <- date()
    print(paste(c("Program start - ",date),collapse=""))
    file <- file.path(outputdir, "RData", paste0("start.", regionName ,".RData"))
    if(generateInputOnly==FALSE)
        save(date, file = file)


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
    T <- out$T
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
        phase = phase,
        T = T
    )
    pos <- out$pos
    gen <- out$gen
    phase <- out$phase
    L <- out$L
    T <- out$T
    inRegionL <- out$inRegionL



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
        originalRegionName = originalRegionName
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
    generate_or_refactor_input(regenerateInput = regenerateInput, bundling_info = bundling_info, L = L, pos = pos, T = T, bam_files = bam_files, cram_files = cram_files, reference = reference, iSizeUpperLimit = iSizeUpperLimit, bqFilter = bqFilter, chr = chr, outputdir = outputdir, N = N, downsampleToCov = downsampleToCov, sampleNames = sampleNames, inputdir = inputdir, useSoftClippedBases = useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, generateInputOnly = generateInputOnly, environment = environment, nCores = nCores, save_sampleReadsInfo = save_sampleReadsInfo)





    ##
    ## if necessary, shrink BAMs, but only if regenerateInput = FALSE
    ##
    shrinkReads(N = N, nCores = nCores, sampleRange = sampleRange, originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL, environment = environment, regenerateInput = regenerateInput, inputBundleBlockSize = inputBundleBlockSize)



    ##
    ## if we want, re-intersect with new set of positions
    ##
    if(is.na(subsetSNPsfile)==FALSE & subsetSNPsfile!="NA")
    {
        print(paste("Subsetting SNPs from file ",subsetSNPsfile,sep=""))
        ## load in positions
        out=subsetSNPsFunction(N=N,subsetSNPsfile=subsetSNPsfile,regionName=regionName,tempdir=tempdir,L=L,environment=environment,nCores=nCores,outputdir=outputdir)
        print(paste("Back to main ",subsetSNPsfile,sep=""))
        keep=out$keep
        T=as.integer(sum(keep))
        L=L[keep]
        pos=pos[keep,]
        gen=gen[keep,]
        print(paste("Done subsetting section ",subsetSNPsfile,sep=""))
        save(T,L,keep,pos,gen,file=paste(outputdir,"debug/files.RData",sep=""))
        print(T)
    }



    ##
    ## downsample the number of samples
    ##
    if(downsampleSamples < 1)
    {
        out <- downsample_the_samples(
            downsampleSamplesKeepList = downsampleSamplesKeepList,
            downsampleSamples = downsampleSamples,
            highCovInLow = highCovInLow,
            nCores = nCores,
            blockSize = blockSize,
            N =N,
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
    if(downsampleFraction < 1)
        downsampleToFraction(N=N,nCores=nCores,downsampleFraction=downsampleFraction,regionName=regionName,tempdir=tempdir,environment=environment, bundling_info = bundling_info)




    ##
    ## build alleleCount
    ##
    alleleCount <- buildAlleleCount(T=T,N=N,nCores=nCores,regionName=regionName,tempdir=tempdir,environment=environment, bundling_info = bundling_info)



    ##
    ## toggle read aware
    ##
    if(readAware==FALSE) {
        print("Split reads")
        splitReadsCompletely(N=N,nCores=nCores,tempdir=tempdir,regionName=regionName,environment=environment)
        print("Done splitting reads")
    }



    ##
    ## if we're so inclined, output input in VCF format
    ##
    if (outputInputInVCFFormat == TRUE) {
        out=outputInputInVCFFunction(outputdir = outputdir, pos = pos, T = T, tempdir = tempdir, N = N, nCores = nCores, regionName = regionName, environment = environment, sampleNames = sampleNames, outputBlockSize = outputBlockSize, bundling_info = bundling_info)
        return(NULL)
    }


    if (method == "pseudoHaploid") {
        out <- initialize_readProbs(
            N = N,
            nCores = nCores,
            pseudoHaploidModel = pseudoHaploidModel,
            tempdir = tempdir,
            bundling_info = bundling_info,
            regionName = regionName,
            initial_min_hapProb = initial_min_hapProb,
            initial_max_hapProb = initial_max_hapProb
        )
    }



    ##
    ## initialize variables
    ##
    out <- initialize_parameters(reference_haplotype_file = reference_haplotype_file, reference_legend_file = reference_legend_file, reference_sample_file = reference_sample_file, reference_populations = reference_populations, reference_phred = reference_phred, reference_iterations = reference_iterations, T = T, K = K, L = L, pos = pos, inputBundleBlockSize = inputBundleBlockSize, nCores = nCores, regionName = regionName, alleleCount = alleleCount, startIterations = startIterations, windowSNPs = windowSNPs, expRate = expRate, nGen = nGen, tempdir = tempdir, outputdir = outputdir, pseudoHaploidModel = pseudoHaploidModel, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, minRate = minRate, maxRate = maxRate, regionStart = regionStart, regionEnd = regionEnd, buffer = buffer)
    sigmaCurrent <- out$sigmaCurrent
    eHapsCurrent <- out$eHapsCurrent
    alphaMatCurrent <- out$alphaMatCurrent
    hapSumCurrent <- out$hapSumCurrent
    priorCurrent <- out$priorCurrent
    reference_panel_SNPs <- out$reference_panel_SNPs



    ##
    ## run EM algorithm here
    ##
    print(paste("Begin EM - ", date()))
    date <- date()
    save(date, file = file.path(outputdir, "RData", paste0("startEM.", regionName, ".RData")))
    print(paste0("Number of SNPs - ", T))
    ## actually run now
    iteration <- max(0, startIterations)
    while(iteration < niterations) {
        iteration <- iteration + 1
        ##
        ## switch methods if appropriate
        ##
        if(is.na(switchModelIteration)==FALSE) {
            if(iteration==switchModelIteration) {
                method="diploid"
                print(paste("Switching from pseudoHaploid to diploid - iteration ",iteration," - ",date(),sep=""))
            }
        }
        ##
        ## fork out and get results
        ##
        ##mcparallel
        ##c <- mccollect(p)
        results <- completeSampleIteration(N=N,tempdir=tempdir,chr=chr,K=K,T=T,priorCurrent=priorCurrent,eHapsCurrent=eHapsCurrent,alphaMatCurrent=alphaMatCurrent,sigmaCurrent=sigmaCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,whatToReturn=whatToReturn,Jmax=Jmax,highCovInLow=highCovInLow,iteration=iteration,method=method,expRate=expRate,minRate=minRate,maxRate=maxRate,niterations=niterations,splitReadIterations=splitReadIterations,shuffleHaplotypeIterations=shuffleHaplotypeIterations,nCores=nCores,L=L,nGen=nGen,emissionThreshold=emissionThreshold,alphaMatThreshold=alphaMatThreshold,gen=gen,outputdir=outputdir,environment=environment,pseudoHaploidModel=pseudoHaploidModel,outputHaplotypeProbabilities=outputHaplotypeProbabilities,switchModelIteration=switchModelIteration,regionName=regionName,restartIterations=restartIterations,refillIterations=refillIterations,hapSumCurrent=hapSumCurrent,outputBlockSize=outputBlockSize, bundling_info = bundling_info, alleleCount = alleleCount, phase = phase, samples_with_phase = samples_with_phase)
        ## save previous results to disk
        eHapsFuture <- results$eHapsFuture
        alphaMatFuture <- results$alphaMatFuture
        sigmaFuture <- results$sigmaFuture
        priorFuture <- results$priorFuture
        hapSumFuture <- results$hapSum
        ## new - get hwe, info, estimatedAlleleFrequency
        hwe <- results$hwe
        info <- results$info
        estimatedAlleleFrequency <- results$estimatedAlleleFrequency
        ## save everything if needed
        if (keepInterimFiles)
            save(
                hapSumFuture,eHapsCurrent,alphaMatCurrent,sigmaCurrent,priorCurrent,eHapsFuture,alphaMatFuture,sigmaFuture,priorFuture,
                file = file.path(outputdir, "RData", paste0("interim.",regionName,".iteration",iteration,".RData"))
            )
        ## perform switchover here
        eHapsCurrent <- eHapsFuture
        alphaMatCurrent <- alphaMatFuture
        sigmaCurrent <- sigmaFuture
        priorCurrent <- priorFuture
        hapSumCurrent <- hapSumFuture
        hapSum <- hapSumCurrent
        if (plotHapSumDuringIterations == TRUE)
            plotHapSum(
                outname = file.path(outputdir, "plots", paste0("hapSum.",regionName,".iteration.",iteration,".jpeg")),
                L = L, K = K, hapSum = hapSum, T = T, N = N
            )
    }


    ##
    ## build final VCF from pieces
    ##
    write_vcf_after_EM(vcf_output_name = vcf_output_name, outputdir = outputdir, regionName = regionName, sampleNames = sampleNames, tempdir = tempdir, nCores = nCores, info = info, hwe = hwe, estimatedAlleleFrequency = estimatedAlleleFrequency, pos = pos, N = N, outputBlockSize = outputBlockSize, reference_panel_SNPs = reference_panel_SNPs, method = method)
    print(paste("Done EM and outputting - ",date()))


    ##
    ## get imputed dosages back in matrix form
    ##
    gen_imp <- NULL
    if(length(highCovInLow) > 0) {
        gen_imp <- get_high_coverage_final_imputed_genotypes(highCovInLow = highCovInLow, tempdir = tempdir, regionName = regionName, T = T)
    }


    ##
    ## now, if there was a buffer, remove the unwanted SNPs
    ##
    if(is.na(regionStart)==FALSE & is.na(regionEnd)==FALSE) {
        out <- remove_buffer_from_variables(L = L,  regionStart = regionStart, regionEnd = regionEnd, pos = pos, gen = gen, phase = phase, alleleCount =  alleleCount, highCovInLow = highCovInLow, gen_imp = gen_imp, alphaMatCurrent = alphaMatCurrent, sigmaCurrent = sigmaCurrent, eHapsCurrent = eHapsCurrent, hapSum = hapSum, hapSumCurrent = hapSumCurrent, info = info, hwe = hwe, estimatedAlleleFrequency = estimatedAlleleFrequency)
        pos <- out$pos
        gen <- out$gen
        phase <- out$phase
        gen_imp <- out$gen_imp
        alleleCount <- out$alleleCount
        L <- out$L
        T <- out$T
        alphaMatCurrent <- out$alphaMatCurrent
        sigmaCurrent <- out$sigmaCurrent
        eHapsCurrent <- out$eHapsCurrent
        hapSumCurrent <- out$hapSumCurrent
        hapSum <- out$hapSum
        priorCurrent <- out$priorCurrent
        info <- out$info
        hwe <- out$hwe
        estimatedAlleleFrequency <- out$estimatedAlleleFrequency
        ##
        ## shrink the VCF
        ##
        print(paste0("Remove buffer from VCF - ",date()))
        file <- get_output_vcf(vcf_output_name,  outputdir, regionName)
        file_temp <- paste0(file, ".temp")
        command <- paste0("gunzip -c ",file," | cat | awk '{if(substr($0, 1, 1)==",'"#"'," || ( $2>=", regionStart, " && $2<=", regionEnd, ")) {print $0}}' | bgzip -c > ",file_temp)
        system(command)
        system(paste0("mv ", file_temp, " ", file))
    }





    ##
    ## save various important objects to disk
    ##
    print(paste0("Save RData objects to disk - ",date()))
    save(info, file = file.path(outputdir, "RData", paste0("info.",regionName,".RData")))
    save(hwe, file = file.path(outputdir, "RData", paste0("hwe.",regionName,".RData")))
    passQC <- info > 0.4 & hwe > 1e-6 # one interpretation of passQC
    save(passQC, file = file.path(outputdir, "RData", paste0("passQC.",regionName,".RData")))
    save(hapSum, file = file.path(outputdir, "RData", paste0("hapSum.",regionName,".RData")))



    ##
    ## save important files from EM output (ie eHaps, niterations)
    ##
    save(
        eHapsCurrent, alphaMatCurrent,
        sigmaCurrent, priorCurrent,
        hapSumCurrent,
        file = file.path(
            outputdir, "RData", paste0("EMfinalestimates.", regionName, ".RData")
        )
    )
    save(
        alleleCount,
        estimatedAlleleFrequency,
        pos, gen, L,
        file = file.path(
            outputdir, "RData", paste0("EMnecessary.", regionName, ".RData")
        )
    )
    save(
        nCores, N, T, chr, K, niterations,
        file = file.path(
            outputdir, "RData", paste0("EMparameters.", regionName, ".RData")
        )
    )
    if(length(highCovInLow) > 0)
        save(
            gen_imp,
            file = file.path(
                outputdir, "RData", paste0("gen_imp.", regionName, ".RData")
            )
        )



    ##
    ## do some plotting
    ##
    print(paste0("Make plots - ", date()))
    plotMetricsForPostImputationQC(iSample=NULL,highCovList=NULL,gen=gen, gen_imp = gen_imp, alleleCount=alleleCount,chr=chr,L=L,estimatedAlleleFrequency=estimatedAlleleFrequency,info=info,outputdir=outputdir,hwe=hwe,regionName=regionName)
    ##
    ## plot hapProbs sum - should tell what states are being used
    ##
    plotHapSum(outname = file.path(outputdir, "plots", paste0("hapProbs.run.",regionName,".jpg")),L=L,K=K,hapSum=hapSum,T=T,N=N)
    ##
    ## plot estimated AF against real 0.1X pileups (get r2 as well)
    ##
    if(sum(passQC)>1)
        plotEstimatedAgainstReal(outputdir=outputdir,alleleCount=alleleCount,estimatedAlleleFrequency=estimatedAlleleFrequency,which=passQC,chr=chr,regionName=regionName)





    ##
    ## clean up and terminate
    ##
    print(paste0("Clean up and end - ",date()))
    warnings()
    ## remove things from tempdir
    if(keepTempDir==FALSE) {
        system(paste("rm ",tempdir,"sample.","*",".input.",regionName,".RData",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
        system(paste("rm ",tempdir,"sample*RData",sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
        system(paste("rmdir ",tempdir,sep=""),ignore.stdout=TRUE,ignore.stderr=TRUE)
    }

    date=date()
    print(paste(c("Program done - ",date),collapse=""))
    save(date, file = file.path(outputdir, "RData", paste0("end.",regionName,".RData")))

    return(NULL)

}





validate_chr <- function(chr) {
    if (length(chr) == 0)
        stop("Please specify chr, the chromosome to impute")
    if(chr=="")
        stop("Please specify chr, the chromosome to impute")
}


validate_nCores <- function(nCores) {
    if (is.numeric(nCores) == FALSE) {
        stop("nCores must be an integer")
    }
    if (nCores < 1)
        stop("nCores must be greater than or equal to 1")
    if (round(nCores) != nCores)
        stop("nCores must be an integer")
}

validate_nGen <- function(nGen) {
    nGen_error <- paste0(
        "Please specify nGen, the estimate of the number of generations ",
        "since founding. Can be approximated by 4 * N_e / K, where",
        "N_e ~= 20,000 for humans"
    )
    if (is.na(nGen))
        stop(nGen_error)
    if( nGen == "")
        stop(nGen_error)
    if (is.numeric(nGen) == FALSE)
        stop("nGen needs to be a numeric")
    if (nGen < 0)
        stop("nGen must be greater than 0")
}

validate_posfile <- function(posfile)
    if(posfile=="")
        stop("Please specify posfile, the file with sites to impute over")

validate_K <- function(K)
    if(K=="")
        stop("Pleasy specify K")

validate_outputdir <- function(outputdir)
    if(outputdir=="")
        stop("Please specify outputdir")

validate_tempdir <- function(tempdir) {
    if (is.na(tempdir) == FALSE) {
        if(tempdir=="")
            stop(
                paste0(
                    "Please specify tempdir. If in doubt, use tempdir = tempdir()"
                )
            )
    }
}

validate_outputBlockSize <- function(outputBlockSize) {
    if(outputBlockSize < 1)
        stop("outputBlockSize must be greater than 1")
    if(round(outputBlockSize) != outputBlockSize)
        stop("outputBlockSize must be an integer")
    return(NULL)
}


validate_downsampleFraction <- function(downsampleFraction) {
    if(downsampleFraction < 0 | downsampleFraction > 1)
        stop("downsampleFraction must be between 0 and 1")
}

validate_downsampleSamples <- function(downsampleSamples) {
    if (is.numeric(downsampleSamples) == FALSE)
        stop("downsampleSamples must be a integer > 0 and <=1")
    if (downsampleSamples <= 0 | downsampleSamples >1)
        stop("downsampleSamples must be a integer > 0 and <=1")
    return(NULL)
}

validate_plotHapSumDuringIterations <- function(plotHapSumDuringIterations)
    if (is.logical(plotHapSumDuringIterations) == FALSE | is.na(plotHapSumDuringIterations))
        stop("plotHapSumDuringIterations must be either TRUE or FALSE")


validate_regionStart_regionEnd_and_buffer <- function(regionStart, regionEnd, buffer) {
    w <- as.integer(is.na(regionStart))+
        as.integer(is.na(regionEnd))+
        as.integer(is.na(buffer))
    if ( w !=0 & w!=3)
        stop("either all three of regionStart, regionEnd and buffer are NA, or none of them are")
    if (is.na(regionStart)==FALSE) {
        if (regionStart<1) stop("regionStart must be >0")
        if (round(regionStart)!=regionStart) stop("regionEnd must be an integer")
    }
    if (is.na(regionEnd)==FALSE) {
        if (regionEnd<1) stop("regionEnd must be >0")
        if (round(regionEnd)!=regionEnd) stop("regionEnd must be an integer")
    }
    if (is.na(buffer)==FALSE) {
        if (buffer < 0)
            stop("buffer must be >=0")
        if (round(buffer)!=buffer) stop("buffer must be an integer")
    }
    if (is.na(regionStart)==FALSE & is.na(regionEnd)==FALSE) {
        if (regionStart>regionEnd) stop("regionStart must be before regionEnd")
    }
}

validate_bamlist_and_cramlist_for_input_generation <- function(regenerateInput, originalRegionName, bamlist, cramlist, regionStart, regionEnd, buffer, regenerateInputWithDefaultValues, reference = NULL) {
    if (regenerateInput == FALSE) {
        if (is.na(originalRegionName) == TRUE)
            stop("if regenerateInput is FALSE (i.e. using existing data), you must supply the original region name (must not be NA) to load the input properly. Also don't forget to supply the position file used to make the original input data")
        if (bamlist != "" || cramlist != "")
            stop("If not regenerating input, please do not supply bamlist")
        if ((is.na(regionStart) | is.na(regionEnd) | is.na(buffer)) & regenerateInputWithDefaultValues == FALSE)
            stop("If regenerateInput is FALSE, please supply original regionStart, regionEnd and buffer used to make the inputs. This is used to subset the pos file appropriately. If you used the default values of regionStart, regionEnd and buffer (NA), please set regenerateInputWithDefaultValues to TRUE")
    } else {
        if (bamlist == "" & cramlist == "")
            stop("If regenerateInput is TRUE, please supply either bamlist or cramlist")
    }
    if (cramlist != "") {
        if (bamlist != "")
            stop("Please specify only one or cramlist or bamlist")
        if (reference == "")
            stop("If you supply a list of CRAM files, you need to supply the reference they were aligned against")
    }
}


validate_inputBundleBlockSize <- function(inputBundleBlockSize, readAware) {
    if (is.na(inputBundleBlockSize) == FALSE) {
        if (readAware == FALSE)
            stop("readAware is FALSE is not currently supported when bundling inputs together")
        if (inputBundleBlockSize < 1)
            stop("inputBundleBlockSize must be an integer greater than 0")
        if (round(inputBundleBlockSize) != inputBundleBlockSize)
            stop("inputBundleBlockSize must be an integer")
    }
}




validate_reference_files <- function(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations) {
    if (reference_haplotype_file == "" &
        reference_legend_file != "")
        stop("If you specify reference_legend_file you must also specify reference_haplotype_file")
    if (reference_haplotype_file != "" &
        reference_legend_file == "")
        stop("If you specify reference_haplotype_file you must also specify reference_legend_file")
    if (reference_sample_file != "" & is.na(reference_populations[1]) == TRUE)
        stop("The reference sample file has been provided but you are not specifying any populations to retain. If you want to use all reference haplotypes, please omit this file")
    if (reference_sample_file != "" & file.exists(reference_sample_file) == FALSE)
        stop(paste0("Cannot find reference_sample_file:", reference_sample_file))
    if (reference_legend_file != "" & file.exists(reference_legend_file) == FALSE)
        stop(paste0("Cannot find reference_legend_file:", reference_legend_file))
    if (reference_haplotype_file != "" & file.exists(reference_haplotype_file) == FALSE)
        stop(paste0("Cannot find reference_haplotype_file:", reference_haplotype_file))
  if (reference_haplotype_file != "" & niterations == 1)
      stop ("You have selected niterations = 1 with reference haplotypes. Since the reference panel may be missing sites in the posfile, and since these sites will be initialized at random, please use at least two iterations to ensure genotypes at these sites will not be noise")
}


validate_refillIterations <- function(refillIterations, niterations) {
    if (is.na(refillIterations[1]) == FALSE) {
        x <- niterations - refillIterations
        if (sum (x < 5) > 0)
            stop(paste0("The parameter niterations is set to ", niterations, " and refillIterations is set to c(", paste(refillIterations, collapse = ","), "). refillIterations (trying to refill infrequently used ancestral haplotypes) has a small but beneficial improvement on imputation performance, but due to methods involve, will initially lead to a decrease in imputation performance for the next few EM iterations. Please set niterations > (max(refillIterations) + 4), or if you want to use a small number of iterations, it is recommended you disable refillIterations by setting this parameter equal to NA"))
    }
}

validate_shuffleHaplotypeIterations <- function(shuffleHaplotypeIterations, niterations) {
    if (is.na(shuffleHaplotypeIterations[1]) == FALSE) {
        x <- niterations - shuffleHaplotypeIterations
        if (sum (x < 5) > 0)
            stop(paste0("The parameter niterations is set to ", niterations, " and shuffleHaplotypeIterations is set to c(", paste(shuffleHaplotypeIterations, collapse = ","), "). shuffleHaplotypeIterations (trying to reshuffle the ancestral haplotypes to minimize switching among ancestral haplotypes in analyzed samples) has a small but beneficial improvement on imputation performance, but due to methods involve, will initially lead to a decrease in imputation performance for the next few EM iterations. Please set niterations > (max(shuffleHaplotypeIterations) + 4), or if you want to use a small number of iterations, it is recommended you disable shuffleHaplotypeIterations by setting this parameter equal to NA"))
    }
}


validate_hapProb <- function(initial_min_hapProb, initial_max_hapProb) {
    if (is.numeric(initial_min_hapProb) == FALSE)
        stop("initial_min_hapProb must be numeric")
    if (is.numeric(initial_max_hapProb) == FALSE)
        stop("initial_max_hapProb must be numeric")
    if (initial_min_hapProb < 0 | 1 < initial_min_hapProb)
        stop("initial_min_hapProb must be between 0 and 1")
    if (initial_max_hapProb < 0 | 1 < initial_max_hapProb)
        stop("initial_max_hapProb must be between 0 and 1")
}


## make sure the directories exist
initialize_directories <- function(
    tempdir,
    keepTempDir,
    outputdir
) {
    ## not ideal - overriding input parameter
    ## also not ideal - same name as function
    if (is.na(tempdir)) {
        tempdir <- tempdir()
    }
    tempdir <- paste0(tempdir, "/", paste(toupper(letters[sample(26,10,replace=TRUE)]),collapse=""),"/")
    if (keepTempDir == TRUE)
        print(paste0("tempdir=", tempdir))
    out <- sapply(
        c(tempdir,
          outputdir,
          file.path(outputdir, "RData"),
          file.path(outputdir, "debug"),
          file.path(outputdir, "plots"),
          file.path(outputdir, "input")
          ),
        function(dir)
        dir.create(dir, showWarnings = FALSE)
    )
    inputdir <- file.path(outputdir, "input")
    return(
        list(
            tempdir = tempdir,
            inputdir = inputdir
        )
    )
}


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



remove_buffer_from_variables <- function(
    L,
    regionStart,
    regionEnd,
    pos,
    gen = NULL,
    phase = NULL,
    alleleCount,
    highCovInLow = NULL,
    gen_imp = NULL,
    alphaMatCurrent = NULL,
    sigmaCurrent = NULL,
    eHapsCurrent = NULL,
    hapSum = NULL,
    hapSumCurrent = NULL,
    info = NULL,
    hwe = NULL,
    estimatedAlleleFrequency = NULL,
    phase_entropy = NULL,
    strandedness = NULL,
    eHapsUpdate_numer = NULL,
    eHapsUpdate_denom = NULL
) {
    print(paste0("Remove buffer region from output - ",date()))
    ## determine where the region is
    inRegion2 <- L >= regionStart & L <= regionEnd
    inRegion2L <- which(inRegion2)
    ## pos, gen, L
    pos <- pos[inRegion2, ]
    gen <- gen[inRegion2, ]
    phase <- phase[inRegion2, , , drop = FALSE]
    if ( length(highCovInLow) > 0 )
        gen_imp <- gen_imp[inRegion2, ]
    alleleCount <- alleleCount[inRegion2, ]
    L <- L[inRegion2]
    T <- as.integer(nrow(pos))
    ## parameters and hapSum
    x <- inRegion2L
    x <- x[-length(x)]
    alphaMatCurrent <- alphaMatCurrent[x,]
    sigmaCurrent <- sigmaCurrent[x]
    eHapsCurrent <- eHapsCurrent[inRegion2,]
    hapSumCurrent <- hapSumCurrent[inRegion2,]
    hapSum <- hapSum[inRegion2,]
    priorCurrent <- hapSum[1,]/sum(hapSum[1,])
    ##
    info <- info[inRegion2]
    hwe <- hwe[inRegion2]
    estimatedAlleleFrequency <- estimatedAlleleFrequency[inRegion2]
    ## phasing stuff
    phase_entropy <- phase_entropy[inRegion2]
    strandedness <- strandedness[inRegion2, , drop = FALSE]
    eHapsUpdate_numer <- eHapsUpdate_numer[inRegion2, , drop = FALSE]
    eHapsUpdate_denom <- eHapsUpdate_denom[inRegion2, , drop = FALSE]
    return(
        list(
            inRegion2 = inRegion2,
            pos = pos,
            gen = gen,
            phase = phase,
            gen_imp = gen_imp,
            alleleCount = alleleCount,
            L = L,
            T = T,
            alphaMatCurrent = alphaMatCurrent,
            sigmaCurrent = sigmaCurrent,
            eHapsCurrent = eHapsCurrent,
            hapSumCurrent = hapSumCurrent,
            hapSum = hapSum,
            priorCurrent = priorCurrent,
            info = info,
            hwe = hwe,
            estimatedAlleleFrequency = estimatedAlleleFrequency,
            phase_entropy = phase_entropy,
            strandedness = strandedness,
            eHapsUpdate_numer = eHapsUpdate_numer,
            eHapsUpdate_denom = eHapsUpdate_denom
        )
    )
}



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





## get sample names from disk for bam files or cram files
get_sample_names <- function(
    bamlist = "",
    cramlist = "",
    nCores = 1,
    outputdir = NULL,
    regionName = NULL,
    originalRegionName = NULL,
    save = TRUE,
    verbose = TRUE
) {

    bam_files <- NULL
    cram_files <- NULL
    if (bamlist != "" | cramlist != "")
    {
        if (cramlist != "") {
            files <- as.character(read.table(cramlist)[,1])
            cram_files <- files
            file_type = "CRAM"
        } else {
            files <- as.character(read.table(bamlist)[,1])
            bam_files <- files
            file_type = "BAM"
        }
        sampleNames <- get_sample_names_from_bam_or_cram_files(
            files,
            nCores,
            file_type,
            verbose
        )

        if(length(unique(sampleNames))!=length(sampleNames))
            stop("There are repeat sample names")

        ## save sample names
        if (save)
            save(
                sampleNames,
                file = file.path(
                    outputdir, "RData", paste0("sampleNames.", regionName, ".RData"))
            )

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
shrink_region <- function(
    regionStart,
    regionEnd,
    buffer,
    L,
    pos,
    gen,
    phase,
    T
) {
    if ( is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE ) {
        ## need to shrink pos, L, gen
        ## then, need to fix input
        inRegion <- L>=(regionStart-buffer) & L<=(regionEnd+buffer)
        inRegionL <- range((1:length(inRegion))[inRegion]) # this may be used later on
        L <- L[inRegion]
        pos <- pos[inRegion,]
        gen <- gen[inRegion,, drop = FALSE]
        phase <- phase[inRegion, , , drop = FALSE]
        T <- as.integer(nrow(pos))
    } else {
        inRegionL <- c(1, T)
    }
    return(
        list(
            pos = pos,
            gen = gen,
            phase = phase,
            L = L,
            T = T,
            inRegionL = inRegionL
        )
    )
}



get_and_validate_pos_gen_and_phase <- function(
    posfile,
    genfile = "",
    phasefile = "",
    chr,
    verbose = FALSE
) {
    if (verbose)
        print(paste0("Get and validate pos and gen - ", date()))
    pos <- get_and_validate_pos(posfile, chr)
    gen <- get_and_validate_gen(genfile)
    phase <- get_and_validate_phase(phasefile)
    if (is.null(gen) == FALSE)
        if(nrow(gen) != nrow(pos))
            stop(paste0(
                "posfile and genfile must have the same number of SNPs ",
                "but posfile has ", nrow(pos), " SNPs and genfile has ",
                nrow(gen), " SNPs (and 1 header row)"
            ))
    if (is.null(phase) == FALSE)
        if(dim(phase)[1] != nrow(pos))
            stop(paste0(
                "posfile and phasefile must have the same number of SNPs ",
                "but posfile has ", nrow(pos), " SNPs and phasefile has ",
                dim(phase)[1], " SNPs (and 1 header row)"
            ))
    L <- as.integer(as.character(pos[, "POS"]))
    T <- as.integer(length(L))
    if (verbose)
        print(paste0("Done get and validate pos and gen - ", date()))
    return(
        list(
            pos = pos,
            gen = gen,
            phase = phase,
            T = T,
            L = L
        )
    )
}

get_and_validate_gen <- function(genfile) {
    if (genfile == "")
        return(NULL)
    gen2 <- read.table(genfile, header=TRUE, sep="\t")
    gen <- array(0, c(nrow(gen2), ncol(gen2)))
    for(i_col in 1:ncol(gen)) {
        col <- gen2[, i_col]
        validate_gen_col(col, i_col)
        gen[, i_col] <- as.integer(as.character(gen2[, i_col]))
    }
    header <- as.character(unlist(read.table(genfile, sep="\t", nrows = 1)))
    validate_gen_header(header)
    colnames(gen) <- header
    return(gen)
}


validate_gen_header <- function(header) {
    ## only allowed are 0, 1, 2, NA
    m <- match(header, c(0, 1, 2, NA))
    if (sum(is.na(m) == FALSE) > 0) {
        m2 <- which.max(is.na(m) == FALSE)
        stop(paste0(
            "The header for the genfile is either invalid or missing. ",
            "The first invalid entry is in position ", m2, " and is entry ", header[m2]
        ))
    }
}

## can only be 0, 1, 2, NA
validate_gen_col <- function(col, i_col) {
    m <- match(col, c(0, 1, 2, NA))
    if (sum(is.na(m)) > 0) {
        m2 <- which.max(is.na(m))
        stop(paste0(
            "Ineligible entry for genfile column ", i_col, " at position ",
            m2, " with entry ", col[m2], ". Acceptable entries are 0, 1, 2 or NA"
        ))
    }
}

get_and_validate_pos <- function(posfile, chr) {
    pos <- read.table(posfile, header = FALSE, sep = "\t")
    colnames(pos) <- c("CHR", "POS", "REF", "ALT")
    validate_pos(pos, chr)
    return(pos)
}


## validate pos matrix shortly after loading
validate_pos <- function(
    pos,
    chr
    ) {
    if (sum(pos[,1] != chr) > 0) {
        m <- which.max(is.na(match(pos[,1], chr)))
        stop(paste0(
            "pos file needs to be unique to chromosome ",
            "and you supplied chr=", chr, " but ",
            "chromosome ", pos[m, 1], " was observed in row ", m
        ))
    }
    x <- suppressWarnings(is.na(as.integer(as.character(pos[, 2]))))
    if (sum(x) > 0) {
        stop(paste0(
            "pos file column 2 needs to be integer valued ",
            "between 1 and ", .Machine$integer.max, " but ",
            "in row ", which.max(x), " ", pos[which.max(x), 2],
            " was observed"
        ))
    }
    if (sum(diff(pos[, 2]) <= 0) > 0) {
        w <- which.max(diff(pos[, 2]) <= 0)
        stop(paste0(
            "pos file column 2 needs to be sorted on position ",
            "with increasing positions between rows ",
            "but row number ", w, " has position ", pos[w, 2], " ",
            "and row number ", w + 1, " has position ", pos[w + 1, 2]
        ))
    }
    for(col in 3:4) {
        x <- is.na(match(pos[, col], c("A", "C", "G", "T")))
        if (sum(x) > 0) {
            m <- which.max(x)
            y <- c("", "", "ref", "alt")[col]
            stop(paste0(
                "pos file ", y, " column entry ", pos[m, col],
                " in row ", m, " ",
                "contains is not one or A, C, G or T. STITCH is ",
                "only supported for bi-allelic SNPs"
            ))
        }
    }
    return(NULL)
}



validate_vcf_output_name <- function(
    vcf_output_name
) {
    if (is.null(vcf_output_name) == FALSE) {
        if (nchar(vcf_output_name) < 8)
            stop(paste0("vcf_output_name must have at least 8 characters and end with .vcf.gz, and you have supplied vcf_output_name:", vcf_output_name))
        if (substr(vcf_output_name, nchar(vcf_output_name) - 6, nchar(vcf_output_name)) != ".vcf.gz")
            stop(paste0("vcf_output_name must end with .vcf.gz, and you have supplied vcf_output_name:", vcf_output_name))
    }
    return(NULL)
}




validate_phase_header <- function(phasefile) {
    first_row <- as.character(unlist(read.table(phasefile,sep="\t",nrows=1)))
    ## check none are 0|0, 0|1, 1|0, 1|1
    m <- match(first_row, c("0|0", "0|1", "1|0", "1|1"))
    if (sum(is.na(m) == FALSE) > 0) {
        m2 <- which.max(is.na(m) == FALSE)
        stop(paste0(
            "The header for the phasefile is either invalid or missing. ",
            "The first invalid entry is in position ", m2, " and is entry ", first_row[m2]
        ))
    }
}

validate_phase_col <- function(col, i_samp) {
    x <- sapply(col, length)
    if (sum(x != 2) > 0) {
        m <- which.max(x)
        stop(paste0(
            "Unable to split column ", i_samp, " of phasefile at position ", m,
            " with entry ", col[m], " due to lack of field separator |"
        ))
    }
}


get_and_validate_phase <- function(
    phasefile
) {
    if (phasefile == "")
        return(NULL)
    phaseX <- read.table(phasefile, header = TRUE, stringsAsFactors = FALSE)
    validate_phase_header(phasefile)
    phase <- array(0, c(nrow(phaseX), ncol(phaseX), 2))
    for(i_samp in 1:ncol(phase)) {
        col <- strsplit(phaseX[, i_samp], "|", fixed = TRUE)
        validate_phase_col(col, i_samp)
        for(j in 1:2)
            phase[, i_samp, j] <- as.integer(sapply(col, function(x) x[j]))
    }
    if (length(colnames(phaseX)) == 1) {
        dimnames(phase)[[2]] <- list(colnames(phaseX))
    } else {
        dimnames(phase)[[2]] <- colnames(phaseX)
    }
    validate_phase(phase)
    return(phase)
}


validate_phase <- function(phase) {
    if (length(dim(phase)) != 3)
        stop("The phasefile does not have the right number of dimensions")
    if (sum(phase[, , 1] != 0 & phase[, , 1] != 1 & phase[, , 2] != 0 & phase[, , 2] != 1) > 0) {
        for(i_col in 1:dim(phase)[2]) {
            p <- phase[, i_col, ]
            i_row <- which.max(rowSums(p != 0 & p != 1) > 0)
            p2 <- phase[i_row, i_col, ]
            stop(paste0(
                "The phasefile contains entries other than 0 or 1. ",
                "One such entry is in column ", i_col, " and row ", i_row, " ",
                " with value ", paste(p2, collapse = "|")
            ))
        }
    }
}


compare_reference_haps_against_alleleCount <- function(
  alleleCount,
  reference_haps,
  outputdir,
  regionName
) {

    all_cor <- suppressWarnings(
        cor(alleleCount[, 3], rowSums(reference_haps), use = "pairwise.complete.obs")
    )
    low_maf_cor <- NA
    high_maf_cor <- NA

    w <- alleleCount[,3] > 0.05 & alleleCount[,3] < 0.95
    if (sum(w) > 1)
        high_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], rowSums(reference_haps[w, ]), use = "pairwise.complete.obs")
        )

    w <- alleleCount[,3] < 0.05 | alleleCount[,3] > 0.95
    if (sum(w) > 1)
        low_maf_cor <- suppressWarnings(
            cor(alleleCount[w, 3], rowSums(reference_haps[w, ]), use = "pairwise.complete.obs")
        )

    print(paste0("The following correlations are observed between allele frequencies estimated from the sequencing pileup and the reference haplotype counts:"))
    print(paste0(round(all_cor, 3), " for all SNPs"))
    print(paste0(round(high_maf_cor, 3), " for > 5% MAF SNPs"))
    print(paste0(round(low_maf_cor, 3), " for < 5% MAF SNPs"))

    N_haps <- ncol(reference_haps)

    out_plot <- file.path(outputdir, "plots", paste0("alleleFrequency_pileup_vs_reference_haplotypes.", regionName, ".png"))
    print(paste0("A plot of allele frequencies from sequencing pileup vs reference haplotype counts is at:", out_plot))
    png(out_plot, height = 500, width = 500)
    plot(alleleCount[, 3], rowSums(reference_haps) / N_haps, xlab = "Allele frequency from pileup", ylab = "Allele frequency from reference haplotypes")
    dev.off()

    return(NULL)

}





initialize_parameters <- function(
  reference_haplotype_file,
  reference_legend_file,
  reference_sample_file,
  reference_populations,
  reference_phred,
  reference_iterations,
  T,
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
  buffer
) {

    print(paste0("Begin parameter initialization, ", date()))
    a1 <- 0
    b1 <- 1
    eHapsCurrent <- matrix(a1+(b1-a1)*runif(T*K),ncol=K)
    if(startIterations>0) # but, make the rest of them uniform - get info later
        eHapsCurrent[(2*windowSNPs+1):T,]=0.25
    ## initialize using rate assuming 0.5 cM/Mb, T=100
    sigmaCurrent <- exp(-nGen*expRate/100/1000000*diff(L)) # interior is morgans between SNPs
    alphaMatCurrent <- matrix(1/K,nrow=(T-1),ncol=K)
    priorCurrent <- rep(1/K, K)
    hapSumCurrent <- array(1/K, c(T, K))
    reference_panel_SNPs <- array(FALSE, T)

    if (reference_haplotype_file != "") {
        print(paste0("Begin initializing paramters using reference haplotypes, ", date()))

        ## get reference haplotypes matched to posfile
        ## NA's where there are no match
        ## only for populations of interest
        reference_haps <- get_haplotypes_from_reference(
            reference_haplotype_file = reference_haplotype_file,
            reference_legend_file = reference_legend_file,
            reference_sample_file = reference_sample_file,
            reference_populations = reference_populations,
            pos = pos,
            tempdir = tempdir,
            regionName = regionName,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            chr = chr
        )
        reference_panel_SNPs <- is.na(reference_haps[, 1])

        compare_reference_haps_against_alleleCount(
            alleleCount = alleleCount,
            reference_haps = reference_haps,
            outputdir = outputdir,
            regionName = regionName
        )

        if (K > ncol(reference_haps)) {
            ## fill in rest with noise
            print(paste0("You have set K to be more than the number of reference haplotypes. The rest of the K ancestral haplotypes will be filled with noise to start"))
            w <- is.na(eHapsCurrent[, 1])
            eHapsCurrent[w, 1:ncol(reference_haps)] <- reference_haps[w, ]

        } else if (K == ncol(reference_haps)) {
            print(paste0("There are exactly as many reference haplotypes as K. Using these haplotypes directly as the initial estimate of the ancestral haplotypes"))
            ## shrink them from 0 -> e and 1 -> (1-e)
            e <- 0.001
            reference_haps[reference_haps == 0] <- e
            reference_haps[reference_haps == 1] <- (1 - e)
            eHapsCurrent <- reference_haps
        } else {

            N_haps <- ncol(reference_haps)
            reference_bundling_info <- get_bundling_position_information(
                N = N_haps,
                nCores = nCores,
                blockSize = inputBundleBlockSize
            )

            make_fake_sample_reads_from_haplotypes(
                reference_haps = reference_haps,
                nCores = nCores,
                reference_bundling_info = reference_bundling_info,
                tempdir = tempdir,
                N_haps = N_haps,
                reference_phred = reference_phred,
                regionName = regionName
            )

            ## chose some haps at random to fill in eHaps
            cols_to_replace <- sample(1:ncol(reference_haps), K)
            n1 <- nrow(reference_haps)
            n2 <- length(cols_to_replace)
            noise <- matrix(
                runif(n1 * n2, min = 0, max = 1),
                nrow = n1, ncol = n2
            )
            eHapsCurrent <- 0.99 * reference_haps[, cols_to_replace] + 0.01 * noise
            n1 <- sum(is.na(eHapsCurrent[, 1]))
            if (n1 > 0) {
                n2 <- ncol(eHapsCurrent)
                to_replace <- matrix(
                    runif(n1 * n2, min = 0.4, max = 0.6),
                    nrow = n1, ncol = n2
                )
                eHapsCurrent[is.na(eHapsCurrent[, 1]), ] <- to_replace
            }

            out <- run_EM_on_reference_sample_reads(reference_iterations = reference_iterations, sigmaCurrent = sigmaCurrent, eHapsCurrent = eHapsCurrent, alphaMatCurrent = alphaMatCurrent, hapSumCurrent = hapSumCurrent, priorCurrent = priorCurrent, N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, T = T, K = K, L = L, nGen = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, minRate = minRate, maxRate = maxRate, pseudoHaploidModel = pseudoHaploidModel, reference_phred = reference_phred)

            sigmaCurrent <- out$sigmaCurrent
            eHapsCurrent <- out$eHapsCurrent
            alphaMatCurrent <- out$alphaMatCurrent
            hapSumCurrent <- out$hapSumCurrent
            priorCurrent <- out$priorCurrent

            ## add some noise to NA sites - just in case
            ## there are too many of them nearby
            n1 <- sum(is.na(reference_haps[, 1]))
            if (n1 > 0) {
                n2 <- ncol(eHapsCurrent)
                to_replace <- matrix(
                    runif(n1 * n2, min = 0.4, max = 0.6),
                    nrow = n1, ncol = n2
                )
                eHapsCurrent[is.na(reference_haps[, 1]), ] <- to_replace
            }

        }
        print(paste0("Done initializing paramters using reference haplotypes, ", date()))
    }

    print(paste0("Done parameter initialization, ", date()))

    return(
        list(
            sigmaCurrent = sigmaCurrent,
            eHapsCurrent = eHapsCurrent,
            alphaMatCurrent = alphaMatCurrent,
            hapSumCurrent = hapSumCurrent,
            priorCurrent = priorCurrent,
            reference_panel_SNPs = reference_panel_SNPs
        )
    )

}











run_EM_on_reference_sample_reads <- function(
  reference_iterations,
  eHapsCurrent,
  sigmaCurrent,
  alphaMatCurrent,
  priorCurrent,
  hapSumCurrent,
  N_haps,
  nCores,
  reference_bundling_info,
  tempdir,
  regionName,
  T,
  K,
  L,
  nGen,
  emissionThreshold,
  alphaMatThreshold,
  expRate,
  minRate,
  maxRate,
  pseudoHaploidModel,
  reference_phred
) {

  print(paste0("Begin EM using reference haplotypes, ", date()))
  for(iteration in 1:reference_iterations) {
    print(paste0("Reference iteration = ", iteration, ", ", date()))
    out <- single_reference_iteration(N_haps = N_haps, nCores = nCores, reference_bundling_info = reference_bundling_info, tempdir = tempdir, regionName = regionName, T = T, K = K, L = L, nGen  = nGen, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, expRate = expRate, minRate = minRate, maxRate = maxRate, pseudoHaploidModel = pseudoHaploidModel, reference_phred = reference_phred, sigmaCurrent = sigmaCurrent, eHapsCurrent = eHapsCurrent, alphaMatCurrent = alphaMatCurrent, hapSumCurrent = hapSumCurrent, priorCurrent = priorCurrent)
    sigmaCurrent <- out$sigmaSum
    eHapsCurrent <- out$gammaSum
    alphaMatCurrent <- out$alphaMatSum
    hapSumCurrent <- out$hapSum
    priorCurrent <- out$priorSum
  }

  print(paste0("Done EM using reference haplotypes, ", date()))

  return(
    list(
      sigmaCurrent = sigmaCurrent,
      eHapsCurrent = eHapsCurrent,
      alphaMatCurrent = alphaMatCurrent,
      hapSumCurrent = hapSumCurrent,
      priorCurrent = priorCurrent
    )
  )

}





single_reference_iteration <- function(
  N_haps,
  nCores,
  reference_bundling_info,
  tempdir,
  regionName,
  T,
  K,
  L,
  nGen,
  emissionThreshold,
  alphaMatThreshold,
  expRate,
  minRate,
  maxRate,
  pseudoHaploidModel,
  reference_phred,
  sigmaCurrent,
  eHapsCurrent,
  alphaMatCurrent,
  hapSumCurrent,
  priorCurrent
) {

    sampleRanges <- getSampleRange(N_haps, nCores)
    pRgivenH1 <- as.integer(rep(1, T))
    pRgivenH2 <- as.integer(rep(0, T))
    transMatRate_t <- get_transMatRate(
        method = "pseudoHaploid",
        sigmaCurrent
    )
    eHapsCurrent_t <- t(eHapsCurrent)
    alphaMatCurrent_t <- t(alphaMatCurrent)

    out2 <- mclapply(
      sampleRanges,
      mc.cores = nCores,
      function(sampleRange) {

        priorSum <- array(0,K)
        alphaMatSum_t <- array(0, c(K, T - 1))
        gammaSum_t <- array(0, c(K, T, 2))
        hapSum_t <- array(0,c(K, T))

        bundledSampleReads <- NULL
        for(iSample in sampleRange[1]:sampleRange[2]) {

            out <- get_sampleReads_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = reference_bundling_info,
                bundledSampleReads = bundledSampleReads,
                what = "referenceSampleReads"
            )
            sampleReads <- out$sampleReads
            bundledSampleReads <- out$bundledSampleReads

            fbsoL <- forwardBackwardHaploid(
                sampleReads = sampleReads,
                nReads = as.integer(length(sampleReads)),
                Jmax = as.integer(1),
                pi = priorCurrent,
                pRgivenH1 = pRgivenH1,
                pRgivenH2 = pRgivenH2,
                pState = matrix(0, c(2, 2)), # ? deprecated ?
                eHaps_t = eHapsCurrent_t,
                alphaMat_t = alphaMatCurrent_t,
                transMatRate_t = transMatRate_t,
                maxDifferenceBetweenReads = 1 / (10 **(-(reference_phred / 10))),
                whatToReturn = as.integer(0),
                suppressOutput = as.integer(1),
                model = as.integer(pseudoHaploidModel)
            )

            which <- fbsoL$gammaUpdate_t[1, , 1] > 0
            gammaSum_t[, which, ] <- gammaSum_t[, which , ] +
                fbsoL$gammaUpdate_t[, which, ]
            alphaMatSum_t <- alphaMatSum_t +
                fbsoL$jUpdate_t
            priorSum <- priorSum+fbsoL$gammaUpdate_t[, 1, 2]


        }

        return(
            list(
                gammaSum_t = gammaSum_t,
                alphaMatSum_t = alphaMatSum_t,
                priorSum = priorSum,
                hapSum_t = hapSum_t
            )
        )
    }
    )

    error_check <- sapply(out2, class) == "try-error"
    if (sum(error_check) > 0) {
        print(out[[which(error_check)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }

    out <- calculate_updates(
        out2 = out2, x3 = sampleRanges, K = K, T = T, N = N_haps,
        nGen = nGen, expRate = expRate, minRate = minRate, maxRate = maxRate,
        emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, L = L
    )
    sigmaSum <- out$sigmaSum
    priorSum <- out$priorSum
    alphaMatSum <- t(out$alphaMatSum_t) ## TRANSPOSE-CLEAN
    gammaSum <- t(out$gammaSum_t) ## TRANSPOSE-CLEAN
    hapSum <- t(out$hapSum_t) ## TRANSPOSE-CLEAN

    return(
        list(
            gammaSum = gammaSum,
            sigmaSum = sigmaSum,
            alphaMatSum = alphaMatSum,
            priorSum = priorSum,
            hapSum = hapSum
        )
    )

}




make_fake_sample_reads_from_haplotypes <- function(
  reference_haps,
  nCores,
  reference_bundling_info,
  tempdir,
  N_haps,
  reference_phred,
  regionName
) {

  print(paste0("Begin convert reference haplotypes to internal input format, ", date()))

  sampleRanges <- getSampleRange(N_haps, nCores)
  non_NA_cols <- which(is.na(reference_haps[ , 1]) == FALSE)

  out <- mclapply(
    sampleRanges,
    mc.cores = nCores,
    function(sampleRange) {
    for(iSample in sampleRange[1]:sampleRange[2]) {
      sampleReads <- lapply(
        non_NA_cols,
        function(x) {
          return(list(
            0, x - 1,
            reference_phred * (2 * reference_haps[x, iSample] - 1),
            x - 1
          ))
        }
      )
      save(
        sampleReads,
        file = file_referenceSampleReads(tempdir, iSample, regionName)
      )
      if (length(reference_bundling_info) > 0) {
        last_in_bundle <- reference_bundling_info$matrix[iSample, "last_in_bundle"]
        if (last_in_bundle == 1) {
          bundle_inputs_after_generation(
            bundling_info = reference_bundling_info,
            iBam = iSample,
            dir = tempdir,
            regionName = regionName,
            what = "referenceSampleReads"
          )
        }
      }
    }
    return(NULL)
  })

  error_check <- sapply(out, class) == "try-error"
  if (sum(error_check) > 0) {
    print(out[[which(error_check)[1]]])
    stop("There has been an error generating the input. Please see error message above")
  }

  print(paste0("End convert reference haplotypes to internal input format, ", date()))

}



validate_reference_sample_file <- function(samples, reference_sample_file) {
    for(col in c("POP", "SEX")) {
        if (is.na(match(col, colnames(samples))))
            stop(
                paste0(
                    "Cannot find column ", col, " in reference_sample_file:",
                    reference_sample_file
                )
            )
    }
    x <- samples[, "SEX"] != "male" & samples[, "SEX"] != "female"
    if (sum(x) > 0)
        stop(
            "Reference samples must have SEX male or female ",
            "and in row ", which.max(x), " there is entry ",
            samples[which.max(x), "SEX"], " from file:", reference_sample_file
        )
}



get_reference_colClasses <- function(
    reference_sample_file,
    reference_populations,
    chr
) {
    colClasses <- NA
    if (reference_sample_file != "") {
        reference_samples <- read.table(reference_sample_file, header = TRUE)
        validate_reference_sample_file(reference_samples, reference_sample_file)
        keep_samples <- which(
            is.na(
                match(
                    rep(as.character(reference_samples[, "POP"]), each = 2),
                    reference_populations
                )
            ) == FALSE
        )
        print(paste0("Keeping ", length(keep_samples) / 2, " out of ", nrow(reference_samples), " reference samples from populations:", paste0(reference_populations, collapse = ",")))
        colClasses <- rep("NULL", nrow(reference_samples) * 2)
        colClasses[keep_samples] <- "integer"
    }
    return(colClasses = colClasses)
}


print_and_validate_reference_snp_stats <- function(
    pos_snps,
    legend_snps,
    both_snps
) {
    ## print out stats
    print(paste0("In the region to be imputed plus buffer there are the following number of SNPs:"))
    print(paste0(length(pos_snps), " from the posfile"))
    print(paste0(length(legend_snps), " from the reference haplotypes"))
    print(paste0(length(both_snps), " in the intersection of the posfile and the reference haplotypes"))
    if (length(both_snps) == 0)
        stop("There are 0 SNPs that intersect between the posfile and the reference haplotypes. Please troubleshoot to see if there is an error, or re-run without reference haplotypes")
}




load_reference_legend <- function(
    legend_header,
    reference_legend_file,
    regionStart,
    regionEnd,
    buffer
) {
    ## extract reference
    r <- match("position", legend_header)
    if (is.na(regionStart) == FALSE) {
        extract_condition <- paste0(
            "$", r, " >= ", regionStart - buffer, " && ",
            "$", r, "<=", regionEnd + buffer
        )
    } else {
        extract_condition <- "NR>1" ## do not load the header
    }
    legend_and_range <- system(
        paste0(
            "gunzip -c ", reference_legend_file, " | ", "awk '{if(", extract_condition,  ") {print NR\" \"$0}}'"
        ), intern = TRUE
    )
    a <- strsplit(legend_and_range, " ")
    col <- as.integer(sapply(a, function(x) x[1]))
    legend <- t(sapply(a, function(x) x[-1]))
    colnames(legend) <- legend_header
    legend_position <- as.integer(legend[, "position"])
    position_in_haps_file <- as.integer(sapply(a, function(x) x[1])) - 1
    return(
        list(
            legend = legend,
            legend_position = legend_position,
            position_in_haps_file = position_in_haps_file
        )
    )
}


## extract the relevant rows from the reference haplotypes file
## do this using awk to avoid loading way too many sites
extract_validate_and_load_haplotypes <- function(
    legend,
    pos,
    reference_haplotype_file,
    position_in_haps_file,
    regionName,
    tempdir,
    colClasses
) {

    print(paste0("Extract reference haplotypes, ",date()))
    temp_haps_file <- paste0(tempdir, ".haps.", regionName, ".txt.gz")
    temp_extract_file <- paste0(tempdir, ".extract.", regionName, ".txt")

    pos_snps <- paste(pos[,2], pos[,3], pos[,4], sep = "-")
    legend_snps <- paste(
        legend[, "position"], legend[,"a0"], legend[, "a1"],
        sep = "-"
    )

    both_snps <- intersect(legend_snps, pos_snps)

    print_and_validate_reference_snp_stats(pos_snps, legend_snps, both_snps)

    ## only extract those we need
    ## note - legend_snps validated - so no duplicate legend SNPs
    hap_snps_to_extract <- is.na(match(legend_snps, both_snps)) == FALSE
    hap_snps_position_in_pos <- is.na(match(pos_snps, both_snps)) == FALSE

    ## extract rows from haplotypes file
    cat(
        position_in_haps_file[hap_snps_to_extract],
        file = temp_extract_file, sep = "\n"
    )

    awk_command <- paste0(
        "awk '{if (NR==FNR) {a[$1]=$1; next} ",
        "if (FNR in a) {print $0}}' "
    )
    out <- system(
        paste0(
            "gunzip -c ", reference_haplotype_file, " | ",
            awk_command, temp_extract_file, " - | ",
            "gzip > ", temp_haps_file
        )
    )

    print(paste0("Load reference haplotypes, ",date()))
    haps <- read.table(
        temp_haps_file,
        colClasses = colClasses,
        sep = " ",
        na.strings = "-"
    )
    system(paste0("rm ", temp_haps_file))

    haps <- as.matrix(haps)

    haps <- remove_NA_columns_from_haps(haps)
    validate_haps(haps, both_snps)

    haps_both <- array(NA, c(nrow(pos), ncol(haps)))
    haps_both[hap_snps_position_in_pos, ] <- haps

    return(haps_both)
}


remove_NA_columns_from_haps <- function(haps) {
    na_sum <- colSums(is.na(haps))
    x <- na_sum != 0 & na_sum != nrow(haps)
    if (sum(x) > 0)
        stop(
            "The missing haplotype entry '-' for a reference male sample on chromosome X does not make up the entirety of a sample haplotype for at least one sample for the region of interest"
        )
    if (sum(na_sum == 0) == 0)
        stop("There are no viable haplotypes from the reference samples for the region of interest")
    if (sum(na_sum == nrow(na_sum)) > 0)
        print(paste0("Removing ", X, " out of ", Y, " male haplotypes from the reference haplotypes", date()))
    haps <- haps[, na_sum == 0, drop = FALSE]
    return(haps)
}

## central function controlling the loading of reference haplotypes
get_haplotypes_from_reference <- function(
  reference_haplotype_file,
  reference_legend_file,
  reference_sample_file,
  reference_populations,
  pos,
  tempdir = tempdir(),
  regionName = "test",
  regionStart = NA,
  regionEnd = NA,
  buffer = NA,
  chr
) {

    print(paste0("Begin get haplotypes from reference, ", date()))

    print(paste0("Load and validate reference legend header, ",date()))
    legend_header <- as.character(unlist(read.table(reference_legend_file, nrow = 1, sep = " ")))
    validate_legend_header(legend_header)

    print(paste0("Load and validate reference legend, ",date()))
    colClasses <- get_reference_colClasses(
        reference_sample_file,
        reference_populations,
        chr
    )
    out <- load_reference_legend(
        legend_header,
        reference_legend_file,
        regionStart,
        regionEnd,
        buffer
    )
    legend <- out$legend
    legend_position <- out$legend_position
    position_in_haps_file <- out$position_in_haps_file

    validate_reference_legend(legend)

    haps <- extract_validate_and_load_haplotypes(
        legend,
        pos,
        reference_haplotype_file,
        position_in_haps_file,
        regionName,
        tempdir,
        colClasses
    )

    print(paste0("Succesfully extracted ", ncol(haps), " haplotypes from reference data, ", date()))
    print(paste0("End get haplotypes from reference, ", date()))

    return(haps)

}



validate_reference_legend <- function(
    legend
) {
    if (nrow(legend) == 0) {
        stop(paste0("After extraction, the legend file had no rows. Please check the intersection between the regionStart, regionEnd, buffer and the reference_legend_file:", reference_legend_file))
    }
    legend_snps <- paste(
        legend[, "position"], legend[,"a0"], legend[, "a1"],
        sep = "-"
    )
    t <- table(legend_snps)
    if (sum(t > 1) > 0) {
        ## argh R - get character not factor
        m <- match(names(t[t>1])[1], legend_snps)
        example <- sapply(
            legend[m, c("position", "a0", "a1")],
            as.character
        )
        stop(
            paste0(
                "There are ", sum(t > 1), " duplicate row ids. ",
                "One such example is ",
                paste0(example, collapse = " ")
            )
        )
    }
    return(NULL)
}


validate_haps <- function(
    haps,
    both_snps
) {
    if (nrow(haps) != length(both_snps))
        stop ("Something has gone wrong during processing the reference haplotype and legend files. Please double check the legend file and haplotype file have the same number of entries, other than headers, and there are SNPs in the target region in the reference haplotype file")
    if (sum(is.na(haps)) > 0)
        stop ("The reference_haplotype_file contains entries other than 0 or 1")
    if (sum(haps != 0 & haps != 1) > 0)
        stop ("The reference_haplotype_file contains entries other than 0 or 1")
    if (ncol(haps) <= 1)
        stop ("The reference_haplotype_file contains less than 2 haplotypes. Please use more at least 2 reference haplotypes. If you think this message is an error, please check the reference haplotype file is space separated")
    return(NULL)
}


validate_legend_header <- function(
    legend_header
) {
    if (length(legend_header) < 3)
        stop("The reference_legend_file has fewer than 3 columns. Perhaps it is not a space separated file?")
    if (is.na(match("position", legend_header)))
        stop ("Cannot find 'position' column in reference_legend_file")
    if (is.na(match("a0", legend_header)))
        stop ("Cannot find reference allele 'a0' column in reference_legend_file")
    if (is.na(match("a1", legend_header)))
        stop ("Cannot find alternate allele 'a1' column in reference_legend_file")
    return(NULL)
}




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

  print(paste0("Begin downsample samples to ",100 * downsampleSamples , "% of the original numbers, ", date()))
  keep <- array(FALSE,N)
  keelList <- NULL
  # keep - either those on a list, or the high coverage
  if (is.na(downsampleSamplesKeepList) == FALSE) {
    if (file.exists(downsampleSamplesKeepList)) {
      keepList <- read.table(downsampleSamplesKeepList)
      keep[unlist(keepList)]=TRUE
    }
  }
  if (length(highCovInLow) > 0) {
    keep[highCovInLow] <- TRUE
  }
  # sample remaining samples to keep
  which <- sample(
    1:sum(keep == FALSE),
    N * downsampleSamples - sum(keep)
  )
  keep[keep == FALSE][which] <- TRUE

  new_N <- sum(keep)
  keepL <- which(keep)
  new_sampleNames <- sampleNames[keep]

  if (length(highCovInLow > 0))
    new_highCovInLow <- match(sampleNames[highCovInLow], new_sampleNames)


  new_bundling_info <- get_bundling_position_information(
    N = new_N,
    nCores = nCores,
    blockSize = inputBundleBlockSize
  )

  sampleRanges <- getSampleRange(new_N, nCores)

  f <- function(sampleRange, keepL, tempdir, regionName, bundling_info, new_bundling_info) {
    bundledSampleReads <- NULL
    for(iSample in sampleRange[1]:sampleRange[2])
    {

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

      # hack - if inputBundleBlockSize = NA, then would overwrite
      # so save to slightly different name, then push
      regionNameLocal <- regionName
      if (length(bundling_info) == 0)
        regionNameLocal <- paste0(regionName, ".TEMP.")

      save(
        sampleReads,
        file = file_sampleReads(tempdir, newSample, regionNameLocal)
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


  if(environment=="server")
  {
    out2=mclapply(sampleRanges, mc.cores=nCores,FUN=f, keepL= keepL, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info, new_bundling_info = new_bundling_info)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, sampleRanges, fun=f, keepL= keepL, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info, new_bundling_info = new_bundling_info)
    stopCluster(cl)
  }
  error_check <- sapply(out2, class) == "try-error"
  if (sum(error_check) > 0) {
    print(out2[[which(error_check)[1]]])
    stop("There has been an error downsampling the samples. Please see error message above")
  }

  # continuation of hack from above
  # do not want to over-write files
  if (length(bundling_info) == 0) {
    # remove original ones, push new ones
    system(paste0("rm ", file_sampleReads(tempdir, "*", regionName)))
    regionNameLocal <- paste0(regionName, ".TEMP.")
    for(iSample in 1:new_N) {
      system(paste0("mv ",
        file_sampleReads(tempdir, iSample, regionNameLocal),
        " ",
        file_sampleReads(tempdir, iSample, regionName)
      ))
    }
  }


  save(
    sampleNames,
    keepL,
    file = paste0(
      outputdir, "RData/keepL.", regionName, ".RData"
    )
  )

  print(paste0("End downsample samples, ", date()))

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
generate_or_refactor_input <- function(
  regenerateInput,
  bundling_info,
  L,
  pos,
  T,
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
        out <- generate_input(bundling_info = bundling_info, L = L, pos = pos, T = T, bam_files = bam_files, cram_files = cram_files, reference = reference, iSizeUpperLimit = iSizeUpperLimit, bqFilter = bqFilter, chr = chr, N = N, downsampleToCov = downsampleToCov, sampleNames = sampleNames, inputdir = inputdir, useSoftClippedBases = useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, environment = environment, nCores = nCores, save_sampleReadsInfo = save_sampleReadsInfo)
        if(generateInputOnly==TRUE) {
            save(pos,file=paste(outputdir,"RData/pos.",regionName,".RData",sep=""))
            print(paste0("Done generating input - ",date()))
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
                environment = environment
            )
        }
    }
}





get_rebundled_files <- function(inputdir, regionName) {
  files <- system(paste0("ls ", file_bundledSampleReads(inputdir, "*", "*", regionName)), intern = TRUE)
  # get start and end of bundledSampleReads
  a <- unlist(strsplit(files, paste0(".", regionName, ".RData")))
  x <- strsplit(a, paste0( "bundledSamples."))
  if (length(grep("bundledSamples", outputdir)) > 0)
    stop ("Please choose an outputdir name that does not include the term bundledSamples")
  ranges <-  sapply(x, function(x) x[2])
  ranges <- t(sapply(strsplit(ranges, "-"), function(x) as.integer(x)))
  files <- files[order(ranges[, 1])]
  ranges <- ranges[order(ranges[, 1]), ]
  # also check - range is exact and fully spans
  return(list(files = files, ranges = ranges))
}



# if the files already exist, good to go
# otherwise, loop through files using 1 core
# when doesn't exist, load next bunch
# every once in a while, write to disk
#
rebundle_input <- function(
  inputdir,
  regionName,
  bundling_info,
  N,
  tempdir,
  nCores,
  environment
) {

  print(paste0("Rebundle inputs, ", date()))
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
    print(paste0("The previously bundled files are the right size. No need to rebundle"))
    print(paste0("Done rebundling inputs, ", date()))
    return(NULL)
  }

  rebundle_input
  newBundle <- NULL
  bundledSampleReads <- NULL


  # multi-core re-bundling
  sampleRanges <- getSampleRange(N, nCores)

  f <- function(sampleRange, ranges, files, tempdir, regionName, bundling_info) {
    # start with the first file
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
      save(sampleReads, file = file_sampleReads(inputdir, iSample, regionName))
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


  if(environment=="server")
  {
    out2=mclapply(sampleRanges, mc.cores=nCores,FUN=f, ranges = ranges, files  = files, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, sampleRanges, fun=f, ranges = ranges, files  = files, tempdir = tempdir, regionName = regionName, bundling_info = bundling_info)
    stopCluster(cl)
  }
  error_check <- sapply(out2, class) == "try-error"
  if (sum(error_check) > 0) {
    print(out2[[which(error_check)[1]]])
    stop("There has been an error rebundling the input. Please see error message above")
  }

  print(paste0("Done rebundling inputs, ", date()))
  return(NULL)
}



## so were definitely generating input here
## decide whether it is on a server or a cluster
## that's it!
generate_input <- function(
  bundling_info,
  L,
  pos,
  T,
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
    print(paste0("Generate inputs - ", date()))
    chrLength <- NA
    if (is.na(chrStart)) {
        print(paste0("Get chromosome length - ", date()))
        chrLength <- get_chromosome_length(iBam = 1, bam_files, cram_files, chr)
    }
    sampleRanges <- getSampleRange(N = N, nCores = nCores)
    if(environment == "server") {
        out <- mclapply(1:length(sampleRanges), mc.cores=nCores, FUN=loadBamAndConvert_across_a_range,sampleRanges = sampleRanges,bundling_info=bundling_info,L=L,pos=pos,T=T,bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
    }
    if(environment=="cluster") {
        cl = makeCluster(length(sampleRanges), type = "FORK")
        out = parLapply(cl, 1:length(sampleRanges), mc.cores=nCores, fun=loadBamAndConvert_across_a_range,sampleRange=sampleRange,bundling_info=bundling_info, L=L,pos=pos,T=T,bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
        stopCluster(cl)
    }
    if (length(unlist(out)) > 0)  {
        print(out[[which(sapply(out, length) > 0)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }
    print(paste0("Done geneating inputs - ", date()))
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
  T,
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
  #print(paste0("Convert inputs for iCore=", iCore, ", ", date()))
  sampleRange <- sampleRanges[[iCore]]
  for(iBam in sampleRange[1]:sampleRange[2]) {
    loadBamAndConvert(iBam = iBam, L=L,pos=pos,T=T,bam_files = bam_files,cram_files = cram_files,reference=reference,iSizeUpperLimit=iSizeUpperLimit,bqFilter=bqFilter,chr=chr,N=N,downsampleToCov=downsampleToCov,sampleNames=sampleNames,inputdir=inputdir,useSoftClippedBases=useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, chrLength = chrLength, save_sampleReadsInfo = save_sampleReadsInfo)
    # if a bundle iteration, scoop up samples, dump to disk
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
  #print(paste0("Done converting inputs for iCore=", iCore, ", ", date()))
  return(NULL)
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
  #print(paste0("Bundle what=", what, ", iBam=", iBam, ", iCore=", iC, ", iBundle=", iB, ", samples=", s, "-", e, ", ", date()))
  # now, get these samples and write to disk
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
        return(list(pRgivenH1 = pRgivenH1, pRgivenH2 = pRgivenH2, srp = srp, delta = delta))
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
      file = file_bundledSampleReads(dir, s, e, regionName)
    )
  } else if (what == "sampleProbs") {
    bundledSampleProbs <- bundledSampleObject
    save(
      bundledSampleProbs,
      file = file_bundledSampleProbs(dir, s, e, regionName)
    )
  } else if (what == "referenceSampleReads") {
    bundledSampleReads <- bundledSampleObject
    save(
      bundledSampleReads,
      file = file_bundledReferenceSampleReads(dir, s, e, regionName)
    )
  }
  if (remove_files == TRUE) {
    out <- lapply(samples_in_core_bundle, function(iBam) {
      system(paste0("rm ", file_samples(dir, iBam, regionName)))
  })
  }
  return(NULL)
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
  return(list(
    sampleReads = sampleReads,
    bundledSampleReads = bundledSampleReads
  ))
}




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
    delta <- sampleProbs$delta
  }
  return(list(
    pRgivenH1 = pRgivenH1,
    pRgivenH2 = pRgivenH2,
    srp = srp,
    delta = delta,
    bundledSampleProbs = bundledSampleProbs
  ))
}



# for an iteration, find intervals of where there are infrequently
# used haplotypes
getWhereToRefill=function(hapSum,T,hapFreqMin)
{
  # 0.5% frequency for a hap
  hapFreqMin=0.005
  # find runs of a certain length less than this value
  whereToRefill=lapply(1:K,hapSum=hapSum,T=T,hapFreqMin=hapFreqMin,function(k,hapFreqMin,T,hapSum) {
    x=hapSum[,k]
    # those less than the minimum
    y=as.integer(x<hapFreqMin*N)
    # find break points
    z=(y[-1]-y[-length(y)])!=0
    start=c(1,(1:(T-1))[z]+1)
    end=c((1:(T-1))[z],T)
    m=matrix(cbind(start,end,y[start]),ncol=3)
    colnames(m)=c("start","end","insufficient")
    # only keep where we need to redo
    m=matrix(m[m[,3]==1,]    ,ncol=3)
    # now - remove some of them - keep if >100
    if(nrow(m)==1)      return(m)
    #
    # remove those smaller than 100
    #
    m=matrix(m[(m[,2]-m[,1])>100,],ncol=3)
    return(m)
  })
  return(whereToRefill)
}




get_RG_lines_from_SeqLib <- function(file) {
    header_string <- get_header_using_SeqLib(file)
    x <- strsplit(header_string, "\n")[[1]]
    return(x[substr(x, 1, 3) == "@RG"])
}

convert_sam_header_rg_tags_to_sample_name <- function(header) {
    rg_spots <- lapply(header, strsplit, split = "\t")
    if (length(rg_spots) == 0)
        stop(paste0("There is no @RG tag (with sample name) for:", file))
    sm <- sapply(rg_spots, function(x) {
        rg <- x[[1]]
        sm <- substr(rg[substr(rg, 1, 3) == "SM:"], 4, 1000)
        return(sm)
    })
    if (length(unique(sm)) > 1)
        stop(paste0("There is more than one sample name in the header for:", file))
    return(sm[1])
}

get_sample_name_from_bam_file_using_external_samtools <- function(file) {
    header <- system(paste0("samtools view -H ", file, " | grep ^@RG"), intern = TRUE)
    return(convert_sam_header_rg_tags_to_sample_name(header))
}

get_sample_name_from_bam_file_using_SeqLib <- function(file) {
    header <- get_RG_lines_from_SeqLib(file)
    return(convert_sam_header_rg_tags_to_sample_name(header))
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
        print(paste0("Get ", file_type, " sample names - ", date()))

    sampleNames <- mclapply(
        files,
        mc.cores = nCores,
        get_sample_name_from_bam_file_using_SeqLib
    )

    sampleNames <- as.character(unlist(sampleNames))

    if (verbose)
        print(paste0("Done getting ", file_type, " sample names - ", date()))

    return(sampleNames)
}



## in a simple refill, look every 100 SNP. in each one with an average below a threshold, sample from other haplotypes to re-copy
refillSimple <- function(hapSum,T,K,eHapsCurrent,N) {
    hapFreqMin <- 0.005 # percentage
    nSNP <- 100
    if (T < nSNP) {
        r <- list(c(1, T))
    } else {
        r <- lapply(
            1:floor(T/nSNP),
            function(x) c(1 + nSNP * (x-1), nSNP * x)
        )
    }
    a <- t(sapply(r,function(x) {
        colSums(hapSum[x[1]:x[2],])/nSNP
    }))
    b=a<(N*hapFreqMin)
    ## for each k, for each region, fill in
    ## within each continuous region, fill in with respect to frequencies of all other haplotypes
    ##
    for(k in 1:K) {
        ## change into intervals
        z1=b[,k]
        x=(1:(length(z1)-1))[diff(z1)!=0]
        start=c(1,x+1)
        end=c(x,length(z1))
        ## only take "bad" regions
        start=start[z1[start]==TRUE]
        end=end[z1[end]==TRUE    ]
        if(length(start)>0) {
            for(iR in 1:length(start)) {
                ## now - in each region, get sums of haplotypes
                p=colSums(matrix(a[start[iR]:end[iR],],ncol=K)); p=p/sum(p)
                replacement=sample(K,1,prob=p)
                y=r[[start[iR]]][1]:r[[end[iR]]][2]
                eHapsCurrent[y,k]=eHapsCurrent[y,replacement]
            }
        }
    }
    print(paste("Refill infrequently used haplotypes - ",round(100*sum(b==TRUE)/(nrow(b)*ncol(b)),1),"% of regions replaced",sep=""))
    return(list(eHapsCurrent = eHapsCurrent))
}








## the cigar string comes as a long character vector
## reconfigure this into lists to be able to use in C++
## in retrospect, could probably have done this as part
## of other giant C++ function, but oh well
##
## cigarRead is a vector of cigars
## like "101M" or "10M2I3M4D" etc
split_cigar_string <- function(
    cigarRead,
    useSoftClippedBases = FALSE,
    posRead = NULL,
    seqRead = NULL
) {
    ## now in C++
    splitCigarRead <- cpp_cigar_split(unlist(cigarRead))
    ##
    ## deal with soft and hard clipped bases
    ##
    s <- sapply(splitCigarRead, function(x) sum(x[[2]]=="S"))
    h <- sapply(splitCigarRead, function(x) sum(x[[2]]=="H"))
    if(sum(s)>0) {
      if(useSoftClippedBases==TRUE) {
        ## adjust position at the start appropriately if first is S
        w <- sapply(splitCigarRead, function(x) x[[2]][1]=="S")
        ## get position of those that start with an S
        x <- sapply(splitCigarRead[w], function(x) x[[1]][1])
        ## substract the starts back by that
        posRead[w] <- posRead[w] - x
        ## change all the S's to M's
        splitCigarRead <- lapply(splitCigarRead, function(x) {
          x[[2]][x[[2]]=="S"]="M"
          return(x)
        })
      } else {
        ## here we want to get rid of them
        ## we need to effectively hard clip them
        ## the sequence at the end is irrelevant
        w <- sapply(splitCigarRead, function(x) x[[2]][1]=="S")
        if (sum(w) > 0) {
          seqRead[w] <- sapply(which(w), function(i) {
            s1 <- seqRead[i]
            scr <- splitCigarRead[[i]]
            to_remove <- scr[[1]][1]
            return(substr(s1, 1+to_remove,nchar(s1)))
          })
          ## now - remove them from the start
          splitCigarRead[w] <- lapply(splitCigarRead[w], function(x) {
            k <- x[[2]]!="S"
            return(list(x[[1]][k], x[[2]][k]))
          })
        }
      }
    }
    ## done looking at soft clipped bases
    lengthOfSplitCigarRead <- as.integer(unlist(lapply(splitCigarRead,function(x) length(x[[1]])-1))) # 0 BASED
    ## hmm, looks like I can basically ignore hard clipped ones
    ## the src function will skip them
    ## and they aren't part of the alignment
    return(
        list(
            splitCigarRead = splitCigarRead,
            lengthOfSplitCigarRead = lengthOfSplitCigarRead,
            posRead = posRead,
            seqRead = seqRead
        )
    )
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
    chrEnd,
    method_to_use_to_get_raw_data
) {
    if (length(bam_files) > 0) {
        bamName <- bam_files[iBam]
        if (file.exists(bamName) == FALSE)
            stop(paste0("cannot find file:", bamName))
    } else if (length(cram_files) > 0) {
        cramName <- cram_files[iBam]
        if (file.exists(cramName) == FALSE)
            stop(paste0("cannot find file:", cramName))
        if (method_to_use_to_get_raw_data == "Rsamtools") {
            tempBam <- paste0(tempdir, "/sample", iBam, ".bam")
            ## note - chrStart and chrEnd are appropriate
            system(paste0(
                "samtools view -h ", cramName,
                " -T ", reference, " ",
                chr, ":", chrStart, "-", chrEnd, " | ",
                "samtools view -b > ",
                tempBam
            ))
            system(paste0("samtools index ", tempBam))
            bamName <- tempBam
        } else {
            bamName <- cramName
        }
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




## with data out of Rsamtools
## build first sampleReads from C++
## with one entry per line of the BAM file,
## when it intersects a SNP
get_sampleReadsRaw <- function(
    sampleData,
    bqFilter,
    iSizeUpperLimit,
    T,
    L,
    pos,
    useSoftClippedBases
) {

    ## get more info read as well
    mapq <- as.integer(sampleData$mapq)
    isize <- sampleData$isize
    ## figure out which ones we dont want - remove
    keep <- FALSE == (
        mapq < bqFilter |
        ( abs(isize) > iSizeUpperLimit & is.na(isize)==FALSE)
    )
    cigarRead <- sampleData$cigar
    for(l in c("P", "=", "X"))
        keep[grep(l, cigarRead)] <- FALSE
    mapq <- mapq[keep]
    isize <- isize[keep]
    posRead <- as.integer(sampleData$pos)[keep]
    cigarRead <- sampleData$cigar[keep]
    seqRead <- as.character(sampleData$seq)[keep]
    qualRead <- as.character(sampleData$qual)[keep]
    qname <- sampleData$qname[keep]
    strand <- as.character(sampleData$strand)[keep]
    rm(sampleData)

    out <- split_cigar_string(
        cigarRead,
        useSoftClippedBases,
        posRead,
        seqRead
    )
    splitCigarRead <- out$splitCigarRead
    lengthOfSplitCigarRead <- out$lengthOfSplitCigarRead
    posRead <- out$posRead
    seqRead <- out$seqRead

    readLength <- sapply(
        splitCigarRead,
        function(x)
        sum(
            x[[1]][FALSE==is.na(match(x[[2]], c("M","D")))]
        )
    )
    readStart <- posRead # includes this first base
    readEnd <- readStart + readLength - 1 # includes final base

    ##
    ## push through c++ function
    ##
    ref <- as.character(pos[, "REF"])
    alt <- as.character(pos[, "ALT"])
    sampleReadsRaw <- reformatReads(
        mapq = mapq,
        readStart = readStart,
        readEnd = readEnd,
        numberOfReads = as.integer(length(mapq)),
        T = T,
        L = L,
        seqRead = seqRead,
        qualRead = qualRead,
        ref = ref,
        alt = alt,
        splitCigarRead = splitCigarRead,
        lengthOfSplitCigarRead = lengthOfSplitCigarRead,
        bqFilter = bqFilter,
        verbose = as.integer(0)
    )

    return(
        list(
            sampleReadsRaw = sampleReadsRaw,
            qname = qname,
            strand = strand
        )
    )
}




## previous sampleReads comes straight from C++
## and is not linked by reads
## here, match by reads
merge_reads_from_sampleReadsRaw <- function(
    sampleReadsRaw,
    qname,
    strand
) {
    ## wif: 1-based which read it came from
    wif <- sapply(sampleReadsRaw, function(x) x[[5]]) + 1
    qnameUnique <- unique(qname[wif])
    qnameInteger <- match(qname[wif], qnameUnique)
    ord <- order(qnameInteger) - 1 ## 0-based for C++
    qnameInteger_ord <- c(qnameInteger[ord + 1], - 2)
    ## qnameUnique is unique read names of used reads
    ## qnameInteger is integer version of that
    ## ord is 0-based ordered version of qnameInteger
    ## qnameInteger_ord is the ordered version of qnameInteger
    sampleReads <- cpp_read_reassign(
        ord,
        qnameInteger_ord,
        sampleReadsRaw,
        as.integer(0)
    )
    sampleReadsInfo <- data.frame(
        qname = qnameUnique,
        strand = strand[match(qnameUnique, qname)] ## problem if varies
    )
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
    print(
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

## for some given chrStart, chrEnd
## and some desired window size
## choose how to do chunking
determine_loading_windows <- function(
    chrStart,
    chrEnd,
    chrLength,
    width_of_loading_window = 1e6
) {
    if (is.na(chrStart)) {
        chrEnd <- chrLength
        chrStart <- 1
    }
    load_window_start <- seq(chrStart, chrEnd, width_of_loading_window)
    load_window_end <- c(load_window_start[-1] - 1, chrEnd)
    return(
        list(
            start = load_window_start,
            end = load_window_end
        )
    )
}




combineReadsAcrossRegions <- function(sampleReadsAcrossRegions) {

    ## there can be a different number of reads compared to length of qname
    n_regions <- length(sampleReadsAcrossRegions)
    reads_per_region_srr <- sapply(sampleReadsAcrossRegions, function(x) return(length(x$sampleReadsRaw)))
    reads_per_region_qname <- sapply(sampleReadsAcrossRegions, function(x) return(length(x$qname)))
    reads_before_region_srr <- c(0, cumsum(reads_per_region_srr))
    reads_before_region_qname <- c(0, cumsum(reads_per_region_qname))

    n_reads <- sum(reads_per_region_srr)

    if (sum(n_reads) == 0) {

        sampleReadsRaw <- NULL
        qname <- NULL
        strand <- NULL

    } else {

        read_comes_from_region <- unlist(lapply(1:n_regions, function(i_region) return(rep(i_region, reads_per_region_srr[i_region]))))
        sampleReadsRaw <- lapply(1:n_reads, function(iRead) {
            i_region <- read_comes_from_region[iRead]
            iReadInRegion <- iRead - reads_before_region_srr[i_region]
            sampleRead <- sampleReadsAcrossRegions[[i_region]]$sampleReadsRaw[[iReadInRegion]]
            sampleRead[[5]] <- sampleRead[[5]] + reads_before_region_qname[i_region]
            return(sampleRead)
        })
        qname <- unlist(lapply(sampleReadsAcrossRegions, function(x) x$qname))
        strand <- unlist(lapply(sampleReadsAcrossRegions, function(x) x$strand))

    }
    return(
        list(
            sampleReadsRaw = sampleReadsRaw,
            qname = qname,
            strand = strand
        )
    )
}



get_sample_data_from_Rsamtools <- function(
    chr,
    window_start,
    window_end,
    bamName
) {
    ## necessary Rsamtools stuff
    flag <- scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = FALSE, hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA, isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
    ##isNotPrimaryRead = NA
    what <- c("qname","strand","pos","seq","qual","cigar","isize","mapq")
    eval(parse(text = (
        paste0("which = RangesList(\"",chr,"\"=",
               "IRanges(", window_start,",", window_end,"))")
    )))
    param <- ScanBamParam(flag = flag, which = which, what = what)
    idx <- get_index_for_bamName(bamName)
    sampleData <- scanBam(
        file = bamName,
        index = idx,
        param = param
    )[[1]]
    return(sampleData)
}



## takes a lot of inputs
## converts a BAM file for one subject
## into sampleReads, saved into
## file_sampleReads(inputdir, iBam, regionName)
loadBamAndConvert <- function(
    iBam,
    L,
    pos,
    T,
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
    save_sampleReadsInfo = FALSE,
    width_of_loading_window = 1000000, ## what sized chunks to load things in, for RAM reasons
    method_to_use_to_get_raw_data = "SeqLib"
) {

    sampleReadsInfo <- NULL ## unless otherwise created

    if ((iBam %% 100) == 0)
        print(
            paste0(
                "load and convert bam ",
                iBam, " of ",N , ", ", date()
            )
        )

    file_name <- get_bam_name_and_maybe_convert_cram(
        iBam,
        bam_files,
        cram_files,
        reference,
        tempdir,
        chr,
        chrStart,
        chrEnd,
        method_to_use_to_get_raw_data
    )

    loading_windows <- determine_loading_windows(
        chrStart, chrEnd, chrLength, width_of_loading_window
    )

    sampleReadsAcrossRegions <- lapply(1:length(loading_windows$start), function(i_region) {
        window_start <- loading_windows$start[i_region]
        window_end <- loading_windows$end[i_region]

        if (method_to_use_to_get_raw_data == "Rsamtools") {
            sampleData <- get_sample_data_from_Rsamtools(
                chr = chr,
                window_start = window_start,
                window_end = window_end,
                bamName = file_name
            )
        } else if (method_to_use_to_get_raw_data == "SeqLib") {
            sampleData <- get_sample_data_from_SeqLib(
                region = paste0(chr, ":", window_start, "-", window_end),
                file_name = file_name,
                reference = reference
            )
        }
            
        if (length(sampleData$qname) > 0) {
            out <- get_sampleReadsRaw(
                sampleData,
                bqFilter,
                iSizeUpperLimit,
                T,
                L,
                pos,
                useSoftClippedBases
            )
        } else {
            out <- list(sampleReadsRaw = NULL, qname = NULL, strand = NULL)
        }
        return(
            list(
                sampleReadsRaw = out$sampleReadsRaw,
                qname = out$qname,
                strand = out$strand
            )
        )
    })

    out <- combineReadsAcrossRegions(sampleReadsAcrossRegions)
    sampleReadsRaw <- out$sampleReadsRaw
    qname <- out$qname
    strand <- out$strand

    if (length(sampleReadsRaw) == 0) {
        sampleReads <- get_fake_sampleReads(
            name = sampleNames[iBam],
            mess = "has no informative reads"
        )
    } else {
        out <- merge_reads_from_sampleReadsRaw(
            sampleReadsRaw,
            qname,
            strand
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

    if (method_to_use_to_get_raw_data == "Rsamtools")
        if (length(cram_files) > 1)
            system(paste0("rm ", tempdir, "sample", iBam, ".ba*"))

    out <- downsample(
        sampleReads = sampleReads,
        iBam = iBam,
        downsampleToCov = downsampleToCov,
        sampleNames = sampleNames,
        sampleReadsInfo = sampleReadsInfo
    )
    sampleReads <- out$sampleReads
    sampleReadsInfo <- out$sampleReadsInfo

    save(
        sampleReads,
        file = file_sampleReads(inputdir, iBam, regionName)
    )

    if (save_sampleReadsInfo)
        save(
            sampleReadsInfo,
            file = file_sampleReadsInfo(inputdir, iBam, regionName)
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


file_dosages <- function(
  dir,
  iBam,
  regionName
) {
    return(file.path(
        dir,
        paste0("sample.", iBam, ".input.",regionName,".dosage.RData")
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



shrinkReads <- function(
  N,
  nCores,
  sampleRange,
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

    if(regenerateInput == TRUE | regionName == originalRegionName) {

        print(paste0("Copying files onto tempdir, ", date()))
        if (is.na(inputBundleBlockSize)) {
            file_with_files_to_transfer <- file.path(tempdir, "files_to_transfer.txt")
            command1 <- paste0(
                '(cd ', inputdir, ' && find . -name "',
                'sample.*.input.', regionName, '.RData',
                '" > ', file_with_files_to_transfer, ')'
            )
            system(command1)
            command2 <- paste0("rsync -a --files-from=", file_with_files_to_transfer, " ", inputdir, " ", tempdir)
            system(command2)
        } else {
            file_with_files_to_transfer <- file.path(tempdir, "files_to_transfer.txt")
            command1 <- paste0(
                '(cd ', inputdir, ' && find . -name "',
                "bundledSamples.*-*.", regionName, ".RData",
                '" > ', file_with_files_to_transfer, ')'
            )
            system(command1)
            command2 <- paste0("rsync -a --files-from=", file_with_files_to_transfer, " ", inputdir, " ", tempdir)
            system(command2)
        }
        print(paste0("Done copying files onto tempdir, ", date()))
    } else {

        print(paste0("Begin shrink reads, ", date()))
        sampleRanges <- getSampleRange(N = N, nCores = nCores)
        if(environment == "server") {
            out <- mclapply(sampleRanges,FUN=shrinkReads_on_range, mc.cores=nCores, originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL)
        }
        if(environment == "cluster") {
            cl <- makeCluster(nCores, type = "FORK")
            out <- parLapply(cl, sampleRanges,fun=shrinkReads,originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL)
            stopCluster(cl)
        }
        error_check <- sapply(out, class) == "try-error"
        if (sum(error_check) > 0) {
            print(out[[which(error_check)[1]]])
            stop("There has been an error generating the input. Please see error message above")
        }
        print(paste0("End shrink reads, ", date()))

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
            file = file_sampleReads(tempdir, iSample, regionName)
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





### function to downsample BAMs as necessary
downsample <- function(
    sampleReads,
    iBam,
    downsampleToCov,
    sampleNames,
    sampleReadsInfo
) {
    ## first, figure out if reads need to be downsampled
    y <- unlist(lapply(sampleReads,function(x) x[[2]]))
    ## quick counts of each
    z <- table(y)
    z1 <- as.numeric(names(z)[z>downsampleToCov])
    if(length(z1)>0) {
        ## match each against
        z2 <- match(z1,y)
        z3 <- z2[is.na(z2)==FALSE]
        ## get counts of each
        counts=as.numeric(z[match(z1,names(z))])
        ## get positions of these elements
        ## keep a list of those to remove
        toRemove=array(FALSE,length(sampleReads))
        for(i in 1:length(z1)) {
            w=z3[i] # which read position to work on
            c=counts[i] # the counts
            range=w:(w+c-1)# the position range
            ## downsample here
            toRemove[range]=TRUE
            toRemove[sample(range,downsampleToCov,replace=FALSE)]=FALSE
            ## done!
        }
        ## print out
        print(paste0("WARNING - downsample sample ",sampleNames[iBam]," - ",sum(toRemove)," of ",length(sampleReads)," reads removed "))
        ##
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



### split
splitReadsCompletely=function(N,nCores,tempdir,regionName,environment)
{
  # split into 48 cores
  f=function(iSample,tempdir,regionName) {
    load(file=file_sampleReads(tempdir, iSample, regionName))   # split!
    x3=unlist(lapply(sampleReads,function(x) x[[3]]))
    x4=unlist(lapply(sampleReads,function(x) x[[4]]))
    sr=lapply(1:length(x3),function(i) {
      list(0,x4[i],x3[i],x4[i])
    })
    sampleReads=sr[order(x4)]
    save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName))
  }
  if(environment=="server")
  {
    out2=mclapply(1:N,mc.cores=nCores,FUN=f,tempdir=tempdir,regionName=regionName)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, 1:N, fun=f,tempdir=tempdir,regionName=regionName)
    stopCluster(cl)
  }
  # done
}



# convert scaled base quality to probabilities
convertScaledBQtoProbs <- function(x) {
  # for a matrix with 1 column, with base qualities, get the probabilities of the ref or alternate
  # recall < 0 is ref, >1 is alternate
  # output is a matrix, two columns
  n=nrow(x)
  o=array(0,c(n,2))
  e=10**(-abs(x)/10)
  # fill in
  w=cbind(1:n,2-as.integer(x<0))
  o[w]=1-e
  w[,2]=3-w[,2]
  o[w]=e*1/3
  return(o)
}


buildAlleleCount_subfunction <- function(i,T,tempdir,regionName, bundling_info) {
    alleleCount <- array(0,c(T,3))
    bundledSampleReads <- NULL
    for(iSample in i[1]:i[2])
    {
      out <- get_sampleReads_from_dir_for_sample(
        dir = tempdir,
        regionName = regionName,
        iSample = iSample,
        bundling_info = bundling_info,
        bundledSampleReads = bundledSampleReads
      )
      sampleReads <- out$sampleReads
      bundledSampleReads <- out$bundledSampleReads
      # get positions and intensities
      a=unlist(sapply(sampleReads,function(x) x[[3]]))
      b=unlist(sapply(sampleReads,function(x) x[[4]]))
      bqProbs=convertScaledBQtoProbs(matrix(a,ncol=1))
      # y is numeric, z = integer counts (can have 0),
      c1=increment2N(y=as.numeric(bqProbs[,1]),z=as.numeric(b),yT=as.integer(nrow(bqProbs)),xT=as.integer(T-1))
      c2=increment2N(y=as.numeric(bqProbs[,2]),z=as.numeric(b),yT=as.integer(nrow(bqProbs)),xT=as.integer(T-1))
      alleleCount[,1]=alleleCount[,1]+c2 # fixed nov 6 2015 - was backward
      alleleCount[,2]=alleleCount[,2]+c1+c2
    }
    return(list(alleleCount=alleleCount))
}


### build allele count matrix from input RData files
buildAlleleCount <- function(
    T,
    N,
    nCores,
    regionName,
    tempdir,
    environment,
    bundling_info
) {
    x3 <- getSampleRange(N = N, nCores = nCores)
    print(paste0("Generate allele count, ", date()))

    if(environment=="server") {
        out2 <- mclapply(x3,mc.cores=nCores,FUN=buildAlleleCount_subfunction,T=T,tempdir=tempdir,regionName=regionName, bundling_info = bundling_info)
    }
    if(environment=="cluster") {
        cl <- makeCluster(length(x3), type = "FORK")
        out2 <- parLapply(cl, x3, fun=buildAlleleCount_subfunction,T=T,tempdir=tempdir,regionName=regionName, bundling_info = bundling_info)
        stopCluster(cl)
    }

    error_check <- sapply(out2, class) == "try-error"
    if (sum(error_check) > 0) {
        print(out2[[which(error_check)[1]]])
        stop("There has been an error generating the input. Please see error message above")
    }
    alleleCount <- array(0,c(T,3))
    for(i in 1:length(x3)) {
        alleleCount <- alleleCount + out2[[i]][["alleleCount"]]
    }
    alleleCount[, 3] <- alleleCount[, 1] / alleleCount[, 2]

    print(paste0("Quantiles across SNPs of per-sample depth of coverage"))
    print(quantile(alleleCount[,2] / N, prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))
    print(paste0("Done generating allele count, ", date()))

    alleleCount[is.na(alleleCount[,3]),3] <- 0 # not sure why these would exist anymore
    return(alleleCount)
}






downsampleToFraction_a_range <- function(
  sampleRange,
  tempdir,
  regionName,
  downsampleFraction,
  bundling_info
) {
  bundledSampleReads <- NULL
  for(iSample in sampleRange[1]:sampleRange[2])
  {
    out <- get_sampleReads_from_dir_for_sample(
      dir = tempdir,
      regionName = regionName,
      iSample = iSample,
      bundling_info = bundling_info,
      bundledSampleReads = bundledSampleReads
    )
    sampleReads <- out$sampleReads
    bundledSampleReads <- out$bundledSampleReads

    # subset
    keep=runif(length(sampleReads))<downsampleFraction
    # keep 1 no matter what
    if(sum(keep)==0) keep[sample(length(keep),1)]=TRUE
    sampleReads=sampleReads[keep]
    # save result back to disk
    save(sampleReads,file=file_sampleReads(tempdir, iSample, regionName))

    # bundle back together
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


### build allele count matrix from input RData files
downsampleToFraction <- function(N,nCores,downsampleFraction,regionName,tempdir,environment, bundling_info)
{
  print(paste0("Begin downsampling reads, ",date()))
  # do first loop
  sampleRanges <- getSampleRange(N = N, nCores = nCores)
  if(environment=="server")
  {
    out2=mclapply(sampleRanges,mc.cores=nCores,FUN=downsampleToFraction_a_range,tempdir=tempdir,regionName=regionName,downsampleFraction=downsampleFraction, bundling_info = bundling_info)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, sampleRanges, fun=downsampleToFraction_a_range,tempdir=tempdir,regionName=regionName,downsampleFraction=downsampleFraction, bundling_info = bundling_info)
    stopCluster(cl)
  }
  error_check <- sapply(out2, class) == "try-error"
  if (sum(error_check) > 0) {
    print(out2[[which(error_check)[1]]])
    stop("There has been an error downsampling the reads. Please see error message above")
  }
  print(paste0("Done downsampling reads, ",date()))
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
    if (sum(sampleRange == c(1, 1)) == 2)
        return(c(0, 1))
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



# for N, nCores, blockSize
# return a matrix and a list
# -- for the matrix,
# with N rows, one per sample, in sampleNames order
# where first column is what core it is in
# where second column is what bundle in that core it is in
# where third column is position within that bundle
# -- for the list,
# bundlingList[[iCore]][[iBundle]] gives start and end
get_bundling_position_information <- function(
  N,
  nCores,
  blockSize
) {
  if (is.na(blockSize))
    return(NULL)
  # one entry per core, the range
  x3 <- getSampleRange(N = N, nCores = nCores)
  out <- lapply(1:nCores, function(i_core) {
    x2 <- getOutputBlockRange(
      sampleRange = x3[[i_core]],
      outputBlockSize = blockSize
    )
    # for each of these block ranges, get the bundle
    inputBundle <- unlist(lapply(1:(length(x2) - 1), function(i_sub_block) {
      return(rep(i_sub_block, x2[i_sub_block + 1] - x2[i_sub_block]))
    }))
    # within each core and bundle, get the position
    positionWithinBundle <- unlist(lapply(1:(length(x2) - 1), function(i_sub_block) {
      return(1:(x2[i_sub_block + 1] - x2[i_sub_block]))
    }))
    # this is 0-based - transform to 1 based and return
    # also - return a list
    start <- x2[-length(x2)] + 1
    end <- x2[-1]
    bundleList <- lapply(1:length(start), function(i) {
      c(start[i], end[i])
    })
    return(list(
      inputCore = rep(i_core, diff(x3[[i_core]]) + 1),
      inputBundle = inputBundle,
      positionWithinBundle = positionWithinBundle,
      bundleList = bundleList
    ))
  })
  # turn into a matrix
  iCore <- unlist(lapply(out, function(x) x$inputCore))
  iBundle <- unlist(lapply(out, function(x) x$inputBundle))
  iPosition <- unlist(lapply(out, function(x) x$positionWithinBundle))
  bundlePosition = cbind(
    iCore = iCore,
    iBundle = iBundle,
    iPosition = iPosition,
    last_in_bundle = 0
  )
  # also - add "last_in_bundle
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
# if pseudoHaploid, do this, before completeSampleIteration
initialize_readProbs <- function(
  N,
  nCores,
  pseudoHaploidModel,
  tempdir,
  bundling_info,
  regionName,
  initial_min_hapProb,
  initial_max_hapProb
) {
  print(paste0("Initialize readProbs, ", date()))
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
        bundledSampleReads = bundledSampleReads
      )
      sampleReads <- out$sampleReads
      bundledSampleReads <- out$bundledSampleReads
      # now - can make some hapProbs
      out <- get_default_hapProbs(
          pseudoHaploidModel,
          sampleReads,
          initial_min_hapProb,
          initial_max_hapProb
      )
      pRgivenH1 <- out$pRgivenH1
      pRgivenH2 <- out$pRgivenH2
      srp <- out$srp
      delta <- out$delta
      save(
        pRgivenH1, pRgivenH2, srp, delta,
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
    print(out[[which(error_check)[1]]])
    stop("There has been an error generating the input. Please see error message above")
  }
  print(paste0("Done initializing readProbs, ", date()))
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
    srp=unlist(lapply(sampleReads,function(x) x[[2]]))
    delta=0
    if(pseudoHaploidModel==6) {
        delta=array(0,length(srp)) # may or may not use?
        delta=a1+(b1-a1)*runif(length(sampleReads))
    }
    return(list(
        pRgivenH1 = pRgivenH1,
        pRgivenH2 = pRgivenH2,
        srp = srp,
        delta = delta
    ))
}



subset_of_complete_iteration <- function(sampleRange,tempdir,chr,K,T,priorCurrent,eHapsCurrent_t,alphaMatCurrent_t,sigmaCurrent,maxDifferenceBetweenReads,whatToReturn,Jmax,highCovInLow,iteration,method,nsplit,expRate,minRate,maxRate,gen,outputdir,pseudoHaploidModel,outputHaplotypeProbabilities,switchModelIteration,regionName,restartIterations,refillIterations,hapSumCurrentL,outputBlockSize, bundling_info, transMatRate_t, x3, N, shuffleHaplotypeIterations, niterations, L, samples_with_phase, nbreaks, breaks) {

    ## initialize bundling variables
    bundledSampleReads <- NULL
    bundledSampleProbs <- NULL
    ##
    ## determine whether its a special iteration
    ##
    readsSplit <- array(0, N)
    readsTotal <- array(0, N)
    ##
    ## initialize refilling iteration
    ##
    refillList <- NA
    refillMatrixList <- NA
    file <- paste(outputdir,"RData/refillList.",regionName,".iteration.",iteration-1,".RData",sep="")
    if(file.exists(file)) {
      load(file)
      refillMatrixList=lapply(1:K,function(k) return(array(0,c(nrow(refillList[[k]]),K))))
      refillFromList=lapply(1:K,function(k) {
        # yes, we're using "1" for both
        # what we do is look at the start of the rare section
        # look what it was, then what it becomes
        from=sapply(refillList[[k]][,1]-50,function(x) max(1,x))
        to=sapply(refillList[[k]][,1]+50,function(x) min(T,x))
        return(list(from=from,to=to))
      })
    }
    ##
    ## initialize output matrices here
    ##
    priorSum <- array(0,K)
    alphaMatSum_t <- array(0, c(K, T - 1))
    gammaSum_t <- array(0, c(K, T, 2))
    hapSum_t <- array(0,c(K, T))
    fromMat <- array(0, c(nbreaks, K, K))
    restartMatrixList <- as.list(1:N)
    ##
    ## if final iteration - make count matrices
    ##
    hweCount <- NA
    afCount <- NA
    infoCount <- NA
    if (iteration==niterations) {
        hweCount = array(0,c(T,3))
        infoCount = array(0,c(T,2))
        afCount = array(0,T)
    }
    ##
    ## set output block size for VCF pasting
    ##
    list_of_vcf_columns_to_out <- as.list(1:N)
    outputBlockRange <- getOutputBlockRange(
        sampleRange,
        outputBlockSize
    )
    ##
    ## 1 - run forward backwards
    ## 2 - get update pieces
    ## 3 - if necessary add to output
    ##
    whatToReturnOriginal=whatToReturn
    for(iSample in sampleRange[1]:sampleRange[2]) {
        ##
        ##
        ## 1 - run forward backwards
        ##
        ##
        whatToReturn <- whatToReturnOriginal
        if(is.na(match(iSample,highCovInLow))==FALSE)
            whatToReturn=as.integer(2)
        ## load both sample reads and read probabilities
        out <- get_sampleReads_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = iSample,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads
        ##
        ## if pseudoHaploid method, get pRgivenH1, etc.
        ##
        if(method=="pseudoHaploid") {
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
            delta <- out$delta
            bundledSampleProbs <- out$bundledSampleProbs
        }
        ##
        ##
        ## run single sample iteration here
        ##
        ##
        ## run either once (default) or more than once
        nor <- 1 # number of runs - 1 for normal diploid method
        ## run both haplotypes
        if(method=="pseudoHaploid")
            nor=2
        fbsoL=as.list(1:nor)
        for(iNor in 1:nor) {
            if(method=="pseudoHaploid") {
                if (iNor==1) {
                    pRgivenH1L <- pRgivenH1
                    pRgivenH2L <- pRgivenH2
                    deltaL <- delta
                }
                if (iNor==2) {
                    pRgivenH1L <- pRgivenH2
                    pRgivenH2L <- pRgivenH1
                    deltaL <- 1 - delta
                }
                fbsoL[[iNor]]=forwardBackwardHaploid(
                    sampleReads = sampleReads,
                    nReads = as.integer(length(sampleReads)),
                    Jmax = as.integer(Jmax),
                    pi = priorCurrent,
                    pRgivenH1 = pRgivenH1L,
                    pRgivenH2 = pRgivenH2L,
                    pState = matrix(0, c(2, 2)), ## hmm, deprecated?
                    eHaps_t = eHapsCurrent_t,
                    alphaMat_t = alphaMatCurrent_t,
                    transMatRate_t = transMatRate_t,
                    maxDifferenceBetweenReads = as.double(maxDifferenceBetweenReads),
                    whatToReturn = whatToReturn,
                    suppressOutput = as.integer(1),
                    model = as.integer(pseudoHaploidModel))
            }
            if(method=="diploid") {
                fbsoL[[iNor]] <- forwardBackwardDiploid(
                    sampleReads = sampleReads,
                    nReads = as.integer(length(sampleReads)),
                    pi = priorCurrent,
                    eHaps_t = eHapsCurrent_t,
                    alphaMat_t = alphaMatCurrent_t,
                    transMatRate_t = transMatRate_t,
                    Jmax = Jmax,
                    maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                    whatToReturn = whatToReturn,
                    suppressOutput = as.integer(1)
                )
            }
            ##
            ## get update pieces
            ##
            which <- fbsoL[[iNor]]$gammaUpdate_t[1, , 1] > 0
            gammaSum_t[, which, ] <- gammaSum_t[, which , ] +
                fbsoL[[iNor]]$gammaUpdate_t[, which, ]
            alphaMatSum_t <- alphaMatSum_t +
                fbsoL[[iNor]]$jUpdate_t
            priorSum <- priorSum+fbsoL[[iNor]]$gammaUpdate_t[, 1, 2]
        }
        ##
        ##
        ##
        ##
        ## 2 - do various heuristics
        ##
        ##
        ##
        ##
        ##
        if (method == "pseudoHaploid") {
            for(iNor in 1:2)
                fbsoL[[iNor]]$gammaK_t <- fbsoL[[iNor]]$gamma_t
        }else if(method=="diploid") {
            y <- fbsoL[[iNor]]$gamma_t
            gamma_t <- y[1:K, ]
            for(ii in 2:K)
                gamma_t <- gamma_t + y[1:K+(ii-1)*K, ]
            fbsoL[[iNor]]$gammaK_t <- gamma_t
        }
        hapSum_t <- hapSum_t + fbsoL[[iNor]]$gammaK_t
        if(nbreaks>0) {
            ## for each sample, get the changes between every 100
            for(iNor in 1:nor) {
                hp=matrix(t(fbsoL[[1]]$gamma_t[, breaks]) ,ncol = K)
                for(f in 1:nbreaks)
                    fromMat[f,,]=fromMat[f,,] + hp[f,] %*% t(hp[f+1,])
            }
        }
        ##
        ## now - consider whether to re-start - ie fix gammas
        ##
        out <- restartSample(sampleReads = sampleReads, srp = srp, pRgivenH1 = pRgivenH1, pRgivenH2 = pRgivenH2, fbsoL=fbsoL,T=T,eHapsCurrent_t=eHapsCurrent_t,sigmaCurrent=sigmaCurrent,alphaMatCurrent_t=alphaMatCurrent_t,Jmax=Jmax,priorCurrent=priorCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,outputdir=outputdir,iteration = iteration,restartIterations = restartIterations, method = method)
        restartMatrixList[[iSample]] <- out$m
        fbsoL <- out$fbsoL
        ##
        ## update pseudo diploid probabilities
        ##
        if(method=="pseudoHaploid") {
            r=pseudoHaploidUpdate(pRgivenH1=pRgivenH1,pRgivenH2=pRgivenH2,fbsoL=fbsoL,pseudoHaploidModel=pseudoHaploidModel,srp=srp,delta=delta, K = K)
            pRgivenH1=r$pRgivenH1
            pRgivenH2=r$pRgivenH2
            delta=r$delta
            save(srp,pRgivenH1,pRgivenH2,delta,file=paste(tempdir,"sample.",iSample,".readProbs.",regionName,".RData",sep=""))
        }
        ##
        ## check high coverage samples if applicable
        ##
        if(is.na(match(iSample,highCovInLow))==FALSE) {
            dosage <- array(0, T)
            if (method == "pseudoHaploid") {
                for(iNor in 1:nor)
                    dosage <- dosage + (colSums(fbsoL[[iNor]]$gamma_t* eHapsCurrent_t))
            } else {
                dosage <- fbsoL[[iNor]]$dosage
            }
            if(nor==2) dosage <- dosage/2
            save(dosage, file=file_dosages(tempdir, iSample, regionName))
        } # end of check on high coverage sample
        ##
        ## check phasing samples if applicable
        ##
        if (method == "pseudoHaploid" &&
            is.na(match(iSample, samples_with_phase))==FALSE) {
            haplotypes <- array(0, c(T, 2))
            for(i_hap in 1:2)
                haplotypes[, i_hap] <- colSums(fbsoL[[i_hap]]$gamma_t * eHapsCurrent_t)
            save(haplotypes, file=file_haplotypes(tempdir, iSample, regionName))
        }
        ##
        ## do read splitting if the correct iteration
        ##
        if ( nsplit > 0) {
            if(method=="pseudoHaploid") {
                gammaK_t <- fbsoL[[1]]$gammaK_t + fbsoL[[2]]$gammaK_t
            } else {
                gammaK_t <- fbsoL[[1]]$gammaK_t
            }
            out <- findRecombinedReadsPerSample(
                gammaK_t = fbsoL[[1]]$gammaK_t,eHapsCurrent_t=eHapsCurrent_t,K=K,L=L,iSample=iSample,verbose=FALSE,sampleReads=sampleReads,tempdir=tempdir,regionName=regionName
            )
            readsTotal[iSample] <- out$readsTotal
            readsSplit[iSample] <- out$readsSplit
        }
        ##
        ##
        ##
        ##
        ## 3 - output as necessary
        ##
        ##
        ##
        if(iteration==niterations) {
            ## get allele counts for HWE, info calcs, AF
            ## get most likely genotype - add to count
            for(i in 1:length(fbsoL)) {
                if ( method == "diploid") {
                    gp <- fbsoL[[1]]$genProbs
                } else {
                    gp <- array(0, c(T, 3))
                    g10 <- colSums(fbsoL[[1]]$gamma_t * (1-eHapsCurrent_t))
                    g20 <- colSums(fbsoL[[2]]$gamma_t * (1-eHapsCurrent_t))
                    gp[, 1] <- g10 * g20
                    gp[, 2] <- g10 * (1-g20) + (1-g10) * g20
                    gp[, 3] <- (1-g10) * (1-g20)
                }
            }
            ## info counts
            eij <- gp[,2] + 2 * gp[,3]
            fij <- gp[,2] + 4 * gp[,3]
            infoCount[,1] <- infoCount[,1] + eij
            infoCount[,2] <- infoCount[,2] + (fij - eij**2)
            ## do counts for HWE
            w <- get_max_gen_rapid(gp)
            hweCount[w] <- hweCount[w]+1
            ## get counts for allele frequency
            afCount <- afCount + (gp[,2] + 2*gp[,3]) / 2
            ## if pseudo-haploid, get probabilities
            ## disable outputting for now
            if (method == "pseudoHaploid" && 1 == 0) {
                read_proportions <- estimate_read_proportions(
                    sampleReads = sampleReads,
                    pRgivenH1 = pRgivenH1,
                    pRgivenH2 = pRgivenH2,
                    T = T
                )
            } else {
                read_proportions <- NULL
            }
            ##
            ## add column to VCF, and possibly write block to disk
            ##
            list_of_vcf_columns_to_out[[iSample]] <- make_column_of_vcf(gp, read_proportions)
            i_core <- match(sampleRange[1], sapply(x3, function(x) x[[1]]))
            iBlock <- match(iSample, outputBlockRange)
            if (iBlock > 1 & is.na(iBlock)==FALSE) {
                iSample_list <- (outputBlockRange[iBlock-1]+1):outputBlockRange[iBlock]
                write_block_of_vcf(i_core, iBlock, list_of_vcf_columns_to_out, iSample_list = iSample_list, outputdir, regionName)
                list_of_vcf_columns_to_out <- as.list(1:N)
            }
        } ## end of check on whether its the final iteration and to output
    } # end of sample loop
    ##
    ## end loop on samples being processed
    ##

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
            refillMatrixList = refillMatrixList,
            infoCount = infoCount,
            afCount = afCount,
            hweCount = hweCount
        )
    )
}




get_nbreaks <- function(iteration, shuffleHaplotypeIterations, T) {
    nbreaks <- 0
    breaks <- NULL
    rStart <- NULL
    rEnd <- NULL
    if(is.na(match(iteration, shuffleHaplotypeIterations))==FALSE) {
        nSNP <- 100
        start <- round(nSNP / 2) + 1
        if (T >= (start + 100)) {
            breaks <- seq(start, T - nSNP, by = nSNP)
            nbreaks <- length(breaks) - 1
            rStart <- breaks - start + 1
            rEnd <- breaks - start + nSNP
            rEnd[length(rEnd)] <- T
        }
    }
    return(
        list(
            breaks = breaks,
            nbreaks = nbreaks,
            rStart = rStart,
            rEnd = rEnd
        )
    )
}

get_transMatRate <- function(method, sigmaCurrent) {
    if (method == "diploid") {
        x <- sigmaCurrent
        transMatRate_t <- rbind(x ** 2, x * (1 - x), (1 - x) ** 2)
    } else if (method == "pseudoHaploid") {
        transMatRate_t <- rbind(sigmaCurrent, 1 - sigmaCurrent)
    }
    return(transMatRate_t)
}

completeSampleIteration <- function(N,tempdir,chr,K,T,priorCurrent,eHapsCurrent,alphaMatCurrent,sigmaCurrent,maxDifferenceBetweenReads,whatToReturn,Jmax,aSW=NA,bSW=NA,highCovInLow,iteration,method,expRate,minRate,maxRate,niterations,splitReadIterations,shuffleHaplotypeIterations,nCores,L,nGen,alphaMatThreshold,emissionThreshold,gen,outputdir,environment,pseudoHaploidModel,outputHaplotypeProbabilities,switchModelIteration,regionName,restartIterations,refillIterations,hapSumCurrent,outputBlockSize,bundling_info, alleleCount, phase, samples_with_phase
    ) {

    print(paste(c("Start of iteration ",iteration," - ",date()),collapse=""))
    if(is.na(match(iteration,restartIterations))==FALSE  )
        print(paste0("Restart read probabilities ", date()))

    ##
    ## initial set up based on iteration
    ##
    if (niterations > 10)
        if(iteration <= 10)
            Jmax <- as.integer(iteration - 1)
    whatToReturn <- as.integer(1)
    ## final iteration - return everything - or splitting reads, need probabilities
    if (iteration == niterations)
        whatToReturn <- as.integer(2)
    ##
    ##
    ## set up switching iteration
    ##
    out <- get_nbreaks(iteration, shuffleHaplotypeIterations, T)
    nbreaks <- out$nbreaks
    breaks <- out$breaks
    rStart <- out$rStart
    rEnd <- out$rEnd
    ## set up splitting iteration
    nsplit <- 0
    if(is.na(match(iteration,splitReadIterations))==FALSE)
        nsplit <- 1
    ##
    x3 <- getSampleRange(N, nCores)
    ##
    hapSumCurrentL <- hapSumCurrent/rowSums(hapSumCurrent)
    ##
    ## transition matrix amount
    ##
    transMatRate_t <- get_transMatRate(
        method = method,
        sigmaCurrent = sigmaCurrent
    )

    ## ## TRANSPOSE-CLEAN
    alphaMatCurrent_t <- t(alphaMatCurrent)
    eHapsCurrent_t <- t(eHapsCurrent)

    if(environment=="server") {
        out2 <- mclapply(x3,mc.cores=nCores,tempdir=tempdir,chr=chr,K=K,T=T,priorCurrent=priorCurrent,eHapsCurrent_t=eHapsCurrent_t,alphaMatCurrent_t=alphaMatCurrent_t,sigmaCurrent=sigmaCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,whatToReturn=whatToReturn,Jmax=Jmax,highCovInLow=highCovInLow,iteration=iteration,method=method,nsplit=nsplit,expRate=expRate,minRate=minRate,maxRate=maxRate,gen=gen,outputdir=outputdir,pseudoHaploidModel=pseudoHaploidModel,outputHaplotypeProbabilities=outputHaplotypeProbabilities,switchModelIteration=switchModelIteration,regionName=regionName,restartIterations=restartIterations,refillIterations=refillIterations,hapSumCurrentL=hapSumCurrentL,outputBlockSize=outputBlockSize, bundling_info = bundling_info, transMatRate_t = transMatRate_t, x3 = x3, FUN=subset_of_complete_iteration, N = N, shuffleHaplotypeIterations = shuffleHaplotypeIterations, niterations = niterations, L = L, samples_with_phase = samples_with_phase, nbreaks = nbreaks, breaks = breaks)
    }
    if(environment=="cluster") {
        cl <- makeCluster(nCores, type = "FORK")
        out2 <- parLapply(cl, x3, fun=subset_of_complete_iteration,tempdir=tempdir,chr=chr,K=K,T=T,priorCurrent=priorCurrent,eHapsCurrent_t=eHapsCurrent_t,alphaMatCurrent_t = alphaMatCurrent_t,sigmaCurrent=sigmaCurrent,maxDifferenceBetweenReads=maxDifferenceBetweenReads,whatToReturn=whatToReturn,Jmax=Jmax,highCovInLow=highCovInLow,iteration=iteration,method=method,nsplit=nsplit,expRate=expRate,minRate=minRate,maxRate=maxRate,gen=gen,outputdir=outputdir,pseudoHaploidModel=pseudoHaploidModel,outputHaplotypeProbabilities=outputHaplotypeProbabilities,switchModelIteration=switchModelIteration,regionName=regionName,restartIterations=restartIterations,refillIterations=refillIterations,hapSumCurrentL=hapSumCurrentL,outputBlockSize=outputBlockSize, bundling_info = bundling_info, transMatRate_t= transMatRate_t, x3 = x3, N = N, shuffleHaplotypeIterations = shuffleHaplotypeIterations,niterations = niterations, L = L, samples_with_phase = samples_with_phase, nbreaks = nbreaks, breaks = breaks)
        stopCluster(cl)
    }

    ##
    ## done running samples
    ## first, check if there was an error, print it
    ##
    te <- sapply(out2, class) == "try-error"
    if (sum(te) > 0) {
        print(out2[[which(te)[1]]]) # print first error
        stop("An error occured during STITCH. The first such error is above")
    }


    ##
    ## update results
    ##
    out <- calculate_updates(out2 = out2, x3 = x3, K = K, T = T, N = N, nGen = nGen, expRate = expRate, minRate = minRate, maxRate = maxRate, emissionThreshold = emissionThreshold, alphaMatThreshold = alphaMatThreshold, L = L)
    sigmaSum <- out$sigmaSum
    priorSum <- out$priorSum
    alphaMatSum <- t(out$alphaMatSum_t) ## TRANSPOSE-CLEAN
    gammaSum <- t(out$gammaSum_t) ## TRANSPOSE-CLEAN
    hapSum <- t(out$hapSum_t) ## TRANSPOSE-CLEAN


    ##
    ## get other updates
    ##
    out <- calculate_misc_updates(N = N, nsplit = nsplit, nbreaks = nbreaks, x3 = x3, out2 = out2, fromMat = fromMat, K = K)
    restartMatrixList <- out$restartMatrixList
    readsTotal <- out$readsTotal
    readsSplit <- out$readsSplit
    fromMat <- out$fromMat


    ##
    ## calculate HWE, info, EAF if final iteration
    ##
    hwe <- NA
    info <- NA
    estimatedAlleleFrequency <- NA
    if(iteration == niterations) {
        infoCount <- array(0,c(T,2))
        hweCount <- array(0,c(T,3))
        afCount <- array(0,T)
        for(i in 1:length(out2)) {
            infoCount <- infoCount + out2[[i]]$infoCount
            hweCount <- hweCount + out2[[i]]$hweCount
            afCount <- afCount + out2[[i]]$afCount
        }
        thetaHat <- infoCount[,1] / 2 / N
        denom <- 2 * N * thetaHat * (1-thetaHat)
        info <- 1 - infoCount[,2] / denom
        estimatedAlleleFrequency <- afCount / N
        hwe <- generate_hwe_on_counts(hweCount, T, nCores)
    }
    ##
    ##
    ## look at switching
    ##
    if(nbreaks > 0) {
        out <- getBetterSwitchesSimple(fromMat=fromMat,nbreaks=nbreaks,K=K,T=T,eHapsFuture=gammaSum,alphaMatFuture=alphaMatSum,iteration=iteration, breaks = breaks, rStart = rStart, rEnd = rEnd)
        gammaSum <- out$eHapsFuture
        alphaMatSum <- out$alphaMatFuture
    }
    ##
    ## look at refilling - round 1
    ##
    if(is.na(match(iteration, refillIterations)) == FALSE) {
        print(paste("Iteration - ",iteration," - refill infrequently used haplotypes",sep=""))
        out=refillSimple(hapSum=hapSum,T=T,K=K,eHapsCurrent=gammaSum,N=N)
        eHapsCurrent=out$eHapsCurrent
        noise=0.2 # add in noise so things aren't exact
        eHapsCurrent=noise * array(runif(T*K),c(T,K)) + (1-noise)*eHapsCurrent
        priorSum=noise * rep(1/K,K) + (1-noise)*priorSum # restart these as well
        gammaSum=eHapsCurrent
        ## add in noise as well
    }
    ##
    ## look at splitting
    ##
    if(nsplit>0) {
        print(paste("split reads - average ",round(mean(readsSplit))," % ",round(100*mean(readsSplit)/mean(readsTotal),3),",",date(),sep=""))
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
            print(paste0("Printing out per-sample and mean estimates of correlation between provided genotypes from genfile and sample dosages, oriented so that the major allele has dosage 0 and the minor allele has dosage 1"))
        y <- round(currentR2, 3)
        print(paste0(
            "iteration=", iteration, ", sample(r2)=",
            paste(y[-length(y)], collapse = ", "),
            " - mean=", y[length(y)], ", ", date()
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
            print(paste0("Printing out per-sample phase switch error"))
        print(paste0(
            "iteration=", iteration, ", sample(pse %)=",
            paste(round(pse, 2) * 100, collapse = ", "), ", ", date()
        ))
    }
    ##
    ## also, if its a restart iteration, print out some stats
    ##
    if(is.na(match(iteration,restartIterations))==FALSE) {
        restartMatrix <- cbind(rep(1:N,sapply(restartMatrixList,nrow)),matrix(unlist(lapply(restartMatrixList,t)),ncol=5,byrow=TRUE))
        save(restartMatrix,file=paste0(outputdir,"RData/restartMatrix.iteration.",iteration,".",regionName,".RData"))
    }

    return(
        list(
            eHapsFuture = gammaSum,
            sigmaFuture = sigmaSum,
            priorFuture = priorSum,
            alphaMatFuture = alphaMatSum,
            hapSum = hapSum,
            info = info,
            estimatedAlleleFrequency = estimatedAlleleFrequency,
            hwe = hwe
        )
    )

}



calculate_misc_updates <- function(
N, nsplit, nbreaks, x3, out2, fromMat, K
) {
    fromMat <- array(0, c(nbreaks, K, K))
    restartMatrixList <- as.list(1:N)
    ## note - already defined fromMat before the function f above
    readsTotal <- array(0, N)
    readsSplit <- array(0, N)
    for(i in 1:length(x3)) {
        if(nbreaks > 0) {
            fromMat <- fromMat + out2[[i]][["fromMat"]]
        }
        if(nsplit > 0) {
            readsTotal <- readsTotal + out2[[i]][["readsTotal"]]
            readsSplit <- readsSplit + out2[[i]][["readsSplit"]]
        }
        restartMatrixList[x3[[i]][1]:x3[[i]][2]]=out2[[i]]$restartMatrixList[x3[[i]][1]:x3[[i]][2]]
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
    out2, x3, K, T, N,
    nGen, expRate, minRate, maxRate,
    emissionThreshold, alphaMatThreshold, L
) {
    priorSum <- array(0, K)
    alphaMatSum_t <- array(0, c(K, T - 1))
    gammaSum_t <- array(0, c(K, T, 2))
    hapSum_t <- array(0, c(K, T))
    ## sum and sum
    for(i in 1:length(x3)) {
        priorSum <- priorSum + out2[[i]]$priorSum
        alphaMatSum_t <- alphaMatSum_t + out2[[i]]$alphaMatSum_t
        gammaSum_t <- gammaSum_t + out2[[i]]$gammaSum_t
        hapSum_t <- hapSum_t + out2[[i]]$hapSum_t
    }

    ## normalize prior
    priorSum <- priorSum/sum(priorSum)
    ## note - it is possible in simulations not to have sampled the first SNP
    ## in which case there is no posterior information, so set to random
    if(is.na(priorSum[1]))
        priorSum <- rep(1/K,K)
    priorSum[priorSum<0.0001] <- 0.0001
    priorSum <- priorSum / sum(priorSum)
    ## bound above and below by a factor of 1000
    ## normalize sigma
    sigmaSum <- colSums(alphaMatSum_t) / N / 2
    ##sigmaSum=rowSums(alphaMatSum)/2  /nAlphaSum
    sigmaSum <- exp(-sigmaSum) # 100 is useless after start if its cant update
    ## NA values - set to expectation
    which <- is.na(sigmaSum)
    dl <- diff(L)
    ##if(is.na(aSW)==FALSE) dl=diff(L[aSW:bSW+1])
    sigmaSum[which] <- exp(-nGen * expRate / 100 / 1000000*dl)[which]
    ## morgans per SNP assuming T=100, 0.5 cM/Mb
    x1 <- exp(-nGen * minRate * dl/100/1000000) # lower
    x2 <- exp(-nGen * maxRate * dl/100/1000000) # upper
    ## recombination average rate
    ## we estiamte the compound parameter T * sigma_t
    ## so we need to use nGen to bound apppropriately
    sigmaSum[sigmaSum > x1] <- x1[sigmaSum > x1]
    sigmaSum[sigmaSum < x2] <- x2[sigmaSum < x2]

    gammaSum_t <- gammaSum_t[, , 1]/ gammaSum_t[, , 2]
    gammaSum_t[, colSums(is.na(gammaSum_t))>0] <- 0.5
    gammaSum_t[gammaSum_t > (1 - emissionThreshold)] <- (1 - emissionThreshold)
    gammaSum_t[gammaSum_t < emissionThreshold] <- emissionThreshold

    alphaMatSum_t<- alphaMatSum_t + 1
    alphaMatSum_t[alphaMatSum_t < alphaMatThreshold] <- 0
    ## hmm
    y <- alphaMatThreshold * colSums(alphaMatSum_t == 0)
    x <- colSums(alphaMatSum_t)
    for(k in 1:K)
        alphaMatSum_t[k, ] <- (1 - y) * alphaMatSum_t[k, ] / x
    alphaMatSum_t[alphaMatSum_t == 0] <- alphaMatThreshold

    ## now - cols, meaning SNPs, want to
    t1 <- colSums(is.na(alphaMatSum_t)) > 0
    if ( sum(t1) > 0)
        print("WARNING - alphaMat update contains NA")
    alphaMatSum_t[, t1] <- 0.25

    return(
        list(
            sigmaSum = sigmaSum,
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
restartSample <- function(sampleReads, srp, pRgivenH1, pRgivenH2,fbsoL,T,eHapsCurrent_t,sigmaCurrent,alphaMatCurrent_t,Jmax,priorCurrent,maxDifferenceBetweenReads,outputdir,iteration,restartIterations, method) {
    if(is.na(match(iteration, restartIterations)) == TRUE)
        return(list(m = NULL, fbsoL = fbsoL))
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
        transMatRate_t <- get_transMatRate(
            "diploid",
            sigmaCurrent[TR[1]:(TR[2]-1)]
        )
        fbd <- forwardBackwardDiploid(
            sampleReads = sampleReads,
            nReads = as.integer(length(sampleReads)),
            pi = priorCurrent,
            eHaps_t = eHapsCurrent_t[, TR[1]:TR[2]],
            alphaMat_t = alphaMatCurrent_t[, TR[1]:(TR[2]-1)],
            transMatRate_t = transMatRate_t,
            Jmax = Jmax,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            whatToReturn = as.integer(1),
            suppressOutput = as.integer(1)
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
        if(wSX!=wEX & TR[2]!=T) # if they are different, and there's something left to do
        {
            a1=fbsoL[[1]]$gamma[(TR[2]+1):T,]
            a2=fbsoL[[2]]$gamma[(TR[2]+1):T,]
            fbsoL[[2]]$gamma[(TR[2]+1):T,]=a2
            fbsoL[[1]]$gamma[(TR[2]+1):T,]=a1
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







getBetterSwitchesSimple <- function(
    fromMat,
    nbreaks,
    K,
    T,
    eHapsFuture,
    alphaMatFuture,
    iteration,
    breaks,
    rStart,
    rEnd
) {
    ## greedy version
    ## choose, in order, one which means fewest moves
    switchOrder <- t(sapply(1:nbreaks,fromMat=fromMat,K=K,function(ib,fromMat,K) {
        fromMatL=fromMat[ib,,]
        ## columns are to
        ## rows are from
        ## ie entry fromMat[10,2,3] is from breaks[10] hap 2to breaks[10+1] hap 3
        orderL=array(0,K)
        fromL=1:K
        toL=1:K
        fromR=array(0,K)
        toR=array(0,K)
        for(j in 1:(K-1)) {
            ## choose one which would give best diag
            f=(1:ncol(fromMatL))[rowSums(fromMatL==max(fromMatL))==1]
            t=(1:ncol(fromMatL))[colSums(fromMatL==max(fromMatL))==1]
            fromR[j]=fromL[f]
            toR[j]=toL[t]
            ## now - remove from consideration
            fromMatL=fromMatL[-f,-t]
            fromL=fromL[-f]
            toL=toL[-t]
        }
        fromR[K]=fromL
        toR[K]=toL
        ## re-order by fromR
        return(toR[order(fromR)])
    }))
    whichIsBest <- as.integer(apply(switchOrder,1,function(x) sum(x==(1:K)))==K)
    ##
    ## determine the order of subsequent regions
    ##
    ## now - do the unravelling
    tempMat <- array(0,c(nbreaks+1,K))
    tempMat[1, ] <- 1:K
    currentState <- 1:K
    nextState <- 1:K
    for (iBreak in 1:nbreaks) {
        ## how this works - first, start in your current state
        ## next, choose new state from switchOrder. copy that state, then remember where to go next
        currentState=nextState
        for(k in 1:K) {
            ## first - the current value is wherever you are
            tempMat[iBreak+1,k]=switchOrder[iBreak,currentState[k]]
            nextState[k]=tempMat[iBreak+1,k]
            ## now - record the number for where to start next time
        }
    }
    ##
    ## do the shuffling
    ##
    print(paste0("Shuffle haplotypes - Iteration ",iteration," - change ",sum(whichIsBest!=1)," out of ",nbreaks," N=100 SNP intervals"))
    switchCur=1:K
    for(iBreak in 1:(nbreaks+1)) {
      permL=tempMat[iBreak,]
      eHapsFuture[rStart[iBreak]:rEnd[iBreak],] =eHapsFuture[rStart[iBreak]:rEnd[iBreak],permL]
    }
    for(iBreak in 1:nbreaks) {
        if(whichIsBest[iBreak]==0) {
            ## add in some noise around break
            ## add in randomness around where I decide to switch
            nSNPL=20
            y=(rEnd[iBreak]-nSNPL):(rEnd[iBreak]+nSNPL)
            eHapsFuture[y,]=runif((2*nSNPL+1)*K)
            alphaMatFuture[y,]=matrix(1/K,ncol=K,nrow=2*nSNPL+1)
        }
    }
    return(
        list(
            eHapsFuture = eHapsFuture,
            alphaMatFuture = alphaMatFuture,
            switchOrder = switchOrder,
            whichIsBest = whichIsBest,
            tempMat = tempMat
        )
    )
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

  # maximum-ish value
  mid <- floor((nA * nB) / ((nA + nB)))
  if (is.na(mid))
    mid <- 0
  # make odd if it needs to be
  if ((nA%%2) != (mid%%2))
    mid <- mid + 1

  # determine a mid point
  probs <- array(0, max_het + 1) # note - this is 0-based
  probs[mid + 1] <- 1

  # use recurrence relation - going down
  n_het <- mid
  n_hom_alt <- (nBB * 2 + nAB - n_het) / 2
  n_hom_ref <- (n - n_het - n_hom_alt)
  if ((mid - 2) >= min_het) {
    for(het in seq(mid - 2, min_het, -2)) {
      probs[het + 1] <- probs[het + 3] *
        n_het * (n_het-1) / (4 * (n_hom_ref + 1) * (n_hom_alt + 1))
      n_het <- n_het - 2
      n_hom_ref <- n_hom_ref + 1
      n_hom_alt <- n_hom_alt + 1
    }
  }

  # use recurrence relationship - going up
  n_het <- mid
  n_hom_alt <- (nBB * 2 + nAB - n_het) / 2
  n_hom_ref <- (n - n_het - n_hom_alt)
  if ((mid + 2) <= max_het) {
    for(het in seq(mid + 2, max_het, 2)) {
      probs[het + 1] <- probs[het - 1] *
        (4 * n_hom_ref * n_hom_alt) / ( (n_het + 2) * (n_het + 1))
      n_het <- n_het + 2
      n_hom_ref <- n_hom_ref - 1
      n_hom_alt <- n_hom_alt - 1
    }
  }

  all_probs <- probs / sum(probs)
  p <- min(1, sum(all_probs[all_probs <= all_probs[nAB + 1]]))
  return(p)
}




generate_hwe_on_counts <- function(
  hweCount,
  T,
  nCores
) {
    ## few SNPs - make at least 5 per entry
    x <- seq(1, T, by = ceiling(T/nCores))
    if(ceiling(T/nCores)<10)
        x <- seq(1, T, by = ceiling(T/10))
    x2=c(x[-1]-1, T)
    x3=as.list(1:nCores)
    for(i in 1:length(x2))
        x3[[i]]= c(x[i],x2[i])
    x3=x3[unlist(lapply(x3,length))>1] # make sure each one has at least two to load
    out <- mclapply(x3,hweCount=hweCount,function(i, hweCount) {
        return(apply(hweCount[i[1]:i[2], , drop = FALSE],1,calculate_hwe_p))
    })
    out <- unlist(out)
    return(out)
}



# function to get phred scaled likelihood and genotype
# from input RData file
get_pl_and_rd_from_tempdir <- function(
  sampleReads,
  T
) {
  # get genotype likelihoods
  x=unlist(sapply(sampleReads,function(x) x[[4]]))+1
  bq=unlist(sapply(sampleReads,function(x) x[[3]]))
  s=as.integer(bq>0)
  bq=abs(bq)
  e=10^(-bq/10)
    p=cbind(1-e,e*(1/3))
    pA=p[cbind(1:nrow(p),s+1)]
    pB=p[cbind(1:nrow(p),2-s)]
    # hmm, just do for loop? not ideal, but OK?
    g=array(1,c(T,3))
    rd=array(0,c(T,2))
    for(i in 1:length(x)) {
      a=x[i]
      #b=y[i]
      b=pA[i]/2
      c=pB[i]/2
      g[a,]=g[a,] * c(b+b, b+c, c+c) # new
      z <- which.max(c(pA[i], pB[i]))
      rd[a,z] = rd[a,z] + 1
    }
    g=g/rowSums(g) # divide by max
    a=g[,1]
    w=a<g[,2]; a[w]=g[w,2]
    w=a<g[,3]; a[w]=g[w,3]
    g=g/a
    # 0/0:2,0:2:3:0,3,45
    # turn into PLs
    # ex: 1/1:0,1:1:3:32,3,0
    pl=round(10 * -log10(g) )
    return(list(pl = pl, rd = rd))
}


# if want to output genotype probabilities in VCF format
outputInputInVCFFunction <- function(
  outputdir,
  pos,
  T,
  tempdir,
  N,
  nCores,
  regionName,
  environment,
  sampleNames,
  outputBlockSize,
  bundling_info
) {
  #
  # need to loop over samples, very similar to normal output
  # write out chunks
  # except this time, want genotype, allelic depth, pl
  # gt,ad,pl
  #
  #
  #
  #
  # 1 - build the data
  #
  #
  #
  #
  x3 <- getSampleRange(N, nCores)
  f <- function(sampleRange, N, outputBlockSize, tempdir, regionName, T, bundling_info) {
    i_core <- match(sampleRange[1], sapply(x3, function(x) x[[1]]))
    list_of_vcf_columns_to_out <- as.list(1:N)
    outputBlockRange <- getOutputBlockRange(sampleRange,outputBlockSize )
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

      o <- get_pl_and_rd_from_tempdir(sampleReads, T)
      pl <- o$pl
      rd <- o$rd
      list_of_vcf_columns_to_out[[iSample]] <- make_column_of_vcf_from_pl_and_rd(pl, rd)
      # if appropriate - write to disk!
      iBlock <- match(iSample, outputBlockRange)
      if (iBlock > 1 & is.na(iBlock)==FALSE) {
        iSample_list <- (outputBlockRange[iBlock-1]+1):outputBlockRange[iBlock]
        write_block_of_vcf(i_core, iBlock, list_of_vcf_columns_to_out, iSample_list = iSample_list, outputdir, regionName)
        list_of_vcf_columns_to_out <- as.list(1:N)
      }
    }
  }
  #
  # loop over samples - write to disk!
  #
  if(environment=="server")
  {
    out2=mclapply(x3,mc.cores=nCores,FUN=f, N = N, outputBlockSize = outputBlockSize, tempdir = tempdir, regionName = regionName, T = T, bundling_info = bundling_info)
  }
  if(environment=="cluster")
  {
    cl = makeCluster(nCores, type = "FORK")
    out2 = parLapply(cl, x3, fun=f,N = N, outputBlockSize = outputBlockSize, tempdir = tempdir, regionName = regionName, T = T, bundling_info = bundling_info)
    stopCluster(cl)
  }
  error_check <- sapply(out2, class) == "try-error"
  if (sum(error_check) > 0) {
    print(out2[[which(error_check)[1]]])
    stop("There has been an error generating the input. Please see error message above")
  }
  #
  #
  #
  #
  # step 2 - build the VCF itself
  #
  #
  #
  #
  # build header
  # build left side
  # stitch everything together
  print(paste0("Build vcf from input - ",date()))
  output_vcf_header <- paste0(outputdir,"stitch.input.",regionName,".header.vcf")
  output_vcf_left <- paste0(outputdir,"stitch.input.",regionName,".left.vcf.gz")
  output_vcf <- paste0(outputdir,"stitch.input.",regionName,".vcf.gz")
  #
  # do header lines
  #
  header <- paste0(
    '##fileformat=VCFv4.0\n',
    '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n'
  )
  header2 <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(sampleNames, collapse = "\t", sep="\t"), sep="\t")
  #
  # add output
  #
  cat(header, header2, "\n", sep="", file = output_vcf_header)
  #
  # build first part of file
  #
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
  #
  # get number of cores and blocks
  #
  x3 <- getSampleRange(N, nCores)
  files_to_paste <- lapply(1:nCores, function(i_core) {
    nBlocks <- length(getOutputBlockRange(
      sampleRange = x3[[i_core]],
      outputBlockSize
    ))
    return(paste0('<(gunzip -c ',outputdir, "vcf.piece.",i_core,".",2:nBlocks,".",regionName,".txt.gz", ' | cat ) '))
})
  files_to_paste <- paste(unlist(files_to_paste), collapse = " ")
  #
  # merge together
  #
  command <- paste0(
    'bash -c "',
    'paste -d ',shQuote("\t"),' ',
    '<(gunzip -c ',output_vcf_left,' | cat ) ',
    files_to_paste, " | ",
    " cat ",output_vcf_header," - | bgzip -c > ",
      output_vcf
    ,'"'
  )
  #print(command)
  system(
    command
  )
  system(paste0("rm ",output_vcf_header))
  system(paste0("rm ",output_vcf_left))
  #
  # remove
  #
  x3 <- getSampleRange(N, nCores)
  files_to_remove <- unlist(lapply(1:nCores, function(i_core) {
    nBlocks <- length(getOutputBlockRange(
      sampleRange = x3[[i_core]],
      outputBlockSize
    ))
    return(paste0(outputdir, "vcf.piece.",i_core,".",2:nBlocks,".",regionName,".txt.gz"))
  }))
  system(paste0("rm ",paste0(files_to_remove, collapse=" ")))
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
        y <- 2 * out[[isample]]
        ##print(paste0("range(x) = ", range(x, na.rm=TRUE)))
        ##print(paste0("range(y) = ", range(y, na.rm=TRUE)))
        ##x[to_flip] <- 2 - x[to_flip]
        ##y[to_flip] <- 2 - y[to_flip]
        store[isample] <- cor(x, y, use = "complete.obs") ** 2
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


### plot HWE, info, MAF, coloured by discrepency
plotMetricsForPostImputationQC=function(iSample,highCovList,gen,gen_imp,alleleCount,chr,L,estimatedAlleleFrequency,info,outputdir,colour=TRUE,hwe,regionName) {
  if(is.null(iSample)==FALSE)
  {
    m1A=cbind(gen[,iSample]/2,gen_imp[, iSample])
    m1A[alleleCount[,3]>0.5,]=1-m1A[alleleCount[,3]>0.5,]
    dist=abs(m1A[,1]-m1A[,2])
    dist[dist>0.1]=0.1
  } else {
    dist=rep(0,T)
  }
  colfunc <- colorRampPalette(c("blue", "red"))
  col=colfunc(11)[round(dist,2)*100+1]
  ### plot
  jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQC.",regionName,".sample",iSample,".jpg")),height=1000,width=3000,qual=100)
  par(mfrow=c(1,3))
  plot(info,log10(hwe) ,cex.lab=3 ,col=col)
  plot(info,estimatedAlleleFrequency,cex.lab=3,col=col)
  plot(estimatedAlleleFrequency,log10(hwe) ,cex.lab=3,col=col )
  dev.off()
  ### now plot along the chromosome as well
  jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQCChromosomeWide.",regionName,".sample",iSample,".jpg")),height=3000,width=3000,qual=100)
  par(mfrow=c(3,1))
  plot(L,log10(hwe) ,cex.lab=3 ,col=col)
  plot(L,info,cex.lab=3,col=col)
  plot(L,estimatedAlleleFrequency,cex.lab=3,col=col )
  dev.off()
}

### plot hapSum along genome
plotHapSum <- function(
    outname,
    L,
    K,
    hapSum,
    T,
    N
) {
    jpeg(outname, height = 2000, width = 10000, qual = 100)
    colStore <- rep(c("black","red","green","blue"), ceiling(K/4))
    sum <- array(0, T)
    xlim <- range(L)
    ylim <- c(0,1)
    plot(x = L[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE)
    x <- c(L[1], L, L[length(L):1])
    m <- array(0, c(T, K + 1))
    for(i in 1:K)
        m[, i + 1] <- m[, i] + hapSum[, i] / N
    for(i in K:1) {
        polygon(
            x = x, y = c(m[1,i], m[,i+1], m[T:1,i]),
            xlim = xlim, ylim = ylim, col = colStore[i]
        )
    }
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

get_reads_worse_than_50_50 <- function(sampleReads, eHapsCurrent_t, K) {
    whichReads <- sapply(sampleReads,function(x) x[[1]]) > 2
    m <- t(
        sapply(
            sampleReads[whichReads],
            split_function,
            eHapsCurrent_t = eHapsCurrent_t,
            K = K
        )
    )
    w=(1:length(sampleReads))[whichReads][(apply(m,1,which.max)==(K+1))]
    ## can't do first or last read
    ## w=w[w>1 & w<length(sampleReads)]
    return(w)
}



split_a_read <- function(
    sampleReads,
    read_to_split,
    gammaK_t,
    L,
    eHapsCurrent_t,
    K
) {

    did_split <- FALSE
    sampleRead <- sampleReads[[read_to_split]]

    ## get the likely for and after states
    ## by checking, a few reads up and downstream, what changes
    ## note - probably want to change this to per-base, not per-read!
    b <- c(
        max(1, read_to_split - 10),
        min(length(sampleReads), read_to_split + 10)
    )
    a <- c(
        sampleReads[[b[1]]][[2]] + 1,
        sampleReads[[b[2]]][[2]] + 1
    )
    change <- apply(gammaK_t[, a], 1, diff)
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
            sampleRead, 1, best_split, L
        )
        new_read_2 <- get_sampleRead_from_SNP_i_to_SNP_j(
            sampleRead, best_split + 1, sampleRead[[1]] + 1, L
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


get_sampleRead_from_SNP_i_to_SNP_j <- function(
    sampleRead,
    i,
    j,
    L
) {
    bq <- sampleRead[[3]][i:j]
    u <- sampleRead[[4]][i:j]
    ## central SNP in read
    cr <- u[getCentral(L[u + 1])]
    return(
        list(
            j - i,
            cr,
            matrix(bq, ncol = 1),
            matrix(u, ncol = 1)
        )
    )
}




findRecombinedReadsPerSample <- function(
    gammaK_t,
    eHapsCurrent_t,
    K,
    L,
    iSample,
    verbose=FALSE,
    sampleReads,
    tempdir,
    regionName
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
                eHapsCurrent_t,
                K = K
            )
            sampleReads <- out$sampleReads
            count <- count + as.integer(out$did_split)
        } # end of loop on reads
        sampleReads <- sampleReads[order(unlist(lapply(sampleReads,function(x) x[[2]])))]
        save(sampleReads, file = file_sampleReads(tempdir, iSample, regionName))
        if (verbose==TRUE)
            print(paste0("sample ",iSample," readsSplit ",count," readsTotal ",length(sampleReads)))
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
    pRgivenH1,pRgivenH2,fbsoL,pseudoHaploidModel,srp,delta, K
) {
      #
      # once both haplotypes are available - update omega immediately
      #
      # for each read, get new phasing
      # srp comes in with readProbs originally - doesnt change
      # 4 things here
      # 1 - pRgivenH1 - P(R | h1)
      # 2 - pRgivenH2 - P(R | h2)
      # note that
      # 3 - pRandH1 - P(R and H1) = 1/2 P(R| H1)
      # 4 - pRandH2 - P(R and H2) = 1/2 P(R| H2)
      ###
      ### MODEL 1 - divide one by other (y1/(y1+y2)) to get readProbs
      ###
      # readProbs is y1/(y1+y2) from above
      # then probabilities using
      #    eMatHap(iRead,k) = eMatHap(iRead,k) * readProbs(iRead) + (1/Kd)  * (1-readProbs(iRead));
      # then updating using multiply both by readProbs
      # y1=(fbsoL[[1]]$gamma[srp+1,] * fbsoL[[1]]$eMatHap)
      # y2=(fbsoL[[2]]$gamma[srp+1,] * fbsoL[[2]]$eMatHap)
      # readProbs=y1/(y1+y2)
      ###
      ### MODEL 3, 7 - try to predict which haplotype it came from
      ###
      if(pseudoHaploidModel==7 | pseudoHaploidModel==9) {
          ## original flavour
          #x=pRgivenH1/(pRgivenH1+pRgivenH2)
          #y1=(fbsoL[[1]]$eMatHap - (1-x)*pRgivenH2)/x
          #rowSums(fbsoL[[1]]$gamma[srp+1,] * y1)
          ## TRANSPOSE - write this in C++
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
      ###
      ### MODEL 8
      ###
      if(pseudoHaploidModel==8 | pseudoHaploidModel==10)
      {
        x=pRgivenH1/(pRgivenH1+pRgivenH2)
        # get eMatHap without the second part
        y1=fbsoL[[1]]$eMatHapOri # shouldn't change 1 or 2
        pRgivenH1=rowSums(fbsoL[[1]]$gamma[srp+1,] * y1)
        pRgivenH2=rowSums(fbsoL[[2]]$gamma[srp+1,] * y1)
        # bound above below
        pRgivenH1[pRgivenH1<0.001]=0.001
        pRgivenH2[pRgivenH2<0.001]=0.001
        pRgivenH1[pRgivenH1>0.999]=0.999
        pRgivenH2[pRgivenH2>0.999]=0.999
      }
      ###
      ### MODEL 4 -
      ###
      if(pseudoHaploidModel==4)
      {
        # model 3 -
        #x=0.5;
        #for(k=0; k<=K-1; k++)
        #  eMatHap(iRead,k) = x * eMatHap(iRead,k) + (1-x) * pRgivenH2(iRead);
        # update using proportions
        x=0.5;
        # get eMatHap without the second part
        y1=(fbsoL[[1]]$eMatHap - (1-x)*pRgivenH2)/x
        y2=(fbsoL[[2]]$eMatHap - (x)*pRgivenH1)/(1-x)
        pRgivenH1=rowSums(fbsoL[[1]]$gamma[srp+1,] * y1)
        pRgivenH2=rowSums(fbsoL[[2]]$gamma[srp+1,] * y2)
        # bound above below
        pRgivenH1[pRgivenH1<0.001]=0.001
        pRgivenH2[pRgivenH2<0.001]=0.001
        pRgivenH1[pRgivenH1>0.999]=0.999
        pRgivenH2[pRgivenH2>0.999]=0.999
      }
      ###
      ### MODEL 5 - complete proportioning
      ###
      if(pseudoHaploidModel==5)
      {
        # update using proportions
        x=pRgivenH1/(pRgivenH1+pRgivenH2)
        # get eMatHap without the second part
        y1=(fbsoL[[1]]$eMatHap - (1-x)*pRgivenH2)/x
        y2=(fbsoL[[2]]$eMatHap - (x)*pRgivenH1)/(1-x)
        # might have NA's - if so, replace one with other
        y=array(0,dim(y1))
        w=is.na(y1[,1])==FALSE; y[w,]=y1[w,]
        w=is.na(y2[,1])==FALSE; y[w,]=y2[w,]
        a1=rowSums(fbsoL[[1]]$gamma[srp+1,] * y)
        a2=rowSums(fbsoL[[2]]$gamma[srp+1,] * y)
        pRgivenH1=as.integer(runif(length(srp))<a1/(a1+a2))
        pRgivenH2=1-pRgivenH1
      }
  return(list(pRgivenH1=pRgivenH1,pRgivenH2=pRgivenH2,delta=delta))
}


## on final iteration, get imputed dosages
get_high_coverage_final_imputed_genotypes <- function(
    highCovInLow,
    tempdir,
    regionName,
    T
) {
    n_hc <- length(highCovInLow)
    out <- array(NA, c(T, n_hc))
    for(j in 1:n_hc) {
        load(file=file_dosages(tempdir, highCovInLow[j], regionName))
        out[, j] <- dosage
    }
    return(out)
}


get_max_gen_rapid <- function(x) {
  # assume matrix 3 columns >1 row
  z <- rep(1, nrow(x))
  y <- x[,1]
  for(i in 2:3) {
    w <- x[,i]>y
    z[w] <- i
    y[w] <- x[w,i]
  }
  return(cbind(1:nrow(x), z))
}

## write out a single samples worth of VCF entry
make_column_of_vcf <- function(
    gp,
    read_proportions
) {
    ## write out genotype, genotype likelihood, and dosage
    ##GT:GL:DS
    ##FORMAT=<ID=GT:,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
    ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Posterior probability of 0/0, 0/1, and 1/1">
    ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
    ## 1/1:0,0.054,0.946:1.946
    ## add one samples worth of info to a VCF
    z <- get_max_gen_rapid(gp)
    gt <- c("0/0","0/1","1/1")[z[,2]]
    gt[gp[z] < 0.9] <- "./."
    precision <- 3
    str <- paste0(
        gt,":",
        round(gp[,1], precision), ",",
        round(gp[,2], precision), ",",
        round(gp[,3], precision), ":",
        round(gp[,2] + 2 * gp[,3], precision)
    )
    if (is.null(read_proportions) == FALSE)
        str <- paste0(
            str, ":",
            round(read_proportions[, 1], precision), ",",
            round(read_proportions[, 2], precision), ",",
            round(read_proportions[, 3], precision), ",",
            round(read_proportions[, 4], precision)
        )
    return(str)
}


# write out a single samples worth of VCF entry
make_column_of_vcf_from_pl_and_rd <- function(
  pl,
  rd
) {
  # write out genotype, genotype likelihood, and dosage
  # gt, ad, pl
  # GT:AD:PL
  #chr19   3126522 .       G       A       7661.22 .       AC=464;AF=0.935;AN=496;BaseQRankSum=4.232;DP=333;Dels=0.00;FS=2.773;GC=39.90;HRun=0;HaplotypeScore=0.0738;InbreedingCoeff=-0.0083;MLEAC=465;MLEAF=0.938;MQ=42.81;MQ0=23;MQRankSum=-0.576;QD=28.69;ReadPosRankSum=-0.673      GT:AD:DP:GQ:PL       ./.     ./.     ./.     1/1:0,1:1:3:32,3,0      ./.     ./.
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


## write a block of samples to disk
## don't write header or row information
write_block_of_vcf <- function(
  i_core,
  iBlock,
  list_of_vcf_columns_to_out,
  iSample_list,
  outputdir,
  regionName
) {
    print(paste0("Write block of VCF for samples:",iSample_list[1],"-",iSample_list[length(iSample_list)], " - ", date()))
    m <- matrix(sapply(list_of_vcf_columns_to_out[iSample_list],I), ncol = length(iSample_list))
  write.table(
      m,
      file = gzfile(paste0(outputdir, "vcf.piece.",i_core,".",iBlock,".",regionName,".txt.gz")),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
  )
  return(NULL)
}


## specify here the output VCF from the process
get_output_vcf <- function(
    vcf_output_name,
    outputdir,
    regionName
) {
    if (is.null(vcf_output_name)) {
        return(file.path(
            outputdir,
            paste0("stitch.", regionName, ".vcf.gz")
        ))
    } else {
        if (basename(vcf_output_name) == vcf_output_name) {
            return(file.path(
                outputdir,
                vcf_output_name
            ))
        } else {
            return(vcf_output_name)
        }
    }
}

# build the header
# build the left part of the VCF
# paste together the interim VCF files
# gzip the whole thing together
write_vcf_after_EM <- function(
  vcf_output_name,
  outputdir,
  regionName,
  sampleNames,
  tempdir,
  nCores,
  info,
  hwe,
  estimatedAlleleFrequency,
  pos,
  N,
  outputBlockSize,
  reference_panel_SNPs,
  method
) {
    ## set up file names
    print(paste0("Build final vcf - ",date()))
    output_vcf <- get_output_vcf(vcf_output_name,  outputdir, regionName)
    output_vcf_header <- paste0(output_vcf, ".header.gz")
    output_vcf_left <- paste0(output_vcf, ".left.gz")
    ##
    ## do header lines
    ##
    header_line_if_ref_panel_snps <- NULL
    if (sum(reference_panel_SNPs) > 0)
        header_line_if_ref_panel_snps <- '##INFO=<ID=REF_PANEL,Number=.,Type=Integer,Description="Whether a SNP was (1) or was not (0) found in the reference panel during imputation">\n'
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        '##INFO=<ID=INFO_SCORE,Number=.,Type=Float,Description="Info score">\n',
        '##INFO=<ID=EAF,Number=.,Type=Float,Description="Estimated allele frequency">\n',
        '##INFO=<ID=HWE,Number=.,Type=Float,Description="Hardy-Weinberg p-value">\n',
        header_line_if_ref_panel_snps,
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">\n',
        '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Posterior genotype probability of 0/0, 0/1, and 1/1">\n',
        '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">\n'
    )
    if (method == "pseudoHaploid" && 1 == 0) ## disable for now
        header <- paste0(
            header,
            '##FORMAT=<ID=ER1,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 1">\n',
            '##FORMAT=<ID=EA1,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 1">\n',
            '##FORMAT=<ID=ER2,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 2">\n',
            '##FORMAT=<ID=EA2,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 2">\n'
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
    INFO <- paste0(
        "EAF=", round(estimatedAlleleFrequency, 5), ";",
        "INFO_SCORE=", round(info, 5), ";",
        "HWE=", hwe
    )
    if (sum(reference_panel_SNPs) > 0)
        INFO <- paste0(INFO, ";REF_PANEL=", as.integer(reference_panel_SNPs))
    write.table(
        matrix(paste(pos[,1], pos[,2], ".", pos[,3], pos[,4], ".", "PASS", INFO, "GT:GP:DS", sep="\t"), ncol=1),
        file = gzfile(output_vcf_left),
        row.names = FALSE,
        col.names = FALSE,
        sep = "",
        quote = FALSE
    )
    ##
    ## get number of cores and blocks
    ##
    x3 <- getSampleRange(N, nCores)
    files_to_paste <- lapply(1:length(x3), function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(paste0('<(gunzip -c ',outputdir, "vcf.piece.",i_core,".",2:nBlocks,".",regionName,".txt.gz", ' | cat ) '))
    })
    files_to_paste <- paste(unlist(files_to_paste), collapse = " ")
    ##
    ## merge together
    ##
    command <- paste0(
        'bash -c "',
        'paste -d ',shQuote("\t"),' ',
        '<(gunzip -c ',output_vcf_left,' | cat ) ',
        files_to_paste, " | ",
        " cat ",output_vcf_header," - | bgzip -c > ",
        output_vcf
       ,'"'
    )
    system(
        command
    )
    system(paste0("rm ",output_vcf_header))
    system(paste0("rm ",output_vcf_left))
    ##
    ## remove
    ##
    x3 <- getSampleRange(N, nCores)
    files_to_remove <- unlist(lapply(1:length(x3), function(i_core) {
        nBlocks <- length(getOutputBlockRange(
            sampleRange = x3[[i_core]],
            outputBlockSize
        ))
        return(paste0(outputdir, "vcf.piece.",i_core,".",2:nBlocks,".",regionName,".txt.gz"))
    }))
    system(paste0("rm ",paste0(files_to_remove, collapse=" ")))
    return(NULL)
}



estimate_read_proportions <- function(
    sampleReads,
    pRgivenH1,
    pRgivenH2,
    T
) {
    output <- array(0, c(T, 4))
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
