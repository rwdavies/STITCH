* v1.6.7
	* Fix bugs, change seqlib library
* v1.6.6
	* Export more functions for library use
* v1.6.5
	* Be able to use bx tag to link reads
	* Export yet more functions for library use
* v1.6.4
	* Export more functions for library use
* v1.6.3
	* Misc minor bug fixes
* v1.6.2
	* Change internal haplotype format (rhb)
	* Speed up reference haplotype file loading from disk
* v1.6.1
	* Speed up reference haplotype EM considerably
* v1.6.0
	* Introduce new variable S as the number of sets of parameters results are averaged over
* v1.5.7
	* Push CC through to htslib to robustify against compilation issues
* v1.5.6
	* Add reference_shuffleHaplotypeIterations
* v1.5.5
	* Export functions so STITCH can serve as a library
* v1.5.4
	* Speed up analysis of very large BAMs
	* Can output haplotype dosage
	* Change release mechanism
* v1.5.3
	* Significant speedups and some RAM improvements
	* Minor changes to heuristics
* v1.5.2
	* Capability to substantially decrease RAM usage when running with a large number of samples
	* Better error message when imputing few / no SNPs
	* Able to impute when only 2 SNPs, including 1 SNP in central region
* v1.5.1
	* Add method diploid-inbred
	* Improve capabilities of heuristic looking for shuffled haplotypes
	* Fix C++ bug with HWE
	* Internal improvements with testing flexibility and packaging efficiency
* v1.5.0
	* Add option to write to bgen output
	* Make writing bgzipped VCFs more efficient
	* Various other minor speedups
* v1.4.2
	* Introduce new variable maxEmissionMatrixDifference to better control magnitude of differences in emission matrix, and reduce underflow likelihood when using gridWindowSize
	* Apply downsampling to gridded sampleReads to reduce underflow likelihood
	* Allow disabling of plots in cases with poor x11 support
* v1.4.1
	* Fix bug when using gridWindowSize and buffer
* v1.4.0
	* Added option gridWindowSize to perform HMM on blocks of physical size gridWindowSize, which can dramatically speed up analysis for low coverage samples (see Benchmark section, examples)
	* Added benchmark section to GitHub
	* Misc small speedups, internal variable re-naming
* v1.3.7
	* Misc changes to reduce RAM usage
	* Standardize messaging and change to stderr vs stdout
	* Re-enable printing of r2 oriented for major allele during progression
	* More errors and checking for reference legend file
* v1.3.6
	* Allow very strict imputation from reference panels with niterations=1
	* Fix error message printing for files missing RG bam header tag
	* Fix bug when parsing read with cigarString *
	* Require SM entry in RG tag
* v1.3.5
	* Throw an error if there are problems initializing directories
	* Better error messages for missing files
	* Change downsample behaviour to require no more than downsampleToCov at all SNPs for a sample and not just lead SNPs in a read
* v1.3.4
	* Speed up BAM/CRAM input processing by reformatting reads on the fly
* v1.3.3
	* Fix bug that occured when the number of samples is similar to or less than nCores
* v1.3.2
	* Increase likelihood of succesfull compilation by chaging Makevars to compile SeqLib and htslib with the same configuration as R
* v1.3.1
	* Move SeqLib installation into Makevars to harmonize installation configuration with R
	* Fix generateInputOnly to only generate input then stop
	* Fix bug that arose when a sample has reads but no reads meeting the mapping quality or isize threshold requirements
* v1.3.0
	* Use SeqLib instead of Rsamtools to get read inforamtion. This speeds up analysis of BAM files and significantly speeds up analysis of CRAM files
	* Use SeqLib instead of samtools to get sample names
* v1.2.9
	* Add command line wrapper for STITCH to facilitate running from the command line
	* Reduce RAM when working on large regions of the genome (e.g. chromosomes) by loading raw data 1 Mbp at a time
	* Change examples script to only showcase examples and not installing dependencies
* v1.2.8
	* Change output to bgzipped VCF from gzip. Require bgzip to be in PATH
* v1.2.7
	* Speed up C++ functions by re-orienting internal matrices to better use column-major order
	* Fix potential compilation specific bug requiring specification of C++ header iomanip
	* Fix bug about loading reference haplotypes from X chromosome for male reference samples
	* More tests
* v1.2.6
	* Fix bugs that manifest when the number of SNPs is very small (in the tens)
	* More tests
	* Increase efficiency of code, particularly C++ code
* v1.2.5
	* Require samtools to be in PATH
	* More thorough and better tested input validation
* v1.2.4
	* Enable C++11 compilation
	* Fix bug where sample name from bam header was being grabbed from any line with @RG in it and not specifically lines starting with @RG
* v1.2.3
	* Fix bug where the central SNP in a read was random from SNPs in read and not the central SNP by position among SNPs in the read as it ought to have been
	* Fix bug where the final SNP from posfile wasn't being loaded from the sample BAMs and as a result not being imputed
	* Fix bug where reads split into 3+ pieces were not being properly handled (e.g. long reads where sections map to multiple locations)
	* Faster internal handling of cigar string
* v1.2.2
	* Reduce RAM footprint when loading reference haplotypes
	* Crash early if rsync is not in PATH
	* Added option to override default VCF output name. See vcf_output_name
	* Added unit tests under testthat framework
* v1.2.1
	* Change internal system calls to reduce RAM usage
	* Fix bug passing through variable to subfunction
* v1.2.0
	* Can use reference panels in IMPUTE2 format. See example script and reference_* variables
	* Can bundle together inputs to facilitate imputation of very large N. See inputBundleBlockSize
	* Example human data provided to showcase STITCH functionality. See examples script
* v1.1.4
	* Can work off CRAM files or BAM files. To use CRAM files, see cramlist and reference variables, or see examples script
	* Changed GL as genotype likelihood to GP as genotype posterior probability in output VCF
* v1.1.3
	* Changed R example script to work on provided example data
	* Changed default downsampleToCov to 50 to reduce likelihood of overflow at high coverage SNPs
	* Miscellaneous small fixes to BAM conversion script to better handle samples with very few reads
	* Changed default region where reads from BAM are loaded (chrStart and chrEnd) to NA, to be inferred from posfile and the region to be imputed, rather than to grab reads from the whole chromosome
	* Fixed ability to use high coverage validation samples (genfile) when using generateInputOnly and regenerateInput
* v1.1.2
	* Fixed typo in header of outputted VCF
* v1.1.1
	* Fixed bug where samples with no reads on a chromosome gave an old input format
* v1.1.0
	* Changed default output to VCF (added option outputBlockSize to control how it is written)
        * Removed two package dependencies
        * Remove ability to write output to .gen format (remove outputGenFormat)
        * Add change log to README and miscellaneous other changes
* v1.0.1
	* Added ability to use soft clipped bases
* v1.0.0
	* Version used for paper 

