STITCH - Sequencing To Imputation Through Constructing Haplotypes
=================================================================
**__Current Version: 1.6.7__**
Release date: Dec 8,  2022

![Build Status](https://github.com/rwdavies/STITCH/workflows/CI/badge.svg)


Changes in latest version

1. Fix bugs, change seqlib library

For details of past changes please see [CHANGELOG](CHANGELOG.md).

STITCH is an R program for reference panel free, read aware, low coverage sequencing genotype imputation. STITCH runs on a set of samples with sequencing reads in BAM format, as well as a list of positions to genotype, and outputs imputed genotypes in VCF format. 

For the old website, please see https://www.well.ox.ac.uk/~rwdavies/stitch.html

# Table of contents

1. [Installation](#paragraph-installation)
    1. [github](#paragraph-installation-github)
    2. [conda](#paragraph-installation-conda)
    3. [missing libraries](#paragraph-installation-missing-library)
2. [Quick start run](#paragraph-quickstartrun)
3. [Interactive start](#paragraph-interactive-start)
4. [Options and help](#paragraph-optionsandhelp)
5. [Benchmarks](#paragraph-benchmarks)
6. [Examples](#paragraph-examples)
7. [License](#paragraph-license)
8. [Citation](#paragraph-citation)
9. [Testing](#paragraph-testing)
10. [Bug reports](#paragraph-bugreports)
11. [Output format](#paragraph-output-format)
12. [What method to run](#paragraph-what-method)
13. [What method to run](#paragraph-time-ram-memory)
14. [How to choose K and nGen](#paragraph-k-ngen)
15. [About plotting](#paragraph-plots)
16. [About reference panels](#paragraph-reference-panels)

## Installation <a name="paragraph-installation"></a>

STITCH is available to download either through this github repository, or through conda.

### github <a name="paragraph-installation-github"></a>

A simple way to ensure dependencies are installed, and to install a release of STITCH is as follows. First, install R. Then, do the following 
```
git clone --recursive https://github.com/rwdavies/STITCH.git
cd STITCH
./scripts/install-dependencies.sh
cd releases
wget https://github.com/rwdavies/stitch/releases/download/1.6.6/STITCH_1.6.6.tar.gz ## or curl -O
R CMD INSTALL STITCH_1.6.6.tar.gz
```

You can confirm the installation worked using the quick start run below.

To install the latest development code in the repository, use `./scripts/build-and-install.sh`. To install alternative releases, either download other releases from Github like done above, or use the historical `releases` directory.

Note that STITCH as run in the original paper used version 3 of R. However STITCH should work fine with either version 3 or version 4 of R. If you have dependency problems, you can easier post an issue on github, or try the conda installation below.

### conda <a name="paragraph-installation-conda"></a>

STITCH (as r-stitch) can be installed using [conda](https://conda.io/miniconda.html). Full tutorials can be found elsewhere, but briefly, something like this should work
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install r-stitch -c defaults -c bioconda -c conda-forge
source activate
R -e 'library("STITCH")'
```
Note that currently the command like `STITCH.R` is not included with the bioconda installation, so from the command line, you can either run something like `R -e 'library("STITCH"); STITCH(chr="chr19", bamlist="bamlist.txt", posfile="pos.txt", genfile="gen.txt", outputdir="./", K=4, nGen=100, nCores=1)'`, or clone the repo to get `STITCH.R`.

You can confirm the installation worked using the quick start run below.


### missing libraries <a name="paragraph-installation-missing-libraries"></a>

If you experience a problem with installation, you can either try conda above. Alternatively, if you see an error similar to ```error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory```, then either ask your system administrator to install gmp, mpfr and mpc for you, or try running the following before R CMD INSTALL
```
./scripts/install-package-dependencies.sh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`pwd`/install/lib/
```

If you see an error similar to ```configure: error: liblzma not found, please install lzma```, then either ask your system administrator to install lzma or xz for you, or try running the following before R CMD INSTALL
```
./scripts/install-xz.sh
echo "CPPFLAGS += -I`pwd`/install/include" >> ~/.R/Makevars
echo "LDFLAGS += -L`pwd`/install/lib" >> ~/.R/Makevars
```

If you're on Mac you may see an error similar to ```ld: library not found for -lquadmath```, which is related to STITCH C++ compilation using Rcpp. This can be fixed by updating gfortran using a method such as [this](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/). If you experience other compilation issues, please raise an issue. To experiment with configuration options during compilation, you can edit ```STITCH/src/Makevars``` then build a package and install using ```./scripts/build-and-install.sh``` or test using ```./scripts/test-unit.sh```.



## Quick start run <a name="paragraph-quickstartrun"></a>

A quick test on real data can be performed using 
```
# test on CFW mouse data
wget https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
# or curl -O https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
tar -xzvf STITCH_example_2016_05_10.tgz
./STITCH.R --chr=chr19 --bamlist=bamlist.txt --posfile=pos.txt --genfile=gen.txt --outputdir=./ --K=4 --nGen=100 --nCores=1
# if this works the file stitch.chr19.vcf.gz will be created
```




## Interactive start <a name="paragraph-interactive-start"></a>

It is recommended you follow the instructions above, specifically `./scripts/install-dependencies.sh` to install the dependencies, but if you run into problems, or want to install in a more manual fashion, the below should work

1. Install R if not already installed.
2. Install R dependencies parallel, Rcpp and RcppArmadillo from CRAN (using the "install.packages" option within R)
3. Install [bgzip](http://www.htslib.org/) and make it available to your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)). This can be done using a system installation, or doing a local installation and either modifying the PATH variable using code like ```export PATH=/path/to/dir-with-bgzip-binary/:$PATH```, or through R, doing something like ```Sys.setenv( PATH = paste0("/path/to/dir-with-bgzip-binary/:", Sys.getenv("PATH")))```. You'll know samtools is available if you run something like ```system("which bgzip")``` in R and get the path to bgzip
4. Install STITCH. First, download the latest STITCH tar.gz from the releases folder above (or more ideally the releases section of the github page). Second, install by opening R and using install.packages, giving install.packages the path to the downloaded STITCH tar.gz. This should install SeqLib automatically as well.
5. Download example dataset [STITCH_example_2016_05_10.tgz](https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz).
6. Run STITCH. Open R, change your working directory using setwd() to the directory where the example tar.gz was unzipped, and then run ```STITCH(tempdir = tempdir(), chr = "chr19", bamlist = "bamlist.txt", posfile = "pos.txt", genfile = "gen.txt", outputdir = paste0(getwd(), "/"), K = 4, nGen = 100, nCores = 1)```. Once complete, a VCF should appear in the current working directory named stitch.chr19.vcf.gz


## Options and help <a name="paragraph-optionsandhelp"></a>

For a full list of options, in R, query ```?STITCH```, or from the command line, ```STITCH --help```.

For a brief writeup of commonly used variables, see [Options.md](Options.md). To pass vectors using the command line, do something like ```STITCH.R --refillIterations='c(3,40)'``` or ```STITCH.R --reference_populations='c("CEU","GBR")'```.

For help about errors, see the bug reports section.

## Benchmarks <a name="paragraph-benchmarks"></a>

One can see some speed benchmarks in [benchmarks/summarize_benchmarking.md](benchmarks/summarize_benchmarking.md)

## Examples <a name="paragraph-examples"></a>

In the examples directory, there is a script which contains examples using real mouse and human data. One can either run this interactively in R, or run all examples using ```./examples/example.R```.

## License <a name="paragraph-license"></a>

STITCH and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Citation <a name="paragraph-citation"></a>

Davies, R. W., Flint J, Myers S., Mott R. Rapid genotype imputation from sequence without reference panels. *Nat. Genet.* 48, 965-969 (2016)

## Testing <a name="paragraph-testing"></a>

Tests in STITCH are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of STITCH. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Bug reports <a name="paragraph-bugreports"></a>

The best way to get help is to either submit a bug report on GitHub or to consult the forum and mailing list

https://groups.google.com/forum/#!forum/stitch-imputation

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com

## Output format <a name="paragraph-output-format"></a>

STITCH supports writing to both bgzipped vcfs and bgen, see output_format variable


## What method to run <a name="paragraph-what-method"></a>

STITCH can run using one of three "methods" reflecting different underlying statistical and biological models: "diploid", which is the best general method and has the best statistical properties, but has run time proportional to the square of K and so may be slow for large, diverse populations; "pseudoHaploid", which uses statistical approximations that make it less accurate than the diploid method but has run time proportional to K, and so may be suitable for large, diverse populations; and "diploid-inbred", which assumes all samples are completely inbred and as such uses an underlying haplotype based imputation model with run time proportional to K. Note that each of these assumes subjects are diploid, and as such, all methods output diploid genotypes and probabilities.

## Notes on the relationship between run time, RAM and performance <a name="paragraph-time-ram-memory"></a>

STITCH can be run on hundreds of thousands of samples, SNPs, or both. Default parameters are set to give good performance for situations somewhere in the middle. Depending on your application, you may want to tweak default parameters to change how STITCH is run and the relationship between run time, RAM and performance. Here is a brief summary of relevant parameters. See section below for note about K.
* outputSNPBlockSize: STITCH writes out results approximately this many SNPs at a time. Setting this to a larger value will speed up STITCH but use more RAM.
* keepSampleReadsInRAM, inputBundleBlockSize: STITCH converts reads from BAM files into an internal format. These variables control whether all of those are kept in RAM (keepSampleReadsInRAM = TRUE) or not (keepSampleReadsInRAM=FALSE, default) at once. Setting to TRUE decreases runtime but increases RAM usage. inputBundleBlockSize controls whether sampleReads are bundled together to use fewer temporary files on disk. Setting inputBundleBlockSize to a higher integer value will reduce the number of files on the temporary disk and is likely to increase performance, particularly for large sample sizes.
* gridWindowSize: The default gridWindowSize=NA makes imputation run per-SNP, as is standard. Setting this to an integer greater than 0 (e.g. 10000 (base pairs)) bins the genome into physical windows of this size, and runs imputation between those grids. This can considerably speed up imputation of dense regions but will reduce imputation performance.
* S: S controls the number of sets of ancestral haplotypes used and which final results are averaged over. This may be particularly useful for wild or large populations, like humans. S should affect RAM and run time in a roughly linear fashion.

## How to choose K and nGen <a name="paragraph-k-ngen"></a>

A fuller description is given the supplement of the paper given in the citation section, and this is worth a read for anyone planning to use the method in their work.

K is the number of ancestral haplotypes in the model. Larger K allows for more accurate imputation for large samples and coverages, but takes longer and accuracy may suffer with lower coverage. It is usually wise to try a few values of K and assess performance using either external validation, or the distribution of quality scores (e.g. mean / median INFO score). It is likely wise to choose K that both gives you the best performance (accuracy, correlation or quality score distribution) within computational constraints, while also ensuring K is not too large given your sequencing coverage (e.g. try to ensure that each ancestral haplotype gets at least a certain average X of coverage, say 10X, given your number of samples and average depth). 

nGen controls recombination rate between the sequenced samples and the ancestral haplotypes. It is probably fine to set it to 4 * Ne / K given some estimate of effective population size Ne. If you think your population can reasonably approximated as having been founded some number of generations ago / reduced to 2*K that many generations ago, use that generation time estimate. STITCH should be fairly robust to misspecifications of this parameter. 

## About plotting <a name="paragraph-plots"></a>

STITCH generates some plots while running. They are not meant to substitute for a more in depth investigation of imputation performance in your setting, but can be a useful first start to understanding your data and parameter choices. They are described in [plots.md](plots.md), in decreasing order of usefulness. Note that it can be quite useful to set `plot_shuffle_haplotype_attempts=TRUE` to visualize switching in the model.

## About reference panels <a name="paragraph-reference-panels"></a>

STITCH is designed to impute samples without a reference panel. However, STITCH can take as input reference haplotype information. Below the reference panel format is described, and then how certain options affect how the reference haplotype information is used.

The reference panel format used by STITCH is the same as used by QUILT. Please see the QUILT web page [here](https://github.com/rwdavies/QUILT) for specific details of the format. Alternatively, you can download minimal example human data formatted for STITCH from this link `https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_human_reference_example_2018_07_11.tgz` and check out the files `1000GP_Phase3_20.1000000.1100000.legend.gz` and `1000GP_Phase3_20.1000000.1100000.hap.gz`.

Reference panel information is used depending on the following options

`niterations>1`: the reference haplotypes are used to initialize the ancestral haplotypes. After the first iteration, only information from the samples will be used to update the ancestral haplotypes in the EM algorithm.

`niterations==1`: the reference haplotypes are used to initialize the ancestral haplotypes, and the samples are imputed directly from this. Note that in this case, all posfile SNPs must be found in the reference SNPs, as otherwise you would be imputing SNPs with no reference information and hence no information about how to impute them. Note that this condition does not exist for `niterations>1`, as after the first iteration, you fill the ancestral haplotypes with sample information.

Let `nhaps` be the number of reference haplotypes in the haplotype reference file.

`nhaps > K`: A haplotype EM algorithm is run to initialize the ancestral haplotypes using the reference haplotypes.

`nhaps == K`: There are exactly as many reference haplotypes as ancestral haplotypes, and they are used directly.

`nhaps < K`: There are fewer available haplotypes than desired ones, so the available ones are used to fill in the corresponding number of ancestral haplotypes directly, and the remaining ones are filled with noise.

If you have `niterations==1` and `nhaps == K`, each sample will be imputed independently of each other sample without causing batch effects, as they are only imputed from the reference. Note that in general, if you're doing this, you should strongly consider using QUILT if K is large, as this is the specific situation QUILT is designed for. However for small K, STITCH is more accurate, though doesn't directly output phased information.





