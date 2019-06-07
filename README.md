STITCH - Sequencing To Imputation Through Constructing Haplotypes
=================================================================
**__Current Version: 1.5.7__**
Release date: April 17, 2019

[![Build Status](https://img.shields.io/travis/rwdavies/STITCH/master.svg)](https://img.shields.io/travis/rwdavies/STITCH/master.svg)

Changes in latest version

1. Push CC through to htslib to robustify against compilation issues 

For details of past changes please see [CHANGELOG](CHANGELOG.md).

STITCH is an R program for reference panel free, read aware, low coverage sequencing genotype imputation. STITCH runs on a set of samples with sequencing reads in BAM format, as well as a list of positions to genotype, and outputs imputed genotypes in VCF format. 

For the old website, please see https://www.well.ox.ac.uk/~rwdavies/stitch.html

## Installation and quick start on real data example

### Quick start on Linux and Mac

Conda instructions below. Alternatively, STITCH can be installed in a few ways. The simplest way to get a release is as follows. First, install R. Then, do the following 
```
git clone --recursive https://github.com/rwdavies/STITCH.git
cd STITCH
./scripts/install-dependencies.sh
cd releases
wget https://github.com/rwdavies/stitch/releases/download/1.5.7/STITCH_1.5.7.tar.gz ## or curl -O
R CMD INSTALL STITCH_1.5.7.tar.gz
```

A quick test on real data can be performed using 
```
# test on CFW mouse data
wget https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
# or curl -O https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
tar -xzvf STITCH_example_2016_05_10.tgz
./STITCH.R --chr=chr19 --bamlist=bamlist.txt --posfile=pos.txt --genfile=gen.txt --outputdir=./ --K=4 --nGen=100 --nCores=1
# if this works the file stitch.chr19.vcf.gz will be created
```

To install the latest development code in the repository, use `./scripts/build-and-install.sh`. To install alternative releases, either download other releases from Github, or use the historical `releases` directory. 

If you experience a problem with installation, you can either try conda below. Alternatively, if you see an error similar to ```error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory```, then either ask your system administrator to install gmp, mpfr and mpc for you, or try running the following before R CMD INSTALL
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

### Install using conda

STITCH (as r-stitch) can be installed using [conda](https://conda.io/miniconda.html). Full tutorials can be found elsewhere, but briefly, something like this should work
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install r-stitch -c defaults -c bioconda -c conda-forge
source activate
R -e 'library("STITCH")'
```
Note that currently the command like `STITCH.R` is not included with the bioconda installation, so from the command line, you can either run something like `R -e 'library("STITCH"); STITCH(chr="chr19", bamlist="bamlist.txt", posfile="pos.txt", genfile="gen.txt", outputdir="./", K=4, nGen=100, nCores=1)'`, or clone the repo to get `STITCH.R`. 

### Interactive start
1. Install R if not already installed.
2. Install R dependencies parallel, Rcpp and RcppArmadillo from CRAN (using the "install.packages" option within R)
3. Install [bgzip](http://www.htslib.org/) and make it available to your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)). This can be done using a system installation, or doing a local installation and either modifying the PATH variable using code like ```export PATH=/path/to/dir-with-bgzip-binary/:$PATH```, or through R, doing something like ```Sys.setenv( PATH = paste0("/path/to/dir-with-bgzip-binary/:", Sys.getenv("PATH")))```. You'll know samtools is available if you run something like ```system("which bgzip")``` in R and get the path to bgzip
4. Install STITCH. First, download the latest STITCH tar.gz from the releases folder above. Second, install by opening R and using install.packages, giving install.packages the path to the downloaded STITCH tar.gz. This should install SeqLib automatically as well.
5. Download example dataset [STITCH_example_2016_05_10.tgz](https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz).
6. Run STITCH. Open R, change your working directory using setwd() to the directory where the example tar.gz was unzipped, and then run ```STITCH(tempdir = tempdir(), chr = "chr19", bamlist = "bamlist.txt", posfile = "pos.txt", genfile = "gen.txt", outputdir = paste0(getwd(), "/"), K = 4, nGen = 100, nCores = 1)```. Once complete, a VCF should appear in the current working directory named stitch.chr19.vcf.gz

## Help, command line interface and common options

For a full list of options, in R, query ```?STITCH```, or from the command line, ```STITCH --help```. For a brief writeup of commonly used variables, see [Options.md](Options.md). To pass vectors using the command line, do something like ```STITCH.R --refillIterations='c(3,40)'``` or ```STITCH.R --reference_populations='c("CEU","GBR")'```.

## Benchmarks 

One can see some speed benchmarks in [benchmarks/summarize_benchmarking.md](benchmarks/summarize_benchmarking.md)

## Examples

In the examples directory, there is a script which contains examples using real mouse and human data. One can either run this interactively in R, or run all examples using ```./examples/example.R```.

## License

STITCH and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE).

## Testing

Tests in STITCH are split into unit or acceptance run using ```./scripts/test-unit.sh``` and ```./scripts/test-acceptance.sh```. To run all tests use ```./scripts/all-tests.sh```, which also builds and installs a release version of STITCH. To make compilation go faster do something like ```export MAKE="make -j 8"```.

## Citation

Davies, R. W., Flint J, Myers S., Mott R. Rapid genotype imputation from sequence without reference panels. *Nat. Genet.* 48, 965-969 (2016)

## Contact and bug reports

The best way to get help is to either submit a bug report on GitHub or to consult the forum and mailing list

https://groups.google.com/forum/#!forum/stitch-imputation

For more detailed questions or other concerns please contact Robert Davies robertwilliamdavies@gmail.com

## Output format

STITCH supports writing to both bgzipped vcfs and bgen, see output_format variable

## What method to run

STITCH can run using one of three "methods" reflecting different underlying statistical and biological models: "diploid", which is the best general method and has the best statistical properties, but has run time proportional to the square of K and so may be slow for large, diverse populations; "pseudoHaploid", which uses statistical approximations that make it less accurate than the diploid method but has run time proportional to K, and so may be suitable for large, diverse populations; and "diploid-inbred", which assumes all samples are completely inbred and as such uses an underlying haplotype based imputation model with run time proportional to K. Note that each of these assumes subjects are diploid, and as such, all methods output diploid genotypes and probabilities.

## Notes on the relationship between run time, RAM and performance

STITCH can be run on hundreds of thousands of samples, SNPs, or both. Default parameters are set to give good performance for situations somewhere in the middle. Depending on your application, you may want to tweak default parameters to change how STITCH is run and the relationship between run time, RAM and performance. Here is a brief summary of relevant parameters. See section below for note about K.
* outputSNPBlockSize: STITCH writes out results approximately this many SNPs at a time. Setting this to a larger value will speed up STITCH but use more RAM.
* keepSampleReadsInRAM, inputBundleBlockSize: STITCH converts reads from BAM files into an internal format. These variables control whether all of those are kept in RAM (keepSampleReadsInRAM = TRUE) or not (keepSampleReadsInRAM=FALSE, default) at once. Setting to TRUE decreases runtime but increases RAM usage. inputBundleBlockSize controls whether sampleReads are bundled together to use fewer temporary files on disk. Setting inputBundleBlockSize to a higher integer value will reduce the number of files on the temporary disk and is likely to increase performance, particularly for large sample sizes.
* gridWindowSize: The default gridWindowSize=NA makes imputation run per-SNP, as is standard. Setting this to an integer greater than 0 (e.g. 10000 (base pairs)) bins the genome into physical windows of this size, and runs imputation between those grids. This can considerably speed up imputation of dense regions but will reduce imputation performance.

## Note on the selection of K and nGen

A fuller description is given the supplement of the paper given in the [citation](#citation), and this is worth a read for anyone planning to use the method in their work.

K is the number of ancestral haplotypes in the model. Larger K allows for more accurate imputation for large samples and coverages, but takes longer and accuracy may suffer with lower coverage. It is usually wise to try a few values of K and assess performance using either external validation, or the distribution of quality scores (e.g. mean / median INFO score). It is likely wise to choose K that both gives you the best performance (accuracy, correlation or quality score distribution) within computational constraints, while also ensuring K is not too large given your sequencing coverage (e.g. try to ensure that each ancestral haplotype gets at least a certain average X of coverage, say 10X, given your number of samples and average depth). 

nGen controls recombination rate between the sequenced samples and the ancestral haplotypes. It is probably fine to set it to 4 * Ne / K given some estimate of effective population size Ne. If you think your population can reasonably approximated as having been founded some number of generations ago / reduced to 2*K that many generations ago, use that generation time estimate. STITCH should be fairly robust to misspecifications of this parameter. 
