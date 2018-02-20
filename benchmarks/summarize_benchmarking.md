
## Benchmarking

Incremental version releases of STITCH have often focused on miscellaneous speed improvements. The purpose of this page is to ensure performance hasn't regressed, as well as to benchmark those speedups across versions and options on a consistent set of data. This is done first using whole chromosome CFW mice data for speed and accuracy and later on a subset of the data to benchmark more options for speed.

The major take aways are
1. No obvious regression in performance across normal version updates
2. Consistent gradual speed improvements of approximately 1.5-3X versus earlier releases
3. Order of magnitude speedups can be obtained for high K using "gridding" approach, with only very slight losses in accuracy


## Whole chromosome CFW mice test

This analysis uses 2,073 mice at approximately 0.015X coverage at X SNPs using K=4. Performance is assessed using X mice on X array


|Version |Options               |Time (min) |Accuracy |
|:-------|:---------------------|:----------|:--------|
|1.4.0   |gridWindowSize=1000   |25.7       |0.975    |
|1.4.0   |gridWindowSize=10000  |22.1       |0.955    |
|1.4.0   |gridWindowSize=100000 |10.9       |0.955    |
|1.4.0   |                      |45.3       |0.973    |
|1.3.7   |                      |44.3       |0.974    |
|1.3.6   |                      |45.6       |0.969    |
|1.3.5   |                      |45.6       |0.976    |
|1.3.4   |                      |45.5       |0.975    |
|1.3.3   |                      |47.8       |0.974    |
|1.2.5   |                      |63.1       |0.965    |
|1.1.1   |                      |67.5       |0.977    |




## Smaller profiling

These analyses test performance for K=4 or K=20, and for either BAMs or CRAMs. These tests are primary for speed, given accuracy on the large whole chromosome has been performed above, and showcase speed for a more reasonable choice for outbred samples (K=20), and show that performance on CRAM samples caught up to BAM samples with version >1.3.0.

### K=4, BAMs

|Version |Options               |Input |K  |Time (min) |
|:-------|:---------------------|:-----|:--|:----------|
|1.4.0   |gridWindowSize=1000   |BAMS  |4  |0          |
|1.4.0   |gridWindowSize=10000  |BAMS  |4  |0          |
|1.4.0   |gridWindowSize=100000 |BAMS  |4  |0          |
|1.4.0   |                      |BAMS  |4  |0          |
|1.3.7   |                      |BAMS  |4  |0          |
|1.3.6   |                      |BAMS  |4  |0          |
|1.3.5   |                      |BAMS  |4  |0          |
|1.3.4   |                      |BAMS  |4  |0          |
|1.3.3   |                      |BAMS  |4  |0          |
|1.2.5   |                      |BAMS  |4  |0.1        |
|1.1.1   |                      |BAMS  |4  |0.1        |

### K=4, CRAMs

|Version |Options               |Input |K  |Time (min) |
|:-------|:---------------------|:-----|:--|:----------|
|1.4.0   |gridWindowSize=1000   |CRAMS |4  |0          |
|1.4.0   |gridWindowSize=10000  |CRAMS |4  |0          |
|1.4.0   |gridWindowSize=100000 |CRAMS |4  |0          |
|1.4.0   |                      |CRAMS |4  |0          |
|1.3.7   |                      |CRAMS |4  |0          |
|1.3.6   |                      |CRAMS |4  |0          |
|1.3.5   |                      |CRAMS |4  |0          |
|1.3.4   |                      |CRAMS |4  |0          |
|1.3.3   |                      |CRAMS |4  |0          |
|1.2.5   |                      |CRAMS |4  |0.1        |
|1.1.1   |                      |CRAMS |4  |NA         |

### K=20, BAMs

|Version |Options               |Input |K  |Time (min) |
|:-------|:---------------------|:-----|:--|:----------|
|1.4.0   |gridWindowSize=1000   |BAMS  |20 |0          |
|1.4.0   |gridWindowSize=10000  |BAMS  |20 |0          |
|1.4.0   |gridWindowSize=100000 |BAMS  |20 |0          |
|1.4.0   |                      |BAMS  |20 |0          |
|1.3.7   |                      |BAMS  |20 |0          |
|1.3.6   |                      |BAMS  |20 |0          |
|1.3.5   |                      |BAMS  |20 |0          |
|1.3.4   |                      |BAMS  |20 |0          |
|1.3.3   |                      |BAMS  |20 |0          |
|1.2.5   |                      |BAMS  |20 |0.1        |
|1.1.1   |                      |BAMS  |20 |0.1        |

### K=20, CRAMs

|Version |Options               |Input |K  |Time (min) |
|:-------|:---------------------|:-----|:--|:----------|
|1.4.0   |gridWindowSize=1000   |CRAMS |20 |0          |
|1.4.0   |gridWindowSize=10000  |CRAMS |20 |0          |
|1.4.0   |gridWindowSize=100000 |CRAMS |20 |0          |
|1.4.0   |                      |CRAMS |20 |0          |
|1.3.7   |                      |CRAMS |20 |0          |
|1.3.6   |                      |CRAMS |20 |0          |
|1.3.5   |                      |CRAMS |20 |0          |
|1.3.4   |                      |CRAMS |20 |0          |
|1.3.3   |                      |CRAMS |20 |0          |
|1.2.5   |                      |CRAMS |20 |0.1        |
|1.1.1   |                      |CRAMS |20 |NA         |
