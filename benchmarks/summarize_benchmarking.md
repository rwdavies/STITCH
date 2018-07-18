
## Benchmarking

Incremental version releases of STITCH have often focused on miscellaneous speed improvements. The purpose of this page is to ensure performance hasn't regressed, as well as to benchmark those speedups across versions and options on a consistent set of data. This is done first using whole chromosome CFW mice data for speed and accuracy and later on a subset of the data to benchmark more options for speed.

The major take aways are
1. No obvious regression in performance across normal version updates
2. Consistent gradual speed improvements of approximately 1.5-3X versus earlier releases
3. Order of magnitude speedups can be obtained for high K using "gridding" approach, with only very slight losses in accuracy


## Whole chromosome CFW mice test

This analysis uses 2,073 mice at approximately 0.015X coverage at X SNPs using K=4. Performance is assessed using X mice on X array


|Version |Options                |Time (min) |Accuracy |
|:-------|:----------------------|:----------|:--------|
|1.5.2   |output_format=\'bgen\' |33.2       |0.97     |
|1.5.2   |gridWindowSize=1000    |23.1       |0.97     |
|1.5.2   |gridWindowSize=10000   |20.4       |0.966    |
|1.5.2   |gridWindowSize=100000  |19.1       |0.955    |
|1.5.2   |                       |33.9       |0.966    |
|1.5.1   |                       |34.1       |0.977    |
|1.5.0   |                       |34.1       |0.971    |
|1.4.2   |                       |45.3       |0.976    |
|1.4.1   |                       |45.1       |0.975    |
|1.4.0   |                       |45.1       |0.972    |
|1.3.7   |                       |44         |0.978    |
|1.3.6   |                       |45.5       |0.967    |
|1.3.5   |                       |45.3       |0.968    |
|1.3.4   |                       |45.4       |0.97     |
|1.3.3   |                       |47.6       |0.975    |
|1.2.5   |                       |61.8       |0.971    |
|1.1.1   |                       |66.4       |0.972    |



## Smaller profiling

These analyses test performance for K=4 or K=20, and for either BAMs or CRAMs. These tests are primary for speed, given accuracy on the large whole chromosome has been performed above, and showcase speed for a more reasonable choice for outbred samples (K=20), and show that performance on CRAM samples caught up to BAM samples with version >1.3.0.

### K=4, BAMs

|Version |Options                |Input |K  |Time (min) |
|:-------|:----------------------|:-----|:--|:----------|
|1.5.2   |output_format=\'bgen\' |BAMS  |4  |3.2        |
|1.5.2   |gridWindowSize=1000    |BAMS  |4  |2.5        |
|1.5.2   |gridWindowSize=10000   |BAMS  |4  |2.1        |
|1.5.2   |gridWindowSize=100000  |BAMS  |4  |2          |
|1.5.2   |                       |BAMS  |4  |3.4        |
|1.5.1   |                       |BAMS  |4  |3.3        |
|1.5.0   |                       |BAMS  |4  |3.3        |
|1.4.2   |                       |BAMS  |4  |4.1        |
|1.4.1   |                       |BAMS  |4  |4          |
|1.4.0   |                       |BAMS  |4  |3.9        |
|1.3.7   |                       |BAMS  |4  |4.1        |
|1.3.6   |                       |BAMS  |4  |4          |
|1.3.5   |                       |BAMS  |4  |4.1        |
|1.3.4   |                       |BAMS  |4  |4.1        |
|1.3.3   |                       |BAMS  |4  |4.6        |
|1.2.5   |                       |BAMS  |4  |7.1        |
|1.1.1   |                       |BAMS  |4  |6.5        |

### K=4, CRAMs

|Version |Options                |Input |K  |Time (min) |
|:-------|:----------------------|:-----|:--|:----------|
|1.5.2   |output_format=\'bgen\' |CRAMS |4  |3.5        |
|1.5.2   |gridWindowSize=1000    |CRAMS |4  |2.7        |
|1.5.2   |gridWindowSize=10000   |CRAMS |4  |2.3        |
|1.5.2   |gridWindowSize=100000  |CRAMS |4  |2.2        |
|1.5.2   |                       |CRAMS |4  |3.6        |
|1.5.1   |                       |CRAMS |4  |3.5        |
|1.5.0   |                       |CRAMS |4  |3.5        |
|1.4.2   |                       |CRAMS |4  |4.3        |
|1.4.1   |                       |CRAMS |4  |4.1        |
|1.4.0   |                       |CRAMS |4  |4.1        |
|1.3.7   |                       |CRAMS |4  |4.2        |
|1.3.6   |                       |CRAMS |4  |4.2        |
|1.3.5   |                       |CRAMS |4  |4.2        |
|1.3.4   |                       |CRAMS |4  |4.2        |
|1.3.3   |                       |CRAMS |4  |4.8        |
|1.2.5   |                       |CRAMS |4  |13.3       |
|1.1.1   |                       |CRAMS |4  |NA         |

### K=20, BAMs

|Version |Options                |Input |K  |Time (min) |
|:-------|:----------------------|:-----|:--|:----------|
|1.5.2   |output_format=\'bgen\' |BAMS  |20 |54.2       |
|1.5.2   |gridWindowSize=1000    |BAMS  |20 |16.2       |
|1.5.2   |gridWindowSize=10000   |BAMS  |20 |5.7        |
|1.5.2   |gridWindowSize=100000  |BAMS  |20 |3.8        |
|1.5.2   |                       |BAMS  |20 |54         |
|1.5.1   |                       |BAMS  |20 |54.5       |
|1.5.0   |                       |BAMS  |20 |54         |
|1.4.2   |                       |BAMS  |20 |60.9       |
|1.4.1   |                       |BAMS  |20 |66.9       |
|1.4.0   |                       |BAMS  |20 |67         |
|1.3.7   |                       |BAMS  |20 |63.1       |
|1.3.6   |                       |BAMS  |20 |59.6       |
|1.3.5   |                       |BAMS  |20 |63.1       |
|1.3.4   |                       |BAMS  |20 |62.9       |
|1.3.3   |                       |BAMS  |20 |57.8       |
|1.2.5   |                       |BAMS  |20 |171.7      |
|1.1.1   |                       |BAMS  |20 |184.2      |

### K=20, CRAMs

|Version |Options                |Input |K  |Time (min) |
|:-------|:----------------------|:-----|:--|:----------|
|1.5.2   |output_format=\'bgen\' |CRAMS |20 |54.5       |
|1.5.2   |gridWindowSize=1000    |CRAMS |20 |15.2       |
|1.5.2   |gridWindowSize=10000   |CRAMS |20 |5.8        |
|1.5.2   |gridWindowSize=100000  |CRAMS |20 |3.9        |
|1.5.2   |                       |CRAMS |20 |54.4       |
|1.5.1   |                       |CRAMS |20 |54.5       |
|1.5.0   |                       |CRAMS |20 |54.7       |
|1.4.2   |                       |CRAMS |20 |60.7       |
|1.4.1   |                       |CRAMS |20 |67.3       |
|1.4.0   |                       |CRAMS |20 |66.9       |
|1.3.7   |                       |CRAMS |20 |62.8       |
|1.3.6   |                       |CRAMS |20 |58.5       |
|1.3.5   |                       |CRAMS |20 |62.9       |
|1.3.4   |                       |CRAMS |20 |62.8       |
|1.3.3   |                       |CRAMS |20 |63.5       |
|1.2.5   |                       |CRAMS |20 |178.4      |
|1.1.1   |                       |CRAMS |20 |NA         |
