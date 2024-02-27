# LongcellPre
LongcellPre is an R pipeline to analyze Nanopore long read sequencing dataset based on 10X single cell sequencing toolkit. This pipeline includes preprocessing to do barcode and unique molecular identifier (UMI) assignment to give an accurate isoform quantification. Based on the isoform quantification from LongcellPre, our another  pipeline Longcell incorporates downstream splicing analysis, including identification of highly variable exons and differential alternative splicing analysis between different cell populations.

## Installation
requires:  
- minimap2: https://github.com/lh3/minimap2
- samtools: http://www.htslib.org/
- bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html

You can install `LongcellPre` by `devtools`:
```
devtools::install_github("yuntianf/Longcellsrc")
devtools::install_github("yuntianf/LongcellPre")
```

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Single cell and spatial alternative splicing analysis with long read sequencing
Yuntian Fu, Heonseok Kim, Jenea I. Adams, Susan M. Grimes, Sijia Huang, Billy T. Lau, Anuja Sathe, Paul Hess, Hanlee P. Ji, Nancy R. Zhang
bioRxiv 2023.02.23.529769; doi: https://doi.org/10.1101/2023.02.23.529769
```
