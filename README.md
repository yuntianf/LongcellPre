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
BiocManager::install(c("Rsamtools", "GenomicRanges", "IRanges", "BSgenome", "GenomicFeatures"))
devtools::install_github("yuntianf/LongcellPre",dependencies=TRUE)
```

## Workflow
The simplist way to run `LongcellPre` is to apply:
```
library(future)
library(Longcellsrc)
library(LongcellPre)

# update your input path here
fastq = "path of your input fastq or fq.gz"
barcodes = "path of your input cell barcode whitelist"
gtf_path = "path of your gtf annotation"
genome_path = "path of your genome annotation"
minimap_bed_path = "path of your bed annotation for minimap2, can be generated from GTF/GFF3 with ‘paftools.js gff2bed anno.gtf’" //unnecessary
genome_name = "the genome name used for mapping, ex. hg38"
toolkit = your 10X sequencing toolkit
work_dir = "The output directory"

# specify the path for those tools
samtools = "samtools"
minimap2 = "minimap2"
bedtools = "bedtools"

RunLongcellPre(fastq = fastq,barcode_path = barcodes,toolkit = toolkit,
               genome_path = genome_path,genome_name = genome_name,
               gtf_path = gtf_path,minimap_bed_path = minimap_bed_path,work_dir = work_dir,
               samtools = samtools, minimap2 = minimap2,bedtools = bedtools,cores = 4, strategy="multicore")
```

We provide a demo data with 200 cells and 3 genes. This data is a subset of the colorectal metastasis sample we used in the paper. The data and corresponding annotations can be downloaded from: 
https://www.dropbox.com/scl/fo/21tw8rrkaancani0fzq3t/AKNHUk06onR2c2dYuB4wXWY?rlkey=1zikug28qr9ziw2cdsgelrm9p&st=ypm9m00i&dl=0

For a more detailed explanation of steps for LongcellPre, please refer to https://github.com/yuntianf/LongcellPre/blob/main/vignettes/LongcellPre.ipynb

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Single cell and spatial alternative splicing analysis with long read sequencing
Yuntian Fu, Heonseok Kim, Jenea I. Adams, Susan M. Grimes, Sijia Huang, Billy T. Lau, Anuja Sathe, Paul Hess, Hanlee P. Ji, Nancy R. Zhang
bioRxiv 2023.02.23.529769; doi: https://doi.org/10.1101/2023.02.23.529769
```
