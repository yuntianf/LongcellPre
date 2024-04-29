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
minimap_bed_path = "path of your bed annotation for minimap2, can be generated from gtf" //unnecessary
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
For a more detailed explanation of steps for LongcellPre, please refer to https://github.com/yuntianf/LongcellPre/blob/main/vignettes/LongcellPre.ipynb

The isoform quantification would be output in a long table with 5 columns, including "cell", "gene", "isoform", "count" and "polyA". Due to the UMI scattering filtering, some values in "count" may not be integer. And "polyA" is the average of polyA exitence with a UMI cluster.
The the gtf isoform annotation is provided, each isoform would be aligned to the canonical isoform and a cell by isoform matrix would also be generated.

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Single cell and spatial alternative splicing analysis with long read sequencing
Yuntian Fu, Heonseok Kim, Jenea I. Adams, Susan M. Grimes, Sijia Huang, Billy T. Lau, Anuja Sathe, Paul Hess, Hanlee P. Ji, Nancy R. Zhang
bioRxiv 2023.02.23.529769; doi: https://doi.org/10.1101/2023.02.23.529769
```
