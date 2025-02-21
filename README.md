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
The singularity image for this package is under construction.

## Workflow
The simplist way to run `LongcellPre` is to apply:
```
Rscript ./LongcellPre/exec/RunLongcellPre.R -f {FASTQ} -b {BARCODES} -t {TOOLKIT} -q {PROTOCOL} -g {GENOME} -n {GENOME_NAME} --gtf {GTF} -o {OUTDIR}
```
The basic parameters for this pipeline are shown here:
```
Pipeline for single cell Nanopore RNA-seq preprocessing

options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        The path for the input fastq file
  -b BARCODE, --barcode BARCODE
                        The path for the cell barcode whitelist
  -t {5,3}, --toolkit {5,3}
                        The toolkit used in sequencing, should be 5 or 3
  -q {10X,VISIUM,Curio,other}, --protocol {10X,VISIUM,Curio,other}
                        The sequencing protocol, ex. '10X', 'VISIUM', 'Curio'
  -g GENOME_PATH, --genome_path GENOME_PATH
                        The path of the genome reference
  -n GENOME_NAME, --genome_name GENOME_NAME
                        the genome name used for mapping, ex. hg38
  --gtf GTF             The path of the gtf annotation
  --gene_bed_path GENE_BED_PATH
                        The path of the gene bed annotation
  --minimap_bed_path MINIMAP_BED_PATH
                        The path of your bed annotation for minimap2, can be
                        generated from GTF/GFF3 with ‘paftools.js gff2bed
                        anno.gtf’
  -o WORK_DIR, --work_dir WORK_DIR
                        The output directory
  --to_isoform TO_ISOFORM
                        A flag to indicate if the cell by isoform matrix
                        should be generated
  -c CORES, --cores CORES
                        The number of cores used for parallization
  -m {sequential,multisession,multicore,cluster}, --mode {sequential,multisession,multicore,cluster}
                        The mode for parallization. The parallization is
                        implemented with future.apply, the feasible modes can
                        be 'sequential','multicore','cluster'
  --minimap2 MINIMAP2   The path of the minimap2
  --samtools SAMTOOLS   The path of the samtools
  --bedtools BEDTOOLS   The path of the bedtools
```
To view all parameters for this pipeline, you can run `Rscript ./LongcellPre/exec/RunLongcellPre.R -h --full`

We provide a demo data with 200 cells and 3 genes. This data is a subset of the colorectal metastasis sample we used in the paper. The data and corresponding annotations can be downloaded from: 
https://www.dropbox.com/scl/fo/21tw8rrkaancani0fzq3t/AKNHUk06onR2c2dYuB4wXWY?rlkey=1zikug28qr9ziw2cdsgelrm9p&st=ypm9m00i&dl=0

For a more detailed explanation of steps for LongcellPre, please refer to https://github.com/yuntianf/LongcellPre/blob/main/vignettes/LongcellPre.ipynb

## output
The output of the LongcellPre pipeline includes:
```
├── annotation
│   ├── exon_gtf.rds: The exon annotation for each gene given the gtf annotation, stored in RDS format.
│   ├── exon_gtf.txt: The exon annotation for each gene given the gtf annotation, stored in tsv format.
│   ├── gene_bed.rds: The annotation non-overlapped exon bins for each gene given the gtf annotation, stored in RDS format.
│   ├── gene_bed.txt: The annotation non-overlapped exon bins for each gene given the gtf annotation, stored in tsv format.
├── bam
│   ├── polish.bam: The mapped bam from reads with adapters trimmed.
│   └── polish.bam.bai
├── BarcodeMatch
│   ├── adapterNeedle.txt: The summary statistics for the Needleman score of the kmer adapter sequence in each reads compared to its know sequence.
│   ├── BarcodeMatchIso.txt: 
│   └── BarcodeMatch.txt
├── out
│   ├── gene: The sparse matrix for the cell by gene matrix.
│   │   ├── barcodes.tsv
│   │   ├── features.tsv
│   │   └── matrix.mtx
│   ├── iso_count.txt: The isoform quantification in each cell (without isoform annotation).
│   ├── isoform: The sparse matrix for the cell by isoform matrix.
│   │   ├── barcodes.tsv
│   │   ├── features.tsv
│   │   └── matrix.mtx
│   ├── reads_annot.csv: The read annotation for the UMI collapsed fastq file.
│   └── UMI_collapsed.fq.gz: The polished reads after UMI collapsion. This can be input for other tools.
└── polish.fq.gz
```

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Single cell and spatial alternative splicing analysis with long read sequencing
Yuntian Fu, Heonseok Kim, Jenea I. Adams, Susan M. Grimes, Sijia Huang, Billy T. Lau, Anuja Sathe, Paul Hess, Hanlee P. Ji, Nancy R. Zhang
bioRxiv 2023.02.23.529769; doi: https://doi.org/10.1101/2023.02.23.529769
```
