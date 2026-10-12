# LongcellPre
LongcellPre is an R pipeline to analyze Nanopore long read sequencing dataset based on 10X single cell sequencing toolkit. This pipeline includes preprocessing to do barcode and unique molecular identifier (UMI) assignment to give an accurate isoform quantification. Based on the isoform quantification from LongcellPre, our another  pipeline Longcell incorporates downstream splicing analysis, including identification of highly variable exons and differential alternative splicing analysis between different cell populations.

## Installation
requires:  
- minimap2: https://github.com/lh3/minimap2
- samtools: http://www.htslib.org/
- bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html

Use a compatible R/Bioconductor environment, a C++17 compiler and zlib.
The command-line entry point also uses Python through R's `argparse` package.
The paper used samtools 1.18; the mapping command requires the modern
`samtools sort -o output.bam` syntax. Removing `-o` is not a fix for modern
samtools ([manual](https://www.htslib.org/doc/samtools-sort.html)).

LongcellPre 1.3.0 includes its C++ backend. A separate installation of
`Longcellsrc` is no longer required. Install this local source with the bundle's
`INSTALL.R`, or from R after installing the declared dependencies:

```r
install.packages("/path/to/LongcellPre_1.3.0.tar.gz", repos = NULL, type = "source")
```

For the published integrated version, one repository installation is sufficient:

```r
devtools::install_github("yuntianf/LongcellPre", dependencies = NA)
```

`dependencies = NA` installs required dependencies without downloading optional
whole-genome BSgenome packages. The pipeline can use your indexed reference FASTA.
If you are using an unpublished local revision, use the source installation above;
GitHub installation retrieves whichever version is currently published there.
Existing `library(LongcellPre)`,
`RunLongcellPre()` calls and command-line arguments are unchanged. An older
Longcellsrc installation can remain for other projects; this package does not
load or use it. Native helper routines are internal to LongcellPre. Code that
directly calls `Longcellsrc::` helpers is outside the LongcellPre public interface.
Use a fresh output directory when switching versions, because the pipeline
checks package versions before reusing cached results.

Use the execution script from the same installation:

```r
system.file("exec", "RunLongcellPre.R", package = "LongcellPre")
```

## Automated checks and source builds

The [Package CI workflow](.github/workflows/ci.yml) runs on pushes, pull requests,
and manual requests. It uses the validated Linux R 4.4 / Bioconductor 3.20
environment, builds the C++ code, checks and installs the source package, and runs
the regression tests, including real mapping and separate-worker checks.
It also compares native barcode, UMI and isoform outputs with saved synthetic
reference results. Longcellsrc must be absent, so the checks cannot accidentally
use the old separate backend.

Failed builds, failed or skipped regression tests, changed native results, and new
or changed package-check warnings fail CI. The exact pre-existing documentation
warning is recorded in `.github/ci/known-warnings.json`; existing notes remain in
the logs. Successful runs provide a source `.tar.gz` and SHA-256 checksum as an
Actions artifact. Logs are retained for 30 days. This delivers a checked source
build; public release publication remains a maintainer action.

CI checks the documented Linux environment. It does not run the large external
demo on every push or certify other platforms and all sequencing datasets. See
[CONTRIBUTING.md](CONTRIBUTING.md) for local checks, release steps, and how to require
CI to pass before merging. After adding the workflow to GitHub, inspect its first
Actions run to confirm the hosted setup as well as the locally tested commands.

## Input preparation

- Supply raw reads containing the library adapter, cell barcode, and UMI.
  Previously adapter-trimmed reads may leave too few tags to match.
- Prefer the barcodes detected in your sample, for example Cell Ranger's filtered
  `barcodes.tsv.gz`. A full chemistry whitelist is accepted but can be much slower;
  it was the key difference identified in [issue #7](https://github.com/yuntianf/LongcellPre/issues/7).
- Match the toolkit to the library: `-t 5` for 5-prime and `-t 3` for 3-prime.
  The provided demo uses `-t 5` and contains 244 cells.
- Keep the mapping FASTA and GTF consistent. Sequence export uses this same FASTA
  and its chromosome names. Lower-level calls without a FASTA can use an already
  installed BSgenome package; analysis never installs reference packages automatically.
- For a demo file whose `.gz` suffix does not match its contents
  ([issue #13](https://github.com/yuntianf/LongcellPre/issues/13)), check with
  `file example.fastq.gz` and `gzip -t example.fastq.gz`. A plain FASTQ can be
  renamed to `.fastq` or genuinely compressed; do not assume the suffix proves
  compression. `python exec/prepare_fastq.py example.fq.gz example.fastq.gz`
  makes a correctly compressed copy without changing the input. The externally
  hosted demo archive has not been replaced.
- PARSE combinatorial barcodes are not implemented. `protocol="other"` configures
  a contiguous barcode/UMI layout; it is not a general multi-round barcode parser.

## Workflow
The simplist way to run `LongcellPre` is to apply:
```
Rscript ./LongcellPre/exec/RunLongcellPre.R -f {FASTQ} -b {BARCODES} -t {TOOLKIT} -q {PROTOCOL} -g {GENOME} -n {GENOME_NAME} --gtf {GTF} -o {OUTDIR}
```
The execution file can also be downloaded from https://github.com/yuntianf/LongcellPre/blob/main/exec/RunLongcellPre.R


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

We provide a demo data with 244 cells and 3 genes. This data is a subset of the colorectal metastasis sample we used in the paper. The data and corresponding annotations can be downloaded from: 
https://www.dropbox.com/scl/fo/21tw8rrkaancani0fzq3t/AKNHUk06onR2c2dYuB4wXWY?rlkey=1zikug28qr9ziw2cdsgelrm9p&st=ypm9m00i&dl=0

The demo data can be processed by:
```
Rscript RunLongcellPre.R -f example.fastq.gz -b barcodes.txt -t 5 -q 10X -g genome.fa -n hg38 --gtf gencode.v39.sub.gtf -o ./demo/ -c 4 -m multicore
```

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
│   └── BarcodeMatchIso.txt: The cell barcode, UMI, and read alignment for each read.
├── out
│   ├── gene: The sparse matrix for the cell by gene matrix.
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── iso_count.txt: The isoform quantification in each cell (without isoform annotation).
│   ├── isoform: The sparse matrix for the cell by isoform matrix.
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── reads_annot.csv: The read annotation for the UMI collapsed fastq file.
│   └── UMI_collapsed.fq.gz: Reference-reconstructed reads after UMI collapsing; see limitations below.
│   └── UMI_collapsed.bam: The mapping result for UMI_collapsed.fq.gz.
│   └── UMI_collapsed.bam.bai
├── arg.log: The log file to record used parameters.
└── polish.fq.gz
```
The isoform quantification output `iso_count.txt` can be directly analyzed by our downstream package `Longcell` (https://github.com/yuntianf/Longcell).

The matrix files are actually written as `barcodes.tsv.gz`, `features.tsv.gz`,
and `matrix.mtx.gz`, with features in rows and cells in columns.
`BarcodeMatch.txt` is removed at successful completion because its information
has been joined into `BarcodeMatchIso.txt`. Use `--keep_intermediates TRUE` to retain it.

`iso_count.txt` contains `cell`, `gene`, `isoform`, `count`, and `polyA`.
`isoform` is a string of genomic exon intervals, not necessarily an annotated
transcript ID. `count` is the final count after UMI clustering and scatter
filtering. Older diagnostic outputs also included `size` (raw reads) and
`cluster` (UMI count before scatter filtering). `polyA` summarizes tail evidence
among contributing reads; it can be fractional and does not mean an exported
sequence retains a poly(A) tail ([issue #2](https://github.com/yuntianf/LongcellPre/issues/2)).

### Meaning of the reconstructed reads

`UMI_collapsed.fq.gz` is generated by extracting reference-genome sequence at the
corrected exon coordinates. Its qualities are synthetic `~` characters, and
each row contributes `round(count)` reads using R's existing rounding rule.
These sequences preserve the inferred exon structure, but do not preserve
sample variants, original base qualities, or original poly(A) tails. Reference
`N` bases can remain. Use original reads for variant analysis
([issue #9](https://github.com/yuntianf/LongcellPre/issues/9)).

The pipeline writes FASTQ and matching read-to-cell CSV in batches within each
gene (default: 10,000 reads or 10 million bases, whichever is reached first).
One read longer than the base limit is emitted alone. Tune `--batch_reads` and
`--batch_bases` to change these bounds. Export does not use parallel workers.
`out/export_qc.tsv` records each gene's exported reads, bases and reference Ns.
The mapping FASTA now supplies reconstruction sequence as well, avoiding a
second reference with potentially different bases. Missing contigs cause an
error instead of silently dropping genes.

The bundled C++ backend replaces the fixed UMI score arrays with a bounded cache and
streams compatible isoform pairs. The final graph and several R intermediates
can still grow quadratically. `--max_graph_edges` defaults to 5,000,000 edges
per cell/gene; exceeding it reports the gene and cell and stops the run. It is
not a total RAM limit. `out/umi_qc.tsv` records successful, filtered and failed
genes. `--error_policy continue` explicitly permits partial results with a
warning; default `stop` prevents unnoticed loss of genes. Large whole-transcriptome
datasets still require representative pilots; reducing workers can reduce
concurrent memory use. The 180-million-read case has not been validated.

Input MD5 checksums, package versions and settings are saved in
`input_manifest.tsv`, `run_manifest.rds` and `provenance.txt`. Hashing reads the
input files once and is included in measured runtime. Changed inputs or settings
require a fresh output directory or `--force_rerun TRUE`. An interrupted run is
marked incomplete and also requires this explicit rerun to avoid mixing outputs.
Alignment commands
require samtools >= 1.3, preserve modern `sort -o` syntax, and validate staged
BAMs before replacing previous BAM/index pairs. Use Linux or WSL for this pipeline.

If `iso_count.txt` already exists, regenerate only the export and its BAM:

```r
LongcellPre::UMI_consensus_out(
  umi_count = read.delim("output/out/iso_count.txt"),
  dir = "output/out",
  genome_path = "genome.fa", genome_name = "hg38",
  gene_bed = readRDS("output/annotation/gene_bed.rds"),
  cores = 4, force_fastq_out = TRUE, force_map = TRUE
)
```

`isoformCount2Reads()` still returns reads in memory by default for callers that
need them. Pass `return_reads=FALSE` and optionally `annotation_file="reads.csv"`
to use its lower-memory mode directly. `--to_isoform FALSE` skips annotated
isoform matrix generation; it does not skip FASTQ/BAM export.

For the isoform quantification formatted as cell by isoform matrix, we recommend the combination of LongcellPre and IsoQuant (https://github.com/ablab/IsoQuant), which has the overall best performance in our benchmark. The input for IsoQuant is `UMI_collapsed.bam` and `reads_annot.csv`:
```
python isoquant.py --reference $GENOME_PATH \
--genedb $GTF \
--bam UMI_collapsed.bam \
--read_group file:reads_annot.csv:0:1:, \
--data_type nanopore -o $OUTDIR \
--report_novel_unspliced true \
--clean_start
```

The paper benchmarked IsoQuant **v3.3**. The command above is the legacy example;
it is not a version-independent transcript-discovery recipe. Current IsoQuant
has `--polya_trimmed {none,stranded,all}` and `--polya_requirement` options
([official documentation](https://ablab.github.io/IsoQuant/cmd.html)).
For reference-reconstructed LongcellPre reads, poly(A) sequence is absent even
when biological tail evidence was recorded. `--polya_trimmed all` marks every
read as polyadenylated; `--polya_requirement never` relaxes evidence requirements.
Either can change novel-transcript calls and should be evaluated explicitly,
especially for one- and two-exon transcripts. Check the installed version's
`--full_help`, record the options, and validate calls against the original reads.
These options have not been benchmarked here and should not be described as the
settings used for Figure 2 without recovering the original benchmark scripts
([issue #17](https://github.com/yuntianf/LongcellPre/issues/17)).

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Fu Y, Kim H, Roy S, et al. Single cell and spatial alternative splicing analysis with Nanopore long read sequencing[J]. Nature Communications, 2025, 16(1): 6654.
```
