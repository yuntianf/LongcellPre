# LongcellPre

LongcellPre processes Nanopore single-cell RNA sequencing reads into gene and
isoform counts. It assigns cell barcodes, corrects splice sites, and collapses
unique molecular identifiers (UMIs). The resulting isoform counts can be used
with [Longcell](https://github.com/yuntianf/Longcell) for downstream alternative
splicing analysis.

**Version 1.3.0 includes the native backend. You only need to install LongcellPre;
Longcellsrc is no longer required.** See the [version history](NEWS.md) for changes.

## Installation

Run the pipeline on Linux, or inside WSL on Windows. You will need:

- R and a compatible [Bioconductor](https://bioconductor.org/install/) release.
  LongcellPre 1.3.0 has been tested with R 4.4.3 and Bioconductor 3.20.
- A C++17 compiler and zlib development headers to install the package from source.
- [minimap2](https://github.com/lh3/minimap2),
  [samtools](https://www.htslib.org/) version 1.3 or later, and
  [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html),
  available on your `PATH`.
- Python 3 for the command-line interface.

In R, install LongcellPre and its required R dependencies:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("yuntianf/LongcellPre", dependencies = NA, update = FALSE)
library(LongcellPre)
```

[BiocManager](https://bioconductor.github.io/BiocManager/reference/install.html)
selects dependency repositories for your R/Bioconductor version. The installation
above skips optional genome data packages: you supply the reference genome as a
FASTA file when running the pipeline.

## Prepare your inputs

- **Reads:** a FASTQ file, plain or gzip-compressed, containing the library adapter,
  cell barcode, and UMI. Use reads before adapter/barcode trimming.
- **Cell barcodes:** one barcode per line, preferably the cells detected in your
  sample, such as Cell Ranger's filtered `barcodes.tsv.gz`. Cell Ranger suffixes
  such as `-1` are accepted. A full chemistry whitelist is also accepted but can
  substantially increase runtime.
- **Reference genome:** a FASTA file for the assembly used in your experiment.
  The pipeline uses this reference for both mapping and reconstructed sequences.
- **Annotation:** a GTF file from the same assembly, with chromosome names that
  match the FASTA. A GTF is needed to produce the annotated gene and isoform matrices.

For 10X libraries, choose `-t 5` for the 5-prime preset or `-t 3` for the 3-prime
preset. Both use 16-base cell barcodes; the UMI lengths are 10 and 12 bases,
respectively. Check that the preset matches your library chemistry.

Presets also exist for `VISIUM` and `Curio`. Custom contiguous barcode/UMI layouts
require the R interface with `protocol = "other"` and explicit adapter/layout
settings. PARSE combinatorial barcodes are not currently supported.

## Run LongcellPre

### Command line

In a Linux or WSL terminal, locate the script from your installed package and run
it with your input paths:

```sh
LONGCELL_SCRIPT=$(Rscript -e 'cat(system.file("exec", "RunLongcellPre.R", package = "LongcellPre", mustWork = TRUE))')

Rscript "$LONGCELL_SCRIPT" \
  -f reads.fastq.gz -b barcodes.tsv.gz \
  -t 5 -q 10X \
  -g genome.fa -n hg38 --gtf annotation.gtf \
  -o output -c 4 -m multicore
```

The pipeline creates the output directory. The main options are:

| Option | Meaning |
| --- | --- |
| `-f` | Input FASTQ file. |
| `-b` | Cell barcode file. |
| `-t` | Library preset: `5` or `3`. |
| `-q` | Protocol; defaults to `10X`. |
| `-g` | Reference genome FASTA. |
| `-n` | Reference assembly name, such as `hg38` or `mm39`; this does not download a genome. |
| `--gtf` | Gene and transcript annotation. |
| `-o` | Output directory. Use a separate directory for each sample. |
| `-c`, `-m` | Number of workers and parallel mode. Use `-c 1 -m sequential` for a single worker. |

For all options, run:

```sh
Rscript "$LONGCELL_SCRIPT" --help --full
```

If the external tools are outside your `PATH`, provide their executable paths
with `--minimap2`, `--samtools`, and `--bedtools`.

### R

Alternatively, run the same analysis from R:

```r
library(LongcellPre)

RunLongcellPre(
  fastq_path = "reads.fastq.gz",
  barcode_path = "barcodes.tsv.gz",
  toolkit = 5,
  protocol = "10X",
  genome_path = "genome.fa",
  genome_name = "hg38",
  gtf_path = "annotation.gtf",
  work_dir = "output",
  cores = 4,
  mode = "multicore"
)
```

Use `?RunLongcellPre` for R function documentation.

### Demo

The [demo dataset](https://www.dropbox.com/scl/fo/21tw8rrkaancani0fzq3t/AKNHUk06onR2c2dYuB4wXWY?rlkey=1zikug28qr9ziw2cdsgelrm9p&dl=0)
contains reads for three genes from the colorectal metastasis sample used in the
paper. Download and extract the files, then run the following from the directory
containing `example.fq.gz`, `barcodes.txt`, `genome.fa`, and `gencode.v39.sub.gtf`:

```sh
LONGCELL_SCRIPT=$(Rscript -e 'cat(system.file("exec", "RunLongcellPre.R", package = "LongcellPre", mustWork = TRUE))')

Rscript "$LONGCELL_SCRIPT" \
  -f example.fq.gz -b barcodes.txt \
  -t 5 -q 10X \
  -g genome.fa -n hg38 --gtf gencode.v39.sub.gtf \
  -o demo-output -c 4 -m multicore
```

Some copies of `example.fq.gz` contain plain FASTQ despite the `.gz` extension.
LongcellPre detects the file contents and reads it correctly, with a warning;
no conversion is needed for this command.

In the validated demo run, the gene matrix contains 3 genes and 244 cells, and
the annotated isoform matrix contains 21 isoforms and 244 cells. Each matrix has
a total count of 3,918.

## Outputs

Paths below are relative to the output directory:

| Output | Contents |
| --- | --- |
| `out/iso_count.txt` | Isoform counts before assignment to annotated transcripts; input for downstream Longcell analysis. |
| `out/gene/` | Annotated gene count matrix. |
| `out/isoform/` | Annotated transcript/isoform count matrix. |
| `out/UMI_collapsed.fq.gz` | Reads reconstructed from the reference at the corrected exon coordinates. |
| `out/UMI_collapsed.bam`, `out/UMI_collapsed.bam.bai` | Alignments of the reconstructed reads and their index. |
| `out/reads_annot.csv` | Reconstructed read names (`qname`) and their cell barcodes (`annot`). |
| `out/umi_qc.tsv`, `out/export_qc.tsv` | Gene-level processing status and sequence export statistics. |
| `BarcodeMatch/BarcodeMatchIso.txt` | Per-read barcode, UMI, and alignment information. |
| `polish.fq.gz`, `bam/polish.bam` | Adapter-trimmed reads and their alignments. |
| `annotation/` | Processed gene and exon annotations. |
| `input_manifest.tsv`, `run_manifest.rds`, `provenance.txt` | Input fingerprints, run settings, and package versions. |
| `arg.log` | Command-line arguments, when using the command-line interface. |

Each matrix directory contains `matrix.mtx.gz`, `features.tsv.gz`, and
`barcodes.tsv.gz`. **Rows are features; columns are cells.** Both matrix directories
are generated by the annotation step. Setting `--to_isoform FALSE` skips both
matrices; it still writes `iso_count.txt` and the reconstructed FASTQ/BAM.

### Interpreting the counts and reconstructed reads

`iso_count.txt` contains five columns:

| Column | Meaning |
| --- | --- |
| `cell` | Cell barcode. |
| `gene` | Gene identifier. |
| `isoform` | Genomic exon intervals describing the inferred isoform, rather than an annotated transcript ID. |
| `count` | Count after UMI clustering and scatter filtering; it can be fractional. |
| `polyA` | Summary of poly(A) evidence among contributing reads; it can be fractional. |

The matrices are generated after filtering and matching to the supplied
annotation, so their totals and cell sets can differ from `iso_count.txt`.

The reconstructed FASTQ preserves inferred exon structure using **reference
sequence and synthetic quality scores**. It does not retain sample variants,
original base qualities, or original poly(A) tails. Reference `N` bases are
preserved, and each count-table row produces `round(count)` reads. Use the original
reads for variant analysis or inspection of poly(A) tails.

For downstream tools such as [IsoQuant](https://github.com/ablab/IsoQuant), use
`UMI_collapsed.bam` with the read-to-cell mapping in `reads_annot.csv`. The BAM
does not contain `CB` (cell barcode) or `UB` (UMI) tags. Follow the tool's documentation
for your installed version, accounting for the absence of original poly(A) tails
when configuring transcript discovery.

## Reruns and troubleshooting

- **Upgrading or changing inputs:** use a new output directory. To deliberately
  recompute an existing directory, use `--force_rerun TRUE` (R: `force_rerun = TRUE`).
  The same choice is required after an interrupted run.
- **No usable barcodes:** check that reads still contain the adapter/barcode/UMI,
  the library preset matches your chemistry, and the barcode file belongs to the sample.
- **Mapping failure:** check the reported tool log, executable paths, and the
  installed samtools version. LongcellPre requires samtools 1.3 or later.
- **High memory use:** reduce the number of workers. The UMI graph has a default
  limit of 5,000,000 edges per cell/gene (`--max_graph_edges`); this is a graph-size
  limit, not a total RAM limit.
- **Gene processing failure:** the pipeline stops by default. Only use
  `--error_policy continue` if you accept partial results, and inspect
  `out/umi_qc.tsv` for failed genes before interpreting them.

## Version history and support

| Version or tag | Date | Main changes |
| --- | --- | --- |
| 1.3.0 (current `main`) | 2026-10-12 | Combine LongcellPre and Longcellsrc; fix input, splice-correction, and export errors; reduce memory use; add explicit failure reports and rerun checks. |
| v1.2.0 | 2025-03-14 | Add reconstructed-read BAM output; improve annotated isoform matching and rerun controls. |
| v1.1.0 | 2025-01-31 | Add the command-line interface, reconstructed FASTQ export, and demo data. |
| v1.0.0 | 2024-04-30 | First tagged pipeline for barcode assignment, splice correction, and UMI/isoform quantification. |

See [NEWS.md](NEWS.md) for detailed changes and upgrade notes, including untagged
updates and differences between older tag names and package version numbers.

Report problems through [GitHub issues](https://github.com/yuntianf/LongcellPre/issues),
including your package version, command, and error message. Development and
release instructions are in [CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

Fu Y, Kim H, Roy S, et al. **Single cell and spatial alternative splicing analysis
with Nanopore long read sequencing.** *Nature Communications* **16**, 6654 (2025).
[https://doi.org/10.1038/s41467-025-60902-2](https://doi.org/10.1038/s41467-025-60902-2)
