# Version history

This log describes user-visible changes, newest first. Version 1.3.0 is the
current package version on `main`; it does not yet have a `v1.3.0` tag. Earlier
entries use the repository's existing tag names and dates. Some historical tags
differ from the version recorded inside the package; those differences are noted
below.

## 1.3.0 — 2026-10-12

### Installation and compatibility

- Combine LongcellPre and the Longcellsrc native backend into one R package.
  Installing or loading Longcellsrc separately is no longer necessary.
- Preserve the `RunLongcellPre()` entry point and existing command-line arguments.
  Parallel workers load the integrated package automatically.
- Correct required dependency declarations. Building from source requires a
  C++17 compiler, Rcpp, and zlib.

### Input and processing fixes

- Accept sample barcode lists with Cell Ranger suffixes such as `-1`, validate
  barcode lengths against the selected chemistry, and report empty or unusable
  adapter/UMI input clearly.
- Detect FASTQ compression from file contents, including demo files with a
  misleading `.gz` extension.
- Fix splice-site correction for genes with only one or no internal splice
  columns.
- Reduce temporary memory use in barcode matching, UMI scoring, isoform
  comparisons, and distribution of genes across workers.
- Fix native read-extraction bounds and large integer-weight handling.
- Replace silent omission of oversized or failed cell/gene groups with explicit
  errors and gene-level QC. Add `max_graph_edges` and an opt-in
  `error_policy = "continue"` for partial results.

### Output and rerun fixes

- Generate FASTQ quality strings at the correct sequence length, including reads
  longer than 50,000 bases. Handle counts that round to zero and reject invalid
  counts or quality strings.
- Write reconstructed FASTQ and read-to-cell annotations in bounded batches to
  reduce export memory use. Add `batch_reads`, `batch_bases`, and `export_qc.tsv`.
- Use the supplied mapping FASTA for sequence reconstruction. Preserve reference
  `N` bases and report missing chromosomes or invalid exon coordinates.
- Check external tool failures and BAM validity, handle paths containing spaces,
  and rebuild collapsed BAMs when their FASTQ changes or their index is missing.
- Create the output directory before writing the command-line argument log.
- Record inputs, settings, and package versions. Reject reuse of incompatible or
  incomplete outputs unless the user explicitly requests a full rerun.
- Add `keep_intermediates` to retain barcode and coverage intermediate files.

### Documentation and validation

- Rewrite the README around installation, inputs, working command-line/R
  examples, the demo, output interpretation, and troubleshooting. Add this
  version history and separate contributor instructions.
- Clarify that exported reads contain reference-derived sequence and synthetic
  qualities, without original variants or poly(A) tails.
- Add regression tests and GitHub Actions checks for package builds, installation,
  and agreement with saved native-output fixtures.
- Validate three complete demo runs: count values and labels, reconstructed read
  sequences/qualities, and compared alignment records matched the saved upstream
  results. This establishes agreement for the demo, not every dataset.

### Upgrading

Install LongcellPre alone and use a fresh output directory for each analysis
after upgrading. Public pipeline calls remain compatible; code that calls
`Longcellsrc::` helper functions directly is not redirected to LongcellPre.

The earlier **1.2.1 development snapshot**, paired with Longcellsrc 1.1.1, was not
a tagged release. Its fixes are included in 1.3.0 above. PARSE combinatorial
barcode support remains unavailable.

## Updates after v1.2.0 — 2025-07-29 to 2025-09-23

These changes reached `main` without a new release tag and are retained in 1.3.0:

- Correct command-line argument types, UMI-length handling, sequence extraction,
  and adapter-score processing.
- Improve parallel work distribution, read correction, and read filtering.
- Separate UMI counting from annotated isoform assignment and fix the associated
  pipeline and isoform-imputation errors.
- Write genomic coordinates without scientific notation when constructing
  isoform strings.

## v1.2.0 — 2025-03-14

- Add alignment and BAM output for UMI-collapsed, reconstructed reads, together
  with a comma-separated read-to-cell annotation file.
- Improve matching of inferred isoforms to annotated transcripts and correct
  UMI-count and alignment-reading errors.
- Add rerun controls for existing barcode/alignment results.
- Adjust splice-site filtering defaults and fix parallel execution and
  small-reference handling.
- Update dependency declarations and usage documentation.

[Changes from v1.1.0](https://github.com/yuntianf/LongcellPre/compare/v1.1.0...v1.2.0).
The package version recorded at this tag is 1.2.0.

## v1.1.0 — 2025-01-31

- Add the command-line entry point and document pipeline parameters.
- Add reconstructed FASTQ export and the read-to-cell interface for downstream
  transcript analysis.
- Add automatic reference-genome loading for the original BSgenome workflow.
- Improve erroneous-mapping correction and fix processing of single-isoform genes.
- Add demo data and update installation/dependency documentation.

[Changes from v1.0.0](https://github.com/yuntianf/LongcellPre/compare/v1.0.0...v1.1.0).
Although the tag is named v1.1.0, its DESCRIPTION records package version 1.2.0.

## v1.0.0 — 2024-04-30

- First tagged version of the R pipeline, covering cell-barcode assignment,
  read mapping, splice-site correction, UMI counting, and annotated isoform counts.
- Provide the R workflow and configurable library presets.

[Tagged source](https://github.com/yuntianf/LongcellPre/tree/v1.0.0).
Although the tag is named v1.0.0, its DESCRIPTION records package version 1.1.0.
