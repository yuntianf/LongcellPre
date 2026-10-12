# LongcellPre 1.3.0 (local candidate; integrated native backend)

- Bundle the reviewed Longcellsrc 1.1.1 C++ implementation in LongcellPre.
- Remove the Longcellsrc package dependency and GitHub remote; regenerate Rcpp
  bindings and native registration for LongcellPre.
- Load LongcellPre on parallel workers and record the bundled backend in run
  provenance. Existing pipeline functions and command-line arguments are unchanged.
- Require C++17, zlib and Rcpp when building this single source package.
- Retain the bug fixes and algorithms from the paired 1.2.1/1.1.1 revision below.
- Add GitHub Actions CI with a locked Linux environment, source build/install,
  unskipped regression tests, exact native-output fixtures and checked source
  artifacts. Record the existing documentation warning; reject new/changed warnings.

# LongcellPre 1.2.1 (local candidate; paired with Longcellsrc 1.1.1)

- Validate untrimmed FASTQ input and chemistry-specific barcode lengths; accept
  Cell Ranger suffixes and explain mislabeled compression (#5, #6, #7, #13).
- Preserve a matrix dimension for one or zero internal splice columns (#8).
- Export within-gene batches, use the mapping FASTA for reference reconstruction,
  and report reference Ns per gene (#9, #10, #11).
- Fail visibly on gene-level errors; write gene QC, pass a configurable native
  graph edge limit to workers, and remove silent oversized-cell skips (#10–12).
- Check external exit statuses with pipefail, quote paths, validate BAMs and
  stage output replacements; retain modern samtools syntax (#14).
- Record input hashes/settings and prevent reuse of changed or incomplete runs.
- Longcellsrc: reusable barcode postings index, bounded UMI score cache, linear
  alignment workspace, streamed isoform comparisons, 64-bit-capable weights,
  and bounds-safe read extraction (#4, #7, #10–12).
- Full demo equality and performance evidence accompany this local candidate.
  Large patient datasets and Figure 2 have not been reproduced; PARSE (#16)
  is deferred pending the commenter's chemistry details.

## Changes retained from the first local review

- Generate FASTQ qualities at the exact reconstructed sequence length, including
  reads longer than 50,000 bases (#15).
- Qualify sequence concatenation through Biostrings and validate exon bounds (#10).
- Stream the pipeline's reconstructed FASTQ and barcode annotations by gene;
  retain the existing in-memory return mode for direct callers (#10, #11).
- Handle counts that round to zero; reject invalid counts and quality strings.
- Rebuild collapsed BAM output when its FASTQ was regenerated or its index is missing.
- Avoid a read-count-sized temporary vector when assigning genes to workers.
- Declare missing imported dependencies and the Longcellsrc GitHub remote (#1).
- Create a new output directory before writing the command-line argument log.
- Explain sample-specific whitelists, output semantics, version-sensitive
  IsoQuant integration, and current protocol limits (#2, #7, #9, #13, #16, #17).
- Add focused export and pipeline regression tests.
