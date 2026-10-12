# Development and release checks

LongcellPre 1.3.0 contains its native backend in `src/`. Do not add Longcellsrc as
an installation dependency. Keep the existing public R functions and CLI arguments
compatible unless a behavior change is intentional and documented in NEWS.md.

For each upgrade, add a dated entry to [NEWS.md](NEWS.md) describing user-visible
changes, fixes, and any migration steps. Keep changes that are not yet released
under an `Unreleased` heading; when releasing, use the version from DESCRIPTION
and the release date, and create a matching Git tag. Update the README's version
summary and examples when the user workflow changes.

## Run CI locally

Use Linux or WSL, Python 3, and micromamba. From the repository root:

```sh
micromamba create -y -n longcell-ci -f .github/ci/environment-linux-64.lock
micromamba run -n longcell-ci python .github/ci/install-source-deps.py
micromamba run -n longcell-ci python .github/ci/run-checks.py
```

The locked environment reproduces the R 4.4.3 / Bioconductor 3.20 dependency set
used for validation. Four additional R source dependencies have fixed versions
and SHA-256 checksums; installation stops if the downloaded contents differ.
The pipeline tools minimap2, samtools and bedtools are included in the environment.
The optional whole-genome BSgenome data packages are not needed for these tests.

Results are written to ignored `ci-output/`, or the directory specified by
`LONGCELL_CI_OUTPUT`. It contains the built source package, installed-test results,
native-equivalence logs, `R CMD check` diagnostics, CLI help, dependency versions,
session information and a checksum. The checks require Longcellsrc to be absent.
All regression tests must execute; missing mapping tools cannot produce a
successful run by silently skipping tests.

CI currently runs on Ubuntu 24.04 with the locked R environment. It performs a
fresh build/install and exact comparisons of 35 UMI graphs, three barcode
thresholds, and 50 isoform tables from synthetic fixtures, together with the
package regression suite. It does not substitute for the complete external demo,
large-data benchmarks, or platform/chemistry validation.

## Package-check diagnostics

R CMD check errors fail the workflow. New or changed warnings also fail it.
`known-warnings.json` records the exact existing Rd-usage warning after whitespace
normalization; unrelated warnings are not suppressed. Removing that warning is
allowed. Do not expand the baseline to bypass a newly introduced problem.
Existing notes, including global-binding and native stdout/stderr diagnostics,
remain visible in the uploaded check log. PDF manual and vignette builds are
not part of this workflow.

When changing C++ exports, regenerate wrappers with `Rcpp::compileAttributes()`.
When changing roxygen comments, regenerate help and NAMESPACE with roxygen2.
Commit generated `R/RcppExports.R`, `src/RcppExports.cpp`, NAMESPACE and help files
alongside their source changes. Keep compiled objects, editor state and CI outputs
out of new commits.

## GitHub and release workflow

1. Create a branch from the existing repository history, review the changes, and
   commit them. The source package belongs at the repository root, including
   `.github/`, `R/`, `src/`, `tests/`, DESCRIPTION and NAMESPACE.
2. Push the branch and inspect **Package CI** in the Actions tab. Local validation
   does not establish that GitHub-hosted execution or repository Actions settings
   are working; the first hosted run confirms that.
3. In the repository's branch protection or ruleset, require the status check
   **Linux R 4.4 / Bioconductor 3.20** before merging. The workflow alone cannot
   prevent a maintainer from merging failing code. These repository settings are
   not changed by the local source work.
4. After a successful run, download `LongcellPre-source-<commit>` from Actions.
   It contains the tested source archive and checksum. Test the complete demo
   before a scientific release and review the unresolved issue limitations.
5. Tag the reviewed commit and publish the checked archive through the maintainer's
   normal GitHub release process. The workflow uses read-only repository access;
   it does not automatically publish releases, push code, or deploy a service.

Artifacts and logs are retained for 30 days. Actions are pinned to full commit
identifiers; update those pins deliberately and rerun the checks when upgrading.
The environment lock tests reproducibility, not compatibility with every newer
dependency. Review and refresh it when expanding supported R/Bioconductor versions.

Workflow building blocks: [setup-micromamba](https://github.com/mamba-org/setup-micromamba),
[checkout](https://github.com/actions/checkout), and
[upload-artifact](https://github.com/actions/upload-artifact).
