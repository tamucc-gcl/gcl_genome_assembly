# R synteny restoration — 2026-10-28

This supersedes the Matplotlib implementation described in the earlier comparison
contracts. The user's plotting preference is R/ggplot2 wherever practical, now
recorded in AGENTS.md for future changes.

## What changed and why

The previous repair combined legitimate input/error/provenance fixes with two
unnecessary behavioral changes: replacing pafr/gggenomes with custom Python plots,
and adding minimap2 `-c` for base-level alignment. The latter changed the records
subject to the existing 150 kb threshold. Sparse plots were a regression, despite
successful task exits. The length threshold itself had not changed.

The production comparison workflow again uses minimap2 mapping chains without
`-c` or `--eqx`, pafr/ggplot2 dotplots and gggenomes/ggplot2 riparian plots.
The single-index guard, metadata/species pairing, final FASTA/FAI validation,
strict failures and separate alignment/plot task caching remain. There is no
runtime GitHub package install and no Python plotting fallback. The retired
Matplotlib pairwise renderer is removed.

The R environment specifies R 4.5, ggplot2 4.0.1, pafr 0.0.2, gggenomes 1.2.0,
argparse 2.2.5 and dplyr 1.1.4. This combination was solved and installed locally
on Linux; actual R packages executed the rendering checks. Python 3.11 is present
because R argparse invokes Python. The separate Python evidence helper generates
tables and compressed PAFs, never plots. Existing pangenome R environments remain
unchanged.

## Evidence and figure behavior

- Each pair retains the unfiltered mapping PAF, the configured filtered PAF, QC,
  and a small retention table. Review archives exclude both large PAFs.
- Retention reports unfiltered, MAPQ-only, 10 kb, 50 kb, 150 kb and configured
  length thresholds (duplicate length values are combined). Reference and query
  covered bases are interval unions, so overlapping mappings are counted once.
- Coverage denominators are the complete final assemblies, not only plotted
  chromosome-size sequences. These are mapping-chain coordinate spans, not
  proof of base-by-base alignment through the span.
- The QC field is `approx_match_fraction`, explicitly marked
  `measurement_basis=mapping_chain_approximate`. It is not whole-genome identity
  or an exact base-level identity estimate. Missing divergence remains NA.
- Plot subtitles show the effective upstream/plot length and MAPQ thresholds.
- Dotplots retain pafr's mapping geometry, adapting its concatenated offsets to
  the complete final FAI so unaligned sequences remain visible. No artificial
  alignments are added. Forward/reverse mapping geometry is tested.
- Chromosome names label modest-size axes. Fragmented inputs use numeric Mb axes
  without hundreds of separator lines. Riparian links are colored by reference
  sequence, with neutral sequence bars. Labels/legends are limited for dense
  fragmented inputs. Internal sequence IDs include bin prefixes to prevent
  accidental collisions across assemblies; display labels use real final names.
- White backgrounds and 300 dpi exports apply to both figure types.

The pinned pafr version emits ggplot2 deprecation warnings for its old aes/line
syntax. They do not alter output; malformed PAF/gzip/index failures still stop.
No upstream package source was patched.

## Verification

Eight evidence unit tests passed, including prefilter preservation, empty input,
missing divergence, MAPQ handling, overlapping interval union and dictionary errors.
Twelve Nextflow workflow scenarios passed using controlled external minimap2,
QUAST and tidk commands **and real R rendering**. These cover species/sample
selection, disabled/empty/singleton routes, alignment failure/retry, malformed PAF
and wrong dictionaries. The controlled minimap2 asserts that `-c`/`--eqx` are absent.

The R regression script checks actual built plot coordinates, full FAI axes,
inversions, empty/malformed/corrupt inputs and fragmented layouts. Figures were
visually inspected. Both the full-cohort and sample-selected native-check wrappers
passed locally with real R and controlled mapping. The full main.nf input-only
entry compiled/executed. An unchanged normal resume cached comparison alignments,
R plots and other FINAL_VIZ tasks. No expensive assembly or graph run was launched.

Real CTlk/CMat mapping coverage recovery remains a Crest acceptance question;
the local controlled mapping tests do not establish that biological result.

## First Crest check

After syncing, run from the project directory:

```bash
mkdir -p logs comparisons
sbatch gcl_genome_assembly/scripts/comparisons/run_final_comparison_checks.sbatch \
  comparisons/core-operational-01/results/assembly_eligibility.json \
  comparisons/core-final-comparisons-r-01 \
  --comparison_sample Sde-CMat_203
```

This selects the sample by manifest metadata, not a production hard-coded sample
list. It runs the actual production comparison workflow on one finalized pair,
without rebuilding assemblies, QC mappings or the graph. Read
`results/comparison_acceptance.md`, the pair's `.retention.tsv` and both R figures.
Return `comparisons/core-final-comparisons-r-01.review.tar.gz`.

After that check, omit `--comparison_sample` and use a separate output directory
to check all five pairs concurrently. This restoration necessarily invalidates
the changed pairwise alignment/plot tasks; earlier `-c` outputs cannot serve as
the restored mapping chains. Ordinary Nextflow caching is otherwise unchanged.

## Plotting migration scope

These optional pairwise plots are now R-only. Existing Python figures elsewhere
have not silently been rewritten. The reporting pass should migrate the active
misassembly panels/contact plots, mitochondrial circular figure, coverage books
and optional compartment/TAD figures where practical, keeping their accepted
evidence calculations intact. The earlier investigation/legacy chimera scripts
should be checked for active use before porting or retirement. Category 6 owns
that figure migration; category 7 retains native comparison acceptance, remaining
dependency audit and final combined integration.

### Native R comparison pilot accepted

`core-final-comparisons-r-01` passed all five tasks and restored CMat chromosome
synteny: 95.88%/96.44% retained chain-span coverage. Final input dictionaries are
unchanged from the sparse run. Both R figures were inspected. No additional code
fix is indicated; the five-pair cohort check is next. See
[pilot assessment and full-cohort command](core-final-comparisons-r-01-assessment.md).
