# Finalized assembly comparisons — 2026-10-28

**Superseded implementation:** the original Matplotlib/base-level alignment
change below is historical. The current production implementation restores
[pafr/gggenomes in R and mapping-chain synteny](final-comparisons-r-restoration-2026-10-28.md).

Category 7 implementation and local validation are complete for this branch.
Native Crest acceptance remains pending. Categories 2/3/5 remain accepted;
category 6 (unified Markdown reporting) remains open.

## Problems and repairs

| Problem | Production change |
|---|---|
| `all` could compare unrelated species; sample grouping depended on ID suffixes | Pair by explicit sample/taxid metadata; deterministic ordering, unique IDs and same-species pairs only |
| Final comparison inputs did not retain their sequence dictionaries | Pass final FASTA/FAI together; publish a plan, final input inventory, FAI fingerprints and small index copies |
| Runtime unpinned GitHub R installation changed a shared environment | Remove SETUP_PAFR; use separate Python 3.11/minimap2 2.30 and Python 3.11/Matplotlib 3.9.4 environments |
| Missing divergence tags could imply 100% identity | Measure matches / aligned block length directly; missing divergence is NA, with observed-record counts |
| Asynchronous QC could race or hide failures | Stream PAF filtering synchronously under pipefail, then collect validated summaries |
| Malformed/unreadable PAF could become an empty successful plot | Strict parsing, dictionary checks and explicit failure; genuine empty alignments get status tables and explanatory figures |
| Axes hid unaligned final sequences | Construct both axes from final FAI dictionaries, including every sequence meeting the plotting threshold |
| Nested QUAST inputs | Pass the flat assembly and label lists to QUAST_FINAL |
| Disabled comparison sentinel missing | Supply the actual NO_PAIRWISE asset |

Dotplots distinguish forward/reverse orientation. Ribbon plots preserve inversions,
white backgrounds and a shared coordinate scale. Separate alignment/plot tasks
allow plot-only changes to reuse cached alignments. The older R environment and
label remain for accepted pangenome figures; their scripts were not rewritten.

Minimap2 uses a reference-sized single index because multipart indexes produce
incorrect MAPQ ([official minimap2 manual](https://raw.githubusercontent.com/lh3/minimap2/v2.30/minimap2.1)).
Its asm5/10/20 presets remain selectable. These are pairwise assembly comparisons,
not a replacement for the accepted pangenome graph or chromosome harmonization.

## Verification

- Eight Python tests passed: measurements, missing tags, empty/malformed PAF,
  unavailable MAPQ, duplicate summaries, corrupt gzip, dictionary mismatch and plots.
- Twelve actual-workflow cases passed: within-sample/all-species pairing, singleton,
  disabled, empty cohort, no alignments, invalid mode/IDs/species, alignment failure,
  malformed PAF and wrong dictionary. External CLIs were controlled doubles;
  measurement and plotting used the production Python helpers.
- An unchanged normal Nextflow resume cached alignment, plots, QUAST and telomere tasks.
- The full main.nf input-only entry compiled and executed successfully.
- The native-check wrapper passed locally with controlled minimap2, including its
  final-inventory and output reconciliation checks. This does not establish native
  minimap2/conda acceptance.
- Windows plotting tests and visual inspection used the production Matplotlib
  3.9.4 pin. The Linux fixture used Matplotlib 3.10.8 with the available Python 3.14;
  Crest will test the specified Python 3.11/Matplotlib 3.9.4 environment.
- Launcher shell syntax and whitespace checks passed.

## Native Crest check

From the project directory after syncing the repository:

```bash
mkdir -p logs comparisons
sbatch gcl_genome_assembly/scripts/comparisons/run_final_comparison_checks.sbatch \
  comparisons/core-operational-01/results/assembly_eligibility.json \
  comparisons/core-final-comparisons-01
```

This calls the production comparison workflow on accepted final assemblies. It
does not assemble reads or rebuild the graph. Its separate work directory uses
ordinary Nextflow resume. The accepted manifest's FASTA/FAI paths must still exist;
there is no cache recovery or directory search. A normal full pipeline run passes
the final assemblies directly and needs no manually supplied manifest.

Review `results/comparison_acceptance.md` and `results/pairwise_alignments/`.
Return `comparisons/core-final-comparisons-01.review.tar.gz`. The archive includes
plans, inventory/indexes, summaries, statuses, plots, versions, logs and trace;
large PAF files remain on Crest. Expected cohort: ten finalized assemblies and
five within-sample comparisons. The acceptance summary reconciles actual inputs
and comparisons rather than hard-coding those counts.

This check verifies input FAI fingerprints, not an independent rehash of each
whole FASTA. Match fraction describes retained alignments, not whole-genome
identity or biological assembly completeness. Small sequences can be shown by
lowering the plot sequence-size thresholds.

## Remaining operational/report work

Native comparison acceptance, broader dependency/provenance review and the final
combined integration run remain category 7. In particular, existing MitoHiFi
`:master`, FCS image `latest`, and unconstrained dependencies elsewhere still need
review; this change does not claim a globally locked software stack. No blanket
retry policy or cache behavior was changed. The comparison process retries an
external failure once; deterministic evidence/parser/dependency exits terminate.
Optional ignored failures such as pangenome 2D visualization need explicit status
in the category 6 report. Disabled report sections must stay absent.

### Native comparison execution result

`core-final-comparisons-01` completed 13/13 tasks with no failures or retries.
Final input fingerprints and five within-sample pairs reconcile. Sparse retained
alignment evidence and presentation limitations prevent complete visualization
acceptance; see [review and bounded follow-up](core-final-comparisons-01-assessment.md).
No assembly or graph rebuild is indicated.
