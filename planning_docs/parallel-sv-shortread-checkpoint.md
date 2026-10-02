# Work while final HiFi mapping runs

2026-10-02. Comparison-only additions; production defaults unchanged. No pipeline, Python, R or unit tests executed locally. User runs these checks on the cluster. Do not restart the active mapping run.

## 1. SV advisory QC and simplified plots

The existing diagnostics runner now writes `sv_qc.tsv`, `sv_qc_settings.json` and `sv_qc.md`. Every native selected-class event remains present. It checks sequence-boundary/N proximity and direct alignment-child span coverage (unioned, clipped, separately in each assembly). Default flags use 1 kb proximity and 0.5 assigned-span fraction; these are exploratory review thresholds, not calibrated biological acceptance criteria. The Python helper exposes both as arguments for sensitivity checks.

`NO_FLAGS_IN_IMPLEMENTED_CHECKS` is deliberately not PASS, high-confidence, or validated. Read evidence, alignment identity, alternative placements, repeat ambiguity and adjacency validation remain unassessed. CPG/CPL/TDM typically have no direct AL children and cannot be graded by that measure. Failure of the helper fails the diagnostic job; raw calls are never rewritten.

Removed the wide lower panels. Default plots select the largest event per class and advisory status. The optional third runner argument is a comma-separated list of exact event IDs (e.g. INVDP5822,INV2308); unknown IDs fail explicitly. Selection is diagnostic, not a call filter. Existing file numbering can change; consult `largest_candidates.tsv` for IDs.

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_diagnostics.sbatch \
  "$PWD/comparisons/environments/sv-diagnostics-r" \
  comparisons/sv-syri-retry-20261002-070635-1503156
```

Return the `comparisons/sv-diagnostics-<timestamp>-<job>.tar.gz` archive and `logs/sv-diagnostics-<job>.out`. No new alignments or SyRI run required. Acceptance: all native event IDs retained; missing metrics remain unavailable; review a selection of flagged and unflagged intervals before adopting thresholds. Added tests cover interval orientation/union, cross-line N runs, missing evidence, out-of-range coordinates and read-subset integrity.

## 2. Short-read smoke preparation

The subset is the first 2 million paired records, not a random biological sample. It avoids scanning the entire input. Mate names, four-line FASTQ structure and sequence/quality lengths are checked in the selected prefix; the rest of the source is not validated. Input/output provenance and uncompressed subset hashes are recorded. Existing destination directories are refused. A preparation failure leaves a diagnostic partial directory without `subset.json`; use a new destination after fixing it.

From the cluster project root, after syncing:

```bash
mkdir -p logs comparisons
sbatch --partition=normal --time=01:00:00 --mem=2G --cpus-per-task=1 \
  --output=logs/shortread-prepare-%j.out \
  --wrap='python3 gcl_genome_assembly/scripts/comparisons/prepare_shortread_smoke.py raw_fastq/ssl/Sde-CMat_061_1.fq.gz raw_fastq/ssl/Sde-CMat_061_2.fq.gz comparisons/shortread-smoke-01 --sample Sde-CMat_061_smoke --taxid 373251 --pairs 2000000'
```

After preparation succeeds:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/shortread-smoke-01 validate
```

After validation succeeds:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/shortread-smoke-01 assembly
```

The test has separate launch directories (including validation), work and results, and does not use the main assembly resume session. SPAdes uses k=21,33, coverage cutoff off, no isolate mode, and `--only-assembler` (error correction deliberately untested). Redundans retains all three stages with one iteration. Both tasks now inherit production resource, queue, retry and time settings. SPAdes uses 64 CPUs with genome/input-dependent memory; Redundans uses 24 CPUs and 700 GB. Both retain their normal 96-hour task limit; the handler has a four-day limit. These are test settings, not performance promises. Other enabled preparation/decontamination tasks retain their normal resource requirements and shared database behavior.

This first smoke run has QC and post-assembly plots disabled. It tests FASTQ preparation, short-read routing, conditioning, singleton bypass/finalization, eligibility and the always-generated assembly summary. It does not validate the complete QC/reporting stack, assembly quality, SPAdes error correction or every optional branch. Preserve low-quality flags. If the subset yields no useful contigs or skips a downstream operation, that operation is not validated: use a targeted fixture or a larger separate subset. A later final_only/report checkpoint will exercise the remaining reporting path.

SPAdes option reference: https://ablab.github.io/spades/running.html . No production parameter file includes the smoke config.

## 3. Collect short-read results

After completion or failure, from the project root:

```bash
testdir=comparisons/shortread-smoke-01
find "$testdir/results" "$testdir/logs" -type f \
  \( -name '*.tsv' -o -name '*.json' -o -name '*.md' -o -name '*.log' -o -name '*trace*' \) \
  > "$testdir/review.files"
printf '%s\n' "$testdir/subset.json" "$testdir/samples.csv" >> "$testdir/review.files"
tar -czf "$testdir/review.tar.gz" -T "$testdir/review.files"
```

Return that archive and the corresponding `logs/shortread-smoke-<job>.out` and preparation log. For a failed task also return its `.command.sh`, `.command.err` and `.command.out` from the work directory named in the Nextflow log. Keep reads/FASTA/BAM/work directories on the cluster. Check the trace for actual SPAdes/Redundans resource requests and execution rather than assuming the overrides matched.

## Updated launch choice: full reads and production resources

The user prefers getting a usable integration assembly over minimizing resource requests; large nodes are available. Use the full original FASTQs with the reduced assembly settings above, rather than the 2-million-pair prefix. The following supersedes the subset commands for the next run. No input copying is performed, and production defaults remain unchanged. The full dataset can still take substantial time; no runtime or assembly success guarantee is implied.

```bash
python3 gcl_genome_assembly/scripts/comparisons/prepare_shortread_smoke.py \
  raw_fastq/ssl/Sde-CMat_061_1.fq.gz \
  raw_fastq/ssl/Sde-CMat_061_2.fq.gz \
  comparisons/shortread-integration-01 \
  --sample Sde-CMat_061 --taxid 373251 --full-input
sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/shortread-integration-01 validate
# After validation succeeds:
sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/shortread-integration-01 assembly
```

Use `testdir=comparisons/shortread-integration-01` in the collection commands above. This remains an isolated first run. When combining with HiFi inputs, plan the shared launch/work/settings explicitly: this separate work directory does not automatically supply a cache to the main run. Do not promise cross-run reuse before establishing that handoff.