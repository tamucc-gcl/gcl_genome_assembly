# Contamination evidence and reviewed cleaning

## Implementation status

The [assembly pipeline](#run-post-fcs-review-in-the-assembly-pipeline) can generate
post-FCS evidence and apply reviewed removals before scaffolding. It can also
stop at the contig stage. The standalone tools support database setup and report
rebuilds. In the assembly pipeline, `--run_contamination_evidence` enables the
new report and BlobToolKit dataset.

The supported script for this foundation is `py_scripts/contamination_review.py`.
It requires Python 3.11 or later and has no third-party Python dependencies.
It processes one sequence at a time rather than loading a cohort into memory.

## Prepare an editable review sheet

For the real CTlk test, start with [Run the bounded CTlk pilot on Crest](#run-the-bounded-ctlk-pilot-on-crest).
That workflow prepares the review sheet for you. The standalone command below
uses example paths and assembly IDs. Substitute an existing FASTA and a writable
output directory when running it on another assembly.

```bash
python3 gcl_genome_assembly/py_scripts/contamination_review.py prepare \
  --fasta /project/assembly/contigs.fasta \
  --assembly sample_hap1 --stage post_fcs \
  --out /project/review/sample_hap1
```

Outputs are `source.json`, `sequence_metrics.tsv` and `review.tsv`. The template
has one row per source sequence, all selected `NO`. These are **not contamination
calls or suggested removals**. Detection evidence will be added by later tools.
GC fraction uses A/C/G/T bases as its denominator; N and other IUPAC ambiguity
counts are separate. All-ambiguous sequences have GC `NA`.

Review columns:

| Column | Meaning |
| --- | --- |
| `action_id` | Unique within an assembly; use a new ID when adding a manual interval |
| `assembly` | Assembly identity; a cohort sheet may contain multiple assemblies |
| `stage` | Exact coordinate-stage label, initially `post_fcs` |
| `source_sha256` | SHA256 of the exact input FASTA file, including headers/formatting |
| `sequence` | Original FASTA identifier (first header token) |
| `action` | `REMOVE_SEQUENCE` or `REMOVE_INTERVAL` |
| `start0`, `end0` | 0-based, half-open source coordinates |
| `selected` | Explicit `YES` or `NO`; no automatic removal |
| `reviewer`, `reason` | Required for selected actions |

For whole-sequence removal, coordinates must be `[0, sequence length)`. To remove
bases 101 through 200 in familiar 1-based inclusive coordinates, use
`REMOVE_INTERVAL`, `start0=100`, `end0=200`. Copy the source identity/stage/checksum
when adding an interval row; give it a new action ID. Never use coordinates from
a later renamed, scaffolded or finished FASTA against this sheet. Additional
evidence columns are permitted and retained in the applied-action audit.

## Apply a reviewed sheet

Apply the sheet to the exact FASTA used for evidence. For a post-FCS
report, the assembly Nextflow interface accepts `--contamination_review`
as described below. The standalone command is also available.

```bash
python3 gcl_genome_assembly/py_scripts/contamination_review.py apply \
  --fasta /project/assembly/contigs.fasta \
  --assembly sample_hap1 --stage post_fcs \
  --review /project/review/reviewed.tsv \
  --out /project/cleaned/sample_hap1
```

- `retained.fasta`: retained sequence. An internal deletion produces separate
  flanks, never a concatenation across the removed interval. End removals trim.
- `removed.fasta`: all removed sequence; intentionally empty with no selections.
- `coordinate_map.tsv`: original intervals, output piece IDs/coordinates,
  retained/removed status, piece SHA256 and removal action IDs.
- `applied_actions.tsv`: selected rows for this assembly, including added evidence.
- `verification.json`: source/review/output checksums, action count and exact
  source/retained/removed bp accounting.

Changed pieces have deterministic IDs containing their disposition and source
coordinates. Unchanged records retain their headers/sequence. With no selected
actions the retained FASTA is byte-identical to the source. All remaining pieces,
including short ones, are preserved; there is no implicit length-based discard.

Unknown sequences, wrong source checksums/stages, invalid selections, missing
review reasons, overlapping removals and output-ID collisions fail. Adjacent
interval removals are allowed. Removal of the entire assembly is rejected.
Existing output directories are not overwritten. Results are built in a temporary
sibling directory and exposed only after successful validation/application.

The caller validates cohort membership before parallel assembly application;
this per-assembly tool deliberately ignores other assemblies' rows. It does not
carry BUSCO annotations or read mappings onto newly split sequences.

## Prepare parallel protein-search inputs

```bash
python3 gcl_genome_assembly/py_scripts/contamination_review.py search-input \
  --fasta /project/assembly/contigs.fasta \
  --assembly sample_hap1 --stage post_fcs \
  --window-bp 100000 --chunk-bp 5000000 --min-acgt-bp 100 \
  --out /project/search_inputs/sample_hap1
```

`windows.tsv` maps neutral query IDs to source coordinates. `chunks.tsv` lists
bounded search FASTAs and their query counts/bp; these files can be searched in
parallel. Non-overlapping windows cover the entire assembly, including short
terminal windows. Windows below the explicit minimum ACGT count remain in the
table as `insufficient_acgt`, rather than being mistaken for a searched no-hit.
All-ambiguous assemblies emit zero search chunks with explicit accounting.

These are preparation defaults, not validated taxonomic-search sensitivity
settings. Fixed window boundaries and fragmented proteins must be evaluated in
the search pilot. The manifest records source identity, parameters and the window
boundary policy. No nucleotide follow-up selection or biological classification
is performed by this preparation command.

## Pilot the reference-proteome database build on Crest

Use a small combined FASTA of canonical UniProt reference proteome proteins for
the first run. Keep the original UniProt headers, including their `OX=` taxon IDs.
A gzipped FASTA is accepted. Use an absolute path on Crest and a pilot release
name distinct from the full collection. From the assembly project directory:

```bash
module load nextflow/23.10.1 miniconda3
pipeline="$PWD/gcl_genome_assembly"
pilot=/absolute/path/to/pilot_reference_proteomes.fasta.gz
sha=$(sha256sum "$pilot" | cut -d' ' -f1)
nextflow -c "$pipeline/nextflow.config" run \
  "$pipeline/scripts/comparisons/contamination_protein_db.nf" -profile slurm \
  --outdir "$PWD/comparisons/contamination-protein-pilot-01" \
  --contamination_db_root /work/birdlab/databases \
  --proteome_source "$pilot" --proteome_sha256 "$sha" \
  --proteome_release pilot-01 --proteome_format fasta
```

The two mandatory stored outputs are `reference_proteomes.dmnd` and
`protein_database.json`, under
`diamond/uniprot_reference_proteomes/pilot-01` in the selected database root.
The builder uses the shared `ncbi_taxonomy` directory by default and the
existing taxonomy downloader if nodes or names are absent. An explicit
`--taxonomy_dir` selects another location. The manifest records input and
database checksums, taxonomy file checksums, DIAMOND version, sequence totals,
and the taxid recovered in a test search. Repeat the same command to exercise
Nextflow's native `storeDir` reuse.

For the full collection, supply a release-specific local archive or HTTPS URL
with its expected SHA256 and set `--proteome_format tar`. The archive reader
streams canonical `UP…_<taxid>.fasta[.gz]` members and skips CDS and additional
isoform members. Budget temporary space for the archive, combined search inputs
and built index. When taxonomy is intentionally refreshed, use a new storage
release label to preserve the previous database and its recorded snapshot.

## Managed protein databases in the main workflow

Choose either an existing `--contamination_protein_db_dir` bundle or managed
download/build controls:

```bash
--contamination_protein_db_root /absolute/database/root \
--contamination_protein_release 2026_04
```

The managed path uses the official UniProt reference-proteome archive. The release
must match the publisher's HTTPS `RELEASE.metalink`, which pins the archive name,
compressed size and MD5. The completed download is checked against those pins,
also hashed with SHA256, and retained separately from the index build. Later
build attempts reuse this verified download. The build checks its SHA256 again,
streams canonical proteins, constructs taxonomy mappings and performs a real
self-search to check taxid lookup. A matching existing bundle is resolved offline;
the evidence workflow still checks index and taxonomy checksums before searches.
No annotation pipeline is required.

The user-selected root contains:

```text
downloads/uniprot_reference_proteomes/<release>/source/archive
downloads/uniprot_reference_proteomes/<release>/source/source.json
diamond/uniprot_reference_proteomes/<release>/reference_proteomes.dmnd
diamond/uniprot_reference_proteomes/<release>/protein_database.json
```

The last two files form the reusable bundle. Compatible bundles from another
pipeline can be supplied directly; a bare DIAMOND file without the taxonomy and
provenance manifest does not meet this contract. A changed source/checksum under
the same release label is rejected. Use a new label for a new source. Reuse does
not refresh shared taxonomy files; an outdated/incompatible taxonomy must be
updated deliberately rather than dropping unmapped protein records.

Managed protein builds default to a separate reusable taxonomy snapshot at
`<root>/ncbi_taxonomy/<release>`, including `nodes.dmp`, `names.dmp` and
`merged.dmp`. This does not refresh the assembly/FCS taxonomy snapshot. Set
`--contamination_protein_taxonomy_dir` to reuse a compatible shared snapshot.
When reusing an index through `--contamination_protein_db_dir`, the same option
can select its matching taxonomy. All three files control taxonomy `storeDir`
reuse; the protein manifest records their checksums.

Custom local or HTTPS sources also require `--contamination_protein_source`,
`--contamination_protein_sha256` and `--contamination_protein_format fasta|tar`.
This supports archived releases and bounded test FASTAs. Temporary prepared FASTA
and taxon maps use fast gzip compression; the sequences and search settings are
unchanged. A tiny real build confirmed the identical DIAMOND database checksum
before and after this compression-only change.

### Build before benchmarking the full collection

This helper calls the same managed setup workflow as `main.nf`, without starting
the cohort searches:

```bash
pipeline="$PWD/gcl_genome_assembly"
run="$PWD/comparisons/contamination-full-protein-01"
nextflow -c "$pipeline/nextflow.config" run \
  "$pipeline/scripts/comparisons/contamination_managed_protein_db.nf" -profile slurm \
  -w "$run/work" --outdir "$run/build" \
  --contamination_protein_db_root "$run/db" --contamination_protein_release 2026_04
```

The 2026_04 catalog inspected on 20 November 2026 lists 355,492,909,566 compressed
bytes. This is a one-time download; also budget space for the prepared protein
FASTA/taxon map and final index. The final index size and full-database search
runtime must be measured; the small eight-species pilot does not establish them.

### Protein storage/throughput benchmark

The standalone benchmark selects twelve retained 5 Mb chunks across assemblies
and positions from a production trace. It runs the accepted DIAMOND 2.1.14
settings twice on one exclusive node: shared database and copied scratch database.
Both use scratch for temporary files, eight concurrent searches and eight threads
per search. Sensitivity, frameshift handling, range culling, output fields and
reporting thresholds remain unchanged. Normalized hit files must match exactly.

```bash
nextflow -c "$pipeline/nextflow.config" run \
  "$pipeline/scripts/comparisons/contamination_protein_benchmark.nf" -profile slurm \
  -w "$run/benchmark-work" --outdir "$run/benchmark" \
  --benchmark_project "$PWD" \
  --benchmark_trace "$PWD/comparisons/contamination-pipeline-01/results/pipeline/pipeline_trace.txt" \
  --benchmark_protein_db_dir "$run/db/diamond/uniprot_reference_proteomes/2026_04" \
  --benchmark_node crest-c018
```

Use a node with enough scratch for the index plus 100 GiB per simultaneous search;
`--benchmark_node` is optional. The temporary copy is removed on normal completion
or a handled failure. The benchmark records commands, individual runtimes and
peak memory, index staging time, and rough cohort extrapolations by chunk count
and query span at the tested concurrency. It compares warm/cache-influenced runs,
so it cannot establish cold-cache storage performance or scaling across nodes.
Read `<outdir>/benchmark/benchmark.json` and `<outdir>/queries/selection.json`.

The previous panel had 2,088 chunks, 39.9 seconds median, and 36.58 aggregate task
hours (4.57 ideal hours at eight-way concurrency). Full-database forecasts should
use the new measurements. In particular, DIAMOND block size affects results, and
iterative sensitivity is intended for best-hit searches; neither is silently
changed for this multi-hit contamination screen. See the
[DIAMOND performance options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options).

## Run the known-answer search pilot on Crest

From the assembly project root, prepare a fixture and a fresh protein database.
The fixture plants host, bacterial, and mixed sequences, and includes no-hit and
all-ambiguous controls. The bacterial coding sequence is on the reverse strand.

```bash
module load nextflow/23.10.1 miniconda3
pipeline="$PWD/gcl_genome_assembly"
pilot="$PWD/comparisons/contamination-evidence-pilot-01"
taxonomy=/work/birdlab/databases/ncbi_taxonomy
mkdir -p "$pilot"
python3 "$pipeline/scripts/comparisons/prepare_contamination_pilot.py" "$pilot/fixture"
sha=$(sha256sum "$pilot/fixture/proteins.faa" | cut -d' ' -f1)
nextflow -log "$pilot/db.nextflow.log" -c "$pipeline/nextflow.config" run \
  "$pipeline/scripts/comparisons/contamination_protein_db.nf" -profile slurm \
  -w "$pilot/work" --outdir "$pilot/db_run" --contamination_db_root "$pilot/db" \
  --taxonomy_dir "$taxonomy" --proteome_source "$pilot/fixture/proteins.faa" \
  --proteome_sha256 "$sha" --proteome_format fasta --proteome_release known-answer-01
```

Prepare the tiny local BLAST database and the test BAM once:

```bash
conda create -y -p "$pilot/fixture_tools" -c conda-forge -c bioconda \
  'blast=2.16.0' 'samtools=1.21'
conda run -p "$pilot/fixture_tools" makeblastdb -dbtype nucl -parse_seqids \
  -in "$pilot/fixture/nucleotide.fasta" \
  -taxid_map "$pilot/fixture/nucleotide_taxids.tsv" -out "$pilot/fixture/nt"
conda run -p "$pilot/fixture_tools" samtools view -b \
  -o "$pilot/fixture/reads.bam" "$pilot/fixture/reads.sam"
```

Launch the evidence workflow:

```bash
nextflow -log "$pilot/evidence.nextflow.log" -c "$pipeline/nextflow.config" run \
  "$pipeline/scripts/comparisons/contamination_evidence.nf" -profile slurm \
  -w "$pilot/work" --outdir "$pilot/results" \
  --assembly known_answer --fasta "$pilot/fixture/assembly.fasta" \
  --taxid 7955 --taxonomy_dir "$taxonomy" \
  --protein_db_dir "$pilot/db/diamond/uniprot_reference_proteomes/known-answer-01" \
  --nucleotide_db "$pilot/fixture/nt" --coverage_bam "$pilot/fixture/reads.bam" \
  --window_bp 2000 --chunk_bp 2000
```

Open `results/review/README.md`. Compare sequence labels, nucleotide targets,
and the unselected intervals with `fixture/expected.json`. The planted
non-host interval in `mixed` is **[960, 1560)**. The report embeds the R-rendered
GC/coverage plot with span-weighted margins and a cumulative sequence plot. The
BlobToolKit dataset in `results/blobdir` includes MEGABLAST taxonomy where
available, protein evidence as fallback, and the combined review category.
All-ambiguous sequences remain in the evidence tables and are listed in
`blobdir_exclusions.tsv` because BlobToolKit requires measurable GC.
Coverage uses primary reads at MAPQ ≥20, omitting duplicates and QC failures.
Without coverage the static plot shows GC vs sequence length.

## Run the bounded CTlk pilot on Crest

Run from the assembly project root. The helper uses the production trace to
resolve the retained hap1 FCS output and matching contig QC coverage table.
It selects whole contigs across length and GC strata, up to 30 Mb total.
The selection deliberately includes GC extremes, so its counts describe the
pilot subset. The selection manifest records each original contig and the
reason it was included. The protein panel covers three fish, an arthropod,
two fungi and two bacteria. Its source URLs, UniProt release and checksums are
recorded separately.

```bash
module load nextflow/23.10.1 miniconda3
(
  set -euo pipefail
  pipeline="$PWD/gcl_genome_assembly"
  pilot="$PWD/comparisons/contamination-ctlk-pilot-01"
  taxonomy=/work/birdlab/databases/ncbi_taxonomy
  prepare="$pipeline/scripts/comparisons/prepare_contamination_real_pilot.py"
  mkdir -p "$pilot"

  ( [[ -f "$pilot/inputs/selection.json" ]] || python3 "$prepare" contigs \
      --assembly Sde-CTlk_104_hap1 --use-trace-coverage --out "$pilot/inputs" ) &
  contigs_pid=$!
  ( [[ -f "$pilot/panel/panel.json" ]] || python3 "$prepare" panel --out "$pilot/panel" ) &
  panel_pid=$!
  wait "$contigs_pid"
  wait "$panel_pid"

  sha=$(sha256sum "$pilot/panel/proteins.faa.gz" | cut -d' ' -f1)
  nextflow -log "$pilot/db.nextflow.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_protein_db.nf" -profile slurm \
    -w "$pilot/work" --outdir "$pilot/db_run" --contamination_db_root "$pilot/db" \
    --taxonomy_dir "$taxonomy" --proteome_source "$pilot/panel/proteins.faa.gz" \
    --proteome_sha256 "$sha" --proteome_format fasta --proteome_release ctlk-panel-01

  nt=()
  if [[ -n "${CONTAMINATION_NT_DB:-}" ]]; then
    nt=(--nucleotide_db "$CONTAMINATION_NT_DB")
  fi
  nextflow -log "$pilot/evidence.nextflow.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_evidence.nf" -profile slurm \
    -w "$pilot/work" --outdir "$pilot/results" \
    --assembly Sde-CTlk_104_hap1 --taxid 373251 --taxonomy_dir "$taxonomy" \
    --fasta "$pilot/inputs/contigs.fasta" --coverage_tsv "$pilot/inputs/coverage.tsv" \
    --protein_db_dir "$pilot/db/diamond/uniprot_reference_proteomes/ctlk-panel-01" \
    "${nt[@]}"
)
```

The pilot reuses mean depths from production `samtools coverage` defaults.
The FASTA and coverage table must have identical sequence names and lengths
before subsetting. If only the coverage table is unavailable, omit
`--use-trace-coverage` during preparation and `--coverage_tsv` during the
evidence run. The R plot will use GC versus length.

Set `CONTAMINATION_NT_DB` to the absolute prefix of an existing taxonomy-enabled
local BLAST database to enable whole-input MEGABLAST and targeted BLASTN. Without it the report
lists target windows, total target bp, and a FASTA for follow-up. A panel of
eight proteomes can establish that the search and reporting work on real
contigs, but gives limited taxonomic sensitivity and specificity. Use the
protein and nucleotide evidence tables to evaluate individual calls. Review
intervals describe aligned sequence, and the pilot review sheet is bound to
the subset FASTA. Carrying any decision to the full post-FCS FASTA requires a
new source-bound review.

Open `comparisons/contamination-ctlk-pilot-01/results/review/README.md`. For
review, return `inputs/selection.json`, `inputs/selected.tsv`,
`panel/proteomes.tsv`, the Markdown and TSV/JSON files in `results/review`,
and both PNGs in `results/diagnostics`.

## Run the whole post-FCS haplotype with NCBI nt on Crest

Supplying `--nt_db_root` and `--nt_release` enables managed setup. The default
collection is NCBI `nt`. The workflow downloads the exact catalog snapshot
advertised for that date using `update_blastdb.pl`, opens the resulting BLAST v5
database, and stores the complete bundle under
`<nt_db_root>/blast/nt/<nt_release>/database/`. Future runs using that label
reuse the stored database, even after NCBI publishes a newer snapshot. For a
fresh download the date must match the current catalog. The manifest records
the catalog entry, volume inventory, and BLAST version. An explicit
`--nucleotide_db /absolute/path/to/nt` instead uses an existing database.

The full `nt` download requires substantial shared space. Check the catalog
for the currently advertised date and sizes before launching:

```bash
curl -fsSL https://ftp.ncbi.nlm.nih.gov/blast/db/blastdb-metadata-1-1.json |
  python3 -c 'import json, sys; x=next(x for x in json.load(sys.stdin) if x["dbname"]=="nt"); print(x["last-updated"][:10], "compressed GB", round(x["bytes-total-compressed"]/1e9), "expanded GB", round(x["bytes-total"]/1e9))'
```

This command uses the full post-FCS FASTA and coverage paths recorded by the
earlier subset pilot. Set `NT_RELEASE` to the date shown above and set
`PROTEIN_DB_DIR` to the intended shared DIAMOND bundle if available. Otherwise
the run uses the eight-proteome pilot panel and establishes whole-input
execution against `nt` while full protein-database benchmarking remains separate.

```bash
module load nextflow/23.10.1 miniconda3
NT_RELEASE=YYYY-MM-DD
(
  set -euo pipefail
  pipeline="$PWD/gcl_genome_assembly"
  pilot="$PWD/comparisons/contamination-ctlk-pilot-01"
  run="$PWD/comparisons/contamination-ctlk-whole-01"
  selection="$pilot/inputs/selection.json"
  fasta=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source"])' "$selection")
  coverage=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["coverage_source"])' "$selection")
  protein="${PROTEIN_DB_DIR:-$pilot/db/diamond/uniprot_reference_proteomes/ctlk-panel-01}"
  mkdir -p "$run"
  nextflow -log "$run/nextflow.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_evidence.nf" -profile slurm \
    -w "$run/work" --outdir "$run/results" \
    --assembly Sde-CTlk_104_hap1 --stage post_fcs --taxid 373251 \
    --taxonomy_dir /work/birdlab/databases/ncbi_taxonomy \
    --fasta "$fasta" --coverage_tsv "$coverage" --protein_db_dir "$protein" \
    --nt_db_root /work/birdlab/databases --nt_collection nt --nt_release "$NT_RELEASE"
)
```

MEGABLAST searches every eligible window in parallel, including windows with
host-compatible protein hits. Eligible windows contain at least 100 A/C/G/T
bases and retain original sequence coordinates. The default run finishes the
searches after BLASTX and MEGABLAST. The Markdown report prioritizes concordant
and extended non-host evidence. Whole-contig candidates require distinct hits
to one non-host group across at least half the contig and at least 80% of the
distinct non-host span, with host hits covering at most 10% of the contig.
The contig table shows host-aligned and non-host spans as bases and percentages
for those candidates. Full contig summaries are available in the TSV. The interval
table compares representative non-host hits with overlapping host hits from
each search, with accessions, identity, query span and bit score. Protein and
nucleotide agreement is shown separately from agreement between the two
searches of the same nucleotide database. When BLASTN follow-up has been run,
the interval table includes its results for the nominated loci. Window classifications and all hit
and manual review rows are available as TSVs.

Read the planned follow-up query count and base total in
`results/review/README.md` and inspect the exact coordinates in
`results/review/nucleotide_targets.tsv`. To run this plan, resume with
`--targeted_blastn true`. It includes up to one locus per source window, with
500 bp of flanking sequence where available within that window. A locus
qualifies with at least 1 kb of distinct non-host hit coverage or at least
250 bp of overlapping distinct hits from both screens in the same taxonomic
group. The plan prioritizes concordance and limits the total query span with
`--contamination_followup_max_bp` (default 100,000 bp). Targets are packed into
roughly 25 kb search batches. All review actions still require human selection.

Reports preserve the three hit tables separately. The diagnostics include
sequence-span-weighted GC and depth margins and contig counts, span, and N50
by displayed taxonomic group.

### Rebuild a report from completed searches

Use `--previous_report /path/to/results/review` with `--fasta` pointing to the
exact original FASTA and `--taxonomy_dir` pointing to the local taxonomy.
Use a new `--outdir` to keep both reports available. The report checks the FASTA
checksum and uses the published hit tables, coverage, and BLASTN target plan
from that report. This mode runs only reporting and diagnostics.
Click a review ID in the interval shortlist to open its alignment-span plot.

### Search concurrency

| Parameter | Default | Concurrent tasks per Nextflow run |
| --- | ---: | --- |
| `--contamination_blastx_max_forks` | 8 | DIAMOND BLASTX chunks |
| `--contamination_blastn_max_forks` | 8 | MEGABLAST chunks and, separately, targeted BLASTN chunks |

Each search task requests 8 CPUs. BLASTX and MEGABLAST can overlap, so at the
defaults they can run up to 16 search tasks together. Both limits accept positive
integers. The standalone entry point takes one assembly. In the assembly
pipeline, all haplotypes share each process's limit. Separate Nextflow runs
have independent limits.

### Exclusive node workers for MEGABLAST

Set `--contamination_megablast_execution nodes` to search a cohort through
exclusive node allocations. The existing `chunks` mode remains available.
DIAMOND and optional targeted BLASTN retain their separate process limits.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--contamination_megablast_nodes` | 7 | Maximum simultaneous exclusive search workers |
| `--contamination_megablast_slots` | 8 | Simultaneous searches inside each worker |
| `--contamination_megablast_threads` | 8 | Threads per search |
| `--contamination_megablast_scratch` | `/scratch` | Node-local directory used for database and query copies |
| `--contamination_megablast_time` | `96h` | Wall-time limit for each worker and assembly collector |

Each worker requests slots multiplied by threads CPUs, exclusive access, and
all node memory via Slurm `--mem=0`. Configure worker count, scratch directory
and any required node constraints for the local cluster. The driver allocation
must outlast its workers and any retries. Seven 64-CPU workers leave room for
a driver and other tasks on an eight-node allocation limit.

Preparation produces a cohort queue that interleaves assemblies. Each worker
copies and checksum-verifies the complete nucleotide database once, checks
scratch capacity, and dynamically assigns chunks to its available search slots.
The private scratch directory is removed when the worker exits normally.
Workers with insufficient free scratch record the hostname, available capacity
and required space, then release their allocations. Other workers continue
processing the shared queue. A final completion check fails with a missing-chunk
summary if the workers finish without completing the queue. Exclusive allocation
does not remove old scratch files; retained copies from completed benchmarks
can reduce the free space available to a new worker.
The final worker-capacity and chunk-completion summary is published as
`pipeline/contamination_megablast_completion.json`.

Completed chunks are retained under `work/contamination_megablast/` when using
the default work directory. A completion record binds each output to the exact
query checksum, database file metadata and manifest checksum, BLAST version,
threads and search flags, and the hit-file checksum. The input database must be
an immutable snapshot. Changed identities create distinct entries. Corrupt or
incomplete results are recomputed. This is the checkpoint unit inside the larger
Nextflow node task; normal Nextflow caching still handles the surrounding tasks.

Workers coordinate using Linux file locks on shared storage. The Crest acceptance
fixture tests those locks across two physical nodes. An interruption releases
the claims while retaining already committed chunks. Completion records are
written atomically. Assembly collectors release results independently, allowing
reports to begin before the cohort finishes. Empty query sets still produce reports.

The old per-chunk Nextflow tasks have a different execution identity. Their
completed outputs require an explicitly validated one-time import before they
can serve as checkpoints for node workers. Preserve their work directories and
run logs during that transition.

The bounded native acceptance fixture is `tests/contamination_node_fixture.py`.
It creates a tiny nucleotide database and 24 planted queries across two synthetic
assemblies, tests shared filesystem locks with two Slurm ranks, verifies all
searches, then exercises recovery of one missing completion record. It does not
require the production nucleotide database.

## Run post-FCS review in the assembly pipeline

Add these options to the usual assembly command to generate review evidence
and stop before scaffolding:

```bash
--run_decon_contigs true \
--run_contamination_evidence true \
--contamination_protein_db_dir "$protein_bundle" \
--contamination_nucleotide_db "$nt_prefix" \
--stop_after_contig_review true
```

The protein bundle contains the DIAMOND index and its manifest. For a managed
NCBI nucleotide database, replace the local prefix option with
`--contamination_nt_db_root "$db_root" --contamination_nt_release "$nt_release"`.
Sensitive BLASTN is optional via `--contamination_targeted_blastn true`.
Evidence is published under `assembly/contamination/<assembly-id>/`.
By default, read support retains the indexed BAM in Nextflow work and uses it
for coverage summaries and the shortlisted locus plots. If post-FCS mapping QC
is enabled, it uses the same mapping task. Set
`--contamination_read_support false` to generate evidence without the mapping.

For the next run, preserve the original assembly and FCS inputs, add
`-resume --contamination_review "$reviewed_tsv"` and set
`--stop_after_contig_review false`. Use the same Nextflow launch directory and
work directory to resume. The reviewed sheet can contain multiple assembly IDs.
The workflow applies selected rows on the matching source FASTA and publishes
the removal and coordinate audit in each affected assembly's `cleanup/audit/`
directory. Assemblies with zero selected rows pass through using their original
FASTA. Evidence generation can be enabled separately from applying the sheet.

The cohort evidence and contig-stage stop can be checked separately from reviewed
cleanup. A final-assembly evidence pass is a separate follow-up.

### Reviewed-cleanup acceptance

The bounded runner uses the actual `CONTAMINATION_POST_FCS` workflow with synthetic
FASTA files and decision sheets. It does not run FCS, alignment searches, or
scaffolding, and it does not read or alter production assemblies. It checks:

- Whole-sequence removal, an internal split, and terminal trimming against exact
  expected DNA, preserving case and ambiguous bases.
- Both an explicitly unselected assembly and an assembly absent from the sheet:
  their original FASTA paths and bytes must reach the downstream consumer.
- Applied decisions, output checksums, retained/removed base accounting, and
  coordinate partitions that reconstruct every original sequence exactly.
- An unchanged repeat fully cached, a changed review taking effect while the two
  untouched downstream inputs remain cached, all decisions unselected, and the
  cleanup switch disabled.
- Expected failures for a wrong source checksum, wrong coordinate stage, unknown
  assembly, and overlapping selected intervals.

After syncing the pipeline, run from the project directory on Crest:

```bash
module load nextflow/23.10.1 miniconda3
python3 gcl_genome_assembly/scripts/comparisons/contamination_cleanup_acceptance.py run \
  --profile slurm --out "$PWD/comparisons/contamination-cleanup-acceptance-01"
```

The driver must remain running for all nine phases; it can run inside a small
Slurm driver allocation. Each fixture task requests one CPU and 1 GB. Some runs
are intentionally expected to fail; the runner checks the specific rejection
and continues. Overall success is `PASS` in `acceptance_summary.json`. The runner
also writes `<output-directory>.review.tar.gz` without its Nextflow work/cache
directories. Use a new output directory for each acceptance run.

For local testing with Python already available, use `--profile local --no-conda`.
The test entry point is [tests/contamination_cleanup.nf](../tests/contamination_cleanup.nf),
driven by [the acceptance runner](../scripts/comparisons/contamination_cleanup_acceptance.py).

Local foundation checks:

```bash
python3 -m unittest discover -s tests -p test_contamination_review.py -v
nextflow run tests/contamination_grouping.nf
```

These use small, known DNA sequences; no species-specific rules or CTlk coordinates.
The Nextflow grouping fixture uses 2,088 synthetic hit paths across ten haplotypes
and an empty assembly. It checks mixed interpolated/plain-string assembly IDs,
both evidence sources, coverage joins, and unchanged task hashes after ID
normalization. No database downloads or alignment tools are required.
