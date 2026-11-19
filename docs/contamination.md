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

The next Crest checks cover cohort evidence, the contig-stage stop, and reviewed
cleanup with normal Nextflow cache invalidation. A final-assembly evidence pass
is a separate follow-up.

Local foundation checks:

```bash
python3 -m unittest discover -s tests -p test_contamination_review.py -v
```

These use small, known DNA sequences; no species-specific rules or CTlk coordinates.
