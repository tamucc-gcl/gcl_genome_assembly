# Contamination evidence and reviewed cleaning

## Implementation status

The standalone tools cover reviewed cleaning, shared reference database setup,
whole-assembly MEGABLAST, translated protein searches, and targeted BLASTN.
They are **not yet called by `main.nf`**. The search accepts an assembly-matched
coverage BAM or a retained `samtools coverage` table. Main-pipeline coverage
generation and the contig-only stopping control are subsequent work.
The existing FCS workflow is unchanged. Do not enable the currently blocked
`run_blobtools_evidence` flag yet.

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
bases and retain original sequence coordinates. The optional targeted
`blastn` pass (disable with `--targeted_blastn false`) searches windows with
non-host, mixed, or unresolved taxonomic hits. A window with no hits from
either screen stays visible in the no-hit counts. Reports preserve all three
hit tables separately and list candidate bp separately from contig span.
The result includes sequence-span-weighted GC and depth margins, plus
contig counts, span, and N50 by displayed taxonomic group.

## Remaining implementation steps

1. Verify the whole-haplotype run on Crest with the intended databases. Evaluate
   window-boundary sensitivity and the size of the nucleotide target set.
2. Add pipeline-generated HiFi or shotgun coverage, the pre-scaffolding workflow,
   and the general contig-stage stopping control.
   Test zero downstream task submissions and subsequent normal cache reuse.
3. Validate reviewed cleaning/downstream invalidation and the optional final check.

Local foundation checks:

```bash
python3 -m unittest discover -s tests -p test_contamination_review.py -v
```

These use small, known DNA sequences; no species-specific rules or CTlk coordinates.
