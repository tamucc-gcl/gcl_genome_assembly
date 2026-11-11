# Contamination evidence and reviewed cleaning

## Implementation status

The standalone tools cover reviewed cleaning, building a reference-proteome
database, screening a FASTA, and reporting protein and targeted nucleotide
evidence. They are **not yet called by `main.nf`**. The search pilot accepts a
prepared local BLAST nucleotide database and optional assembly-matched coverage
BAM. Shared nucleotide database setup, pipeline-generated coverage, and the
contig-only stopping control are subsequent work.
The existing FCS workflow is unchanged. Do not enable the currently blocked
`run_blobtools_evidence` flag yet.

The supported script for this foundation is `py_scripts/contamination_review.py`.
It requires Python 3.11 or later and has no third-party Python dependencies.
It processes one sequence at a time rather than loading a cohort into memory.

## Prepare an editable review sheet

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
GC/coverage and cumulative sequence plots. The BlobToolKit dataset in
`results/blobdir` includes protein-hit taxonomy and the combined review category.
All-ambiguous sequences remain in the evidence tables and are listed in
`blobdir_exclusions.tsv` because BlobToolKit requires measurable GC.
Coverage uses primary reads at MAPQ ≥20, omitting duplicates and QC failures.
Without coverage the static plot shows GC vs sequence length.

For a real pilot, supply a bounded set of actual post-FCS contigs, the sample
taxid, and a broader taxonomy-aware protein database. Use a BAM mapped to the
exact FASTA if available. The pilot checks that BAM sequence names and lengths
match and does not establish sequence identity from a BAM header alone.
`--nucleotide_db` is the prefix of a local, taxonomy-enabled BLAST database.
Omit it to see the selected targets without running nucleotide searches.
Windows with no protein hit, unresolved taxonomy, non-host or mixed evidence
are targeted. On a whole genome this can be a large fraction of the assembly,
so use a bounded input for the pilot. Protein and nucleotide hit scores are
kept separate. The report uses a configurable host rank (default phylum), and
review intervals reflect alignment evidence rather than proven breakpoints.

## Remaining implementation steps

1. Prepare the shared nucleotide BLAST bundle with native `storeDir`.
   Follow the shared database layout
   (`diamond/<collection>/<release>`, `blast/<collection>/<release>`, existing
   `ncbi_taxonomy`). Read-only inspection of the annotation pipeline confirms its
   `modules/diamond.nf` downloads the current complete Swiss-Prot FASTA and builds
   `swissprot.dmnd` with plain `diamond makedb --in ... --db ...`, stored via
   `storeDir`. That builder does not supply taxonomy mappings/nodes/names or retain
   a database release/checksum manifest; its DIAMOND dependency is not pinned in
   the download/search labels. Annotation uses the index for protein descriptions
   and repeat-library filtering. This identifies the source collection/build recipe,
   not the precise release or binary version of the existing Crest file.
   Leave this annotation database untouched. The agreed reference-proteome
   contamination collection is a separate taxonomy-aware bundle, not an alias for
   complete Swiss-Prot. The user confirms `proteins.dmnd` is unused; exclude it
   from the reuse plan without deleting or modifying it.
2. Run the known mixture and a bounded real-contig pilot on Crest. Evaluate
   window-boundary sensitivity and the size of the nucleotide target set.
3. Add pipeline-generated HiFi or shotgun coverage, the pre-scaffolding workflow,
   and the general contig-stage stopping control.
   Test zero downstream task submissions and subsequent normal cache reuse.
4. Validate reviewed cleaning/downstream invalidation and the optional final check.

Local foundation checks:

```bash
python3 -m unittest discover -s tests -p test_contamination_review.py -v
```

These use small, known DNA sequences; no species-specific rules or CTlk coordinates.
