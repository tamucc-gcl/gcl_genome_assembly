# Contamination evidence and reviewed cleaning

## Implementation status

The standalone preparation/cleaning foundation is implemented. It is **not yet
called by `main.nf`**. Database preparation, taxonomic searches, coverage integration,
classification, reports and the contig-only stopping control are subsequent work.
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

## Agreed next implementation steps

1. Prepare taxonomy-aware reference-proteome DIAMOND and nucleotide BLAST bundles
   with native `storeDir`; preserve existing FCS. Follow the shared database layout
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
2. Search windows, select/record targeted nucleotide follow-up, and incorporate
   source-compatible HiFi or shotgun coverage. Expected-host rank defaults to
   phylum, configurable. Non-host assignments are review evidence, not deletion.
3. Produce a Markdown report, R/ggplot2 figures and BlobToolKit dataset. Test known
   mixtures and real contigs on Crest before main-workflow integration.
4. Add the pre-scaffolding workflow and general contig-stage stopping control.
   Test zero downstream task submissions and subsequent normal cache reuse.
5. Validate reviewed cleaning/downstream invalidation and the optional final check.

Local foundation checks:

```bash
python3 -m unittest discover -s tests -p test_contamination_review.py -v
```

These use small, known DNA sequences; no species-specific rules or CTlk coordinates.
