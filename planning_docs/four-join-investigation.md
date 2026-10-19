# Focused investigation of four retained joins

This is a removable diagnostic diversion. It does not change the production
chimera policy, apply cuts, rebuild assemblies, or change pangenome construction.
The only sample/event selection is the small
`scripts/comparisons/ctlk_four_join_focus.tsv` configuration. The measurement code
does not contain sample names, chromosome pairs, or known breakpoint coordinates.

## Cases and comparison regions

| Assembly | Case | Connection | Role |
|---|---|---|---|
| CTlk hap1 | C02 | chr7 → chr12 | Target |
| CTlk hap1 | C03 | chr4 → chr14 | Target |
| CTlk hap2 | C01 | chr12 → chr7 | Target, including apparent retention support |
| CTlk hap2 | C02 | chr14 → chr4 | Target |
| CTlk hap1 | C01 | chr9 → chr7 | Previously reviewed break benchmark |
| CTlk hap1 | C04 | chr2 → chr5 | Previously reviewed break benchmark |

Benchmarks are comparisons, not assumed ground truth. Both haplotypes of an
individual count once for independent chromosome context.

## Execution

The launcher checks the published v7 evidence against its actual retained
Nextflow task inputs. It verifies FASTA, source FASTA, AGP and assembly-alignment
checksums, BAM identity and sequence dictionaries. Missing cache inputs stop the
diagnostic with a named input; nothing silently falls back to a different run.

Two chromosome-evidence jobs first project the retained CIGAR alignments into
bounded local windows. CIGAR blocks are clipped before the interval sweep so
whole-chromosome interval processing is not repeated. The jobs compare
chromosome-level assignment with unique-locus assignment, preserve competing
alignments, and quantify opposite-side chromosome islands.

Their unique-arm masks require concordant unique placements in at least two
independent eligible individuals. Conflicting evidence remains unassigned.
This is conservative emitted-alignment evidence, not a complete repeat annotation.

After each haplotype's labels finish, its HiFi and Hi-C jobs run concurrently.
Both haplotypes can proceed at once: up to four heavy evidence jobs on separate
Slurm allocations. Resource requests reflect the implementation:

| Job | Count | CPUs per job | Memory | Time limit |
|---|---:|---:|---:|---:|
| Chromosome projection | 2 | 1 | 32 GB | 12 h |
| Targeted competitive HiFi mapping | 2 | 8 | 64 GB | 24 h |
| Hi-C scan | 2 | 1 | 32 GB | 24 h |
| Combined report | 1 | 1 | 12 GB | 6 h |

The projection and pairs scanner are serial programs and request one CPU. HiFi
mapping uses its eight requested threads. The controller checks provenance and
coordinates, orchestrates the jobs and packages the review.

### HiFi measurements

- Reuse the indexed original BAM to extract primary molecules in all selected
  windows, deduplicated across windows.
- Measure original spanning support every 1 kb and at every retained R/gap site.
- Map those selected reads against each complete haplotype separately, retaining
  secondary alignments. This preserves within-haplotype placement competition;
  similarity to the other haplotype is reported separately.
- Report molecule identity, clipping, large indels, competing placement scores,
  support with 1 kb and 5 kb anchors, and unique chromosome-specific bases reached
  on each side within the same uninterrupted molecule chain.
- Estimate the minimum span needed to reach unique anchors and compare it with
  local read lengths. This distinguishes a missing bridge from an interval reads
  cannot realistically bridge. It is not a model of expected bridge counts.

Read extraction is not a new exhaustive search for reads having only secondary
placements in the focal region. Emitted mapping alternatives remain a limitation.
The original graph-path context is preserved; no new hifiasm run is performed.

### Hi-C measurements

- Scan the retained pairs once per haplotype for all selected regions.
- Preserve library identity; a single library is sufficient. Pooled maps are
  supplemental, not independent evidence.
- Lift pairs through the exact source AGP. Reproduce existing 100 kb / 300 kb
  flank measurements exactly before accepting the finer-resolution results.
- Build exact 10, 25, 50 and 100 kb maps from 5 kb internal counts. When coarsening,
  within-bin pairs are counted once rather than twice.
- Scan positions every 10 kb with fixed 100, 300 and 500 kb flanks. Report raw
  cross counts, within-flank counts and undefined ratios explicitly.
- Measure contacts between uniquely assigned arms outside the ambiguous region,
  at 50, 100 and 300 kb of selected sequence per arm. Retain the exact selected
  intervals, including discontinuities and distances.
- Compare those arm contacts with translated, same-chromosome unique-bin
  geometries at matched separation and broadly comparable within-flank counts.
  An expected count is shown only with at least five qualifying control windows.
  Overlapping controls are descriptive, not independent statistical replicates.
- Separately retain gap-free same-chromosome scan controls. Grey observability
  marks do not erase measured zeros; the actual flank counts remain visible.

## Run on Crest after syncing

From the project directory containing `gcl_genome_assembly`, `work`, and
`comparisons`:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_four_join_investigation.sbatch \
  ctlk_analysis \
  comparisons/chimera-observable-contact-v7/results/assembly/chimeras \
  comparisons/core-graph-assignment-01/results/assembly/final \
  comparisons/four-join-investigation-01
```

Arguments are analysis environment, original evidence packet, current final
name-map directory, and diagnostic output directory. The current names are
context only: trimming and extension can shift positions, so the report never
claims an unverified coordinate lift. All measurements and eventual manual cuts
remain tied to the pre-finishing FASTA and its checksum.

Running the same command again resumes the diagnostic through normal Nextflow
caching. Unchanged small input manifests retain their timestamps. No production
assembly or pangenome process is launched by this script.

Main output:

```text
comparisons/four-join-investigation-01/review/README.md
comparisons/four-join-investigation-01/review/sites.tsv
comparisons/four-join-investigation-01/review/manual-cut-worksheet.tsv
comparisons/four-join-investigation-01.review.tar.gz
```

The archive includes the report, plots, compressed local count matrices,
measurements, provenance, trace and failure status. It excludes the work tree,
assembly FASTAs, BAMs, SAMs and extracted read sequences. Published review output
uses symlinks; the portable archive dereferences them.

## Decision after the batch

For each target, assess separately:

1. Whether the chromosome connection has credible continuity.
2. Whether chromosome-labelled islands are uniquely placed or multiply placed.
3. Whether contact absence is unexpected given usable sequence and coverage.
4. Whether there is a defensible exact cut, a range, or no reason to cut.

The worksheet starts with every selection `NO` and every cut coordinate blank.
An R probe or transition midpoint is not automatically a recommended cut.
After review, retain, cut, or identify the precise unresolved evidence for each
case. Incorporate only methods that materially improve those decisions into the
general workflow, then return to the core pangenome assignment issue.

## Validation

`tests/test_join_investigation.py` checks general interval, uniqueness, molecule,
single-library, reverse-orientation, count-conservation and resume behavior.
`tests/setup_join_investigation.py` generates synthetic real FASTA/BAM/PAF/pairs
inputs for an end-to-end Nextflow 23.10.1 run with real samtools and minimap2.
These test software behavior; they do not establish the biology of CTlk joins.

Local validation on 2026-10-19 passed 18 general tests and the real-tool
Nextflow 23.10.1 smoke run. The supported synthetic junction retained 12
unique-arm bridges; the unsupported junction retained zero, with an observable
zero Hi-C cross count and 200 within-arm pairs on each flank at the 100 kb
anchor scale. Exact legacy-resolution replay, all four map resolutions, report
links, unselected cut worksheets and a seven-task cached resume were checked.
The smoke report and logs are in the ignored
`comparisons/four-join-software-validation` directory. CTlk evidence collection
still requires the Crest launch above.
