# Join assessment batch

This is a standalone diagnostic, not a production chimera classifier. It makes
no assembly edits, changes no pipeline defaults and does not invalidate Nextflow
caches. The case-specific registry is data in
`planning_docs/investigations/ctlk-junction-registry.tsv`. The implementation has
no sample-specific coordinate or chromosome rules.

## What this pass implements

1. Read all eleven tracking entries with their exact assessment checksums.
2. Recover physical AGP gaps on the four requested scaffolds from assessment and
   earlier retained stages. Earlier-stage transfer requires identical complete
   scaffold sequence (either orientation), unambiguous correspondence, valid AGP
   layout and literal N-gap verification. Changed earlier scaffolds are explicitly
   unresolved; the script never uses a nearest-gap guess as a cut.
3. Link every tracking interval to all intersecting verified gaps. Select up to
   three nearby comparison gaps per interval, preferring the same placeholder gap
   length. Controls are not presumed correct. Record all recovered gaps separately.
4. Verify retained pre-finishing HiFi BAM provenance and the complete sequence
   dictionary. Extract reads mapped within 100 kb of all requested regions and
   comparison joins, then remap those selected reads competitively against the
   supplied current haplotypes. Retain alternative placements and count distinct
   primary molecules spanning both 1 kb flanks, with a separate CIGAR/MAPQ/NM screen.
5. Generate read-support profiles every 1 kb throughout the broad raw-contig
   intervals and flanks. A profile minimum is diagnostic, not an automatic cut.
6. Stream existing deduplicated UU Hi-C pairs once per assessment assembly and
   lift them through its last-round AGP. Validate source chromosome dictionaries.
   Resolve read-set prefixes to library IDs; separate libraries and retain an
   explicit UNASSIGNED category. Report crossing counts, library-depth scaling,
   external flank contacts and top alternative partner bins, including partners
   elsewhere on the same scaffold. Candidate and control joins use the same windows.
7. Align the same candidate/control windows against both historical haplotypes and
   all supplied current peers, keeping base-level CIGARs and secondary placements.
8. Collect native GFA read-placement/topology records for implicated raw contigs.
   This is not unitig-path reconstruction; purged component name aliases are not guessed.
9. Produce a review decision worksheet keyed to the stable tracker IDs. A supplied
   reviewed table can be validated with `--decisions`, including checksums and
   verified gap location for BREAK entries. No executable cut file is exported.

## What remains a scientific review, not an implemented automatic rule

- Determine whether historical target switches are true misjoins or correct joins
  of previously fragmented chromosomes.
- Assess phase consistency from competing placements and peer correspondence;
  mapping competition alone does not certify phase.
- Calibrate absence of HiFi bridges. Selected reads exclude originally unmapped
  reads, and physical distances across scaffold gaps are unknown. Zero support
  and low MAPQ cannot automatically justify a break.
- Judge narrow raw-contig breakpoints from sequence alignments, support profiles
  and graph records. Broad assignment intervals cannot be cut at their midpoint.
- Evaluate Hi-C after considering distance, mappability and coverage. Per-million
  counts only adjust library size, not these other effects. The scaffolding reads
  are not independent validation. Top alternative bins are not proposed joins.
- Transfer gaps across sequence-edited earlier stages using a further local
  alignment step if full-scaffold identity fails. The existing origin tracer can
  perform local exact-window tracing; this batch does not silently relax it.

## Cluster checks

Sync the new code and the investigation registry. From the project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
bash -n gcl_genome_assembly/scripts/comparisons/run_junction_assessment.sbatch
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

The synthetic tests cover interval identity, half-open bounds, reverse projection,
literal scaffold verification, no nearest-gap snapping, comparison selection,
ambiguous AGP lifts, library identities, pair ordering and decision validation.
They have not been executed locally. The cluster tests and first diagnostic run
are required validation before interpreting measurements.

## Run the complete evidence batch

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_junction_assessment.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  --registry gcl_genome_assembly/planning_docs/investigations/ctlk-junction-registry.tsv \
  --origins comparisons/chimera-origins-20261004-064841-1506290 \
  --results comparisons/chimera-foundation-01 \
  --hifi-results genome_assembly \
  --history comparisons/chimera-history-20261004-083159-1506291 \
  --peer-fasta-dir comparisons/chimera-foundation-01/assembly/scaffold/yahs_round2
```

The origin and history directories must be the original unarchived cluster
directories. Their provenance points to retained FASTAs. If the full HiFi context
was moved, change `--hifi-results`. Missing BAMs/pairs are reported as unavailable,
not zero evidence; mismatched hashes/dictionaries fail. `--prepare-only` runs gap
location and graph inventory without read mapping, contact counting or peer
alignments, if a lightweight first inspection is needed.

The automatic Hi-C layout adapter expects contig-decontamination input for
round 1 and scaffold-correction input for round 2, as in this investigation.
For other layouts (including scaffold decontamination), pass `--hic-inputs`
with a TSV containing `assembly`, `fasta`, and `sha256` for the exact FASTA used
for the last scaffolding round's Hi-C mapping. Paths resolve from the submission
directory. Sequence dictionaries are checked before contact counting; matching
names and lengths alone cannot establish that two FASTAs have identical bases.

This uses eight CPUs on an exclusive node with 96 GB requested. The 24-hour limit
is not a runtime estimate. Targeted read remapping, streaming full pairs and
building peer alignment indexes are the substantial costs. Processes run
sequentially. Existing whole-genome HiFi mappings and assemblies are reused.

## Return and review

The job automatically prints these two files:

```text
comparisons/junction-assessment-TIMESTAMP-JOBID.tar.gz
logs/junction-assessment-JOBID.out
```

The archive excludes FASTAs, FASTQs and SAMs. It includes compact alternative
read placements, support profiles, per-library contact tables, peer PAFs, graph
records, stage projection audits and the decision worksheet. Keep the unarchived
directory: it contains the sequence and alignment material for follow-up.

Review `status.json` first for which evidence was actually available. Then use
`seed_gap_links.tsv` to link J/H IDs to physical gaps. Inspect candidate/control
support and peer/graph context before filling `decisions.tsv`. Record both
supporting and conflicting evidence and a reviewer for KEEP/BREAK decisions.
For a raw-contig cut, a separately reviewed narrow coordinate and the production
manual-break contract are still needed; this batch only validates gap BREAK
entries. Validation never edits an assembly or asserts biological correctness.
