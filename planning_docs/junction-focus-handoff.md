# Focused junction follow-up

This standalone diagnostic reuses job 1506302 and does not change production
defaults, map HiFi reads again, rerun assemblies or emit cuts. Investigation
coordinates live in `planning_docs/investigations/ctlk-junction-focus.tsv`;
the algorithms have no CTlk-specific branches.

## Questions addressed

- J03/J04: do the inserted component's ends have support beyond background,
  and do displaced anchors support a different order or adjacency?
- H01: does its depleted immediate flank become informative farther from the gap?
- H02: assess the two weak joins at 94.575 and 95.222 Mb individually.
- J05: where do the same segments of low-confidence HiFi molecules also align?
  Separate other-haplotype placements from alternatives within the same assembly.

The focus table includes the five exact gaps and the J05 ambiguity interval.
All existing comparison gaps on the focused scaffold are retained. H03/H04 and
the other continuously supported raw-contig transitions are not rerun here.

## Measurements

1. Stream the retained hap1 Hi-C pairs once, accumulating 25/100/250/500 kb flank
   contacts by library and alternative partner bins. Record flank length,
   ambiguous bases, intervening gaps and external contact opportunity.
2. Require the original 100 kb counts to reproduce exactly, including margins and
   total projected UU counts. Check source FASTA/AGP hashes, dictionaries and
   pairs row count. A mismatch fails the run rather than producing mixed results.
3. Extract 10 kb sequence anchors at offsets 0/25/100/250 kb from each gap and
   align them against the same twelve retained current/historical assemblies.
   Retain all emitted placements and summarize pairs with >=80% aligned query
   bases. Show target, orientation and signed separation; show missing flanks and
   truncated combinations explicitly. Comparisons enumerate up to five placements
   per side; raw placements retain the rest. Thresholds are descriptive screens.
4. Decode retained HiFi CIGARs for J05 into original molecule coordinates,
   including strand and hard clips. Report alternatives aligning the same read
   segment, with MAPQ, whole-alignment AS/NM, and query overlap. Select primary
   reads with >=1 kb aligned inside the focal region. The 80% alternative-overlap
   summary is not a calibrated uniqueness threshold. No remapping is performed.

Displaced anchors can cross another join, especially around J03/J04. The anchor
table explicitly records this; such a match cannot alone validate the intervening
joins. Hi-C scales overlap and are not independent replicates. These results do
not correct fully for physical gap length, local distance decay or repeat
mappability, and do not yield a statistical misjoin probability. Self-alignments
are checks, historical assemblies share inputs, and the Hi-C libraries were used
for scaffolding. Alternative contacts identify leads, not authorized new joins.

Graph-overlap reconstruction and any final manual cut remain separate reviewed
steps if these measurements do not resolve the relevant connection.

## Cluster commands

After syncing, from the project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
bash -n gcl_genome_assembly/scripts/comparisons/run_junction_focus.sbatch
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

Six new synthetic helper tests cover clipped/reverse read coordinates, deletions,
overlap union, complete anchors, order/orientation and invalid CIGARs. Tests have
not been executed locally; the cluster run is also required integration validation.

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_junction_focus.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  --packet comparisons/junction-assessment-20261005-084623-1506302 \
  --focus gcl_genome_assembly/planning_docs/investigations/ctlk-junction-focus.tsv
```

Keep the original unarchived packet and all FASTAs/pairs referenced by its
provenance. The job requests eight CPUs, 96 GB and an exclusive node; the 24-hour
limit is a limit, not a runtime estimate. The main costs are streaming hap1 pairs
and constructing the twelve peer indexes. Those operations run sequentially.

Return `comparisons/junction-focus-TIMESTAMP-JOBID.tar.gz` and
`logs/junction-focus-JOBID.out`, as printed by the job. Inspect status.json first,
then contact_scales.tsv, anchors/windows.tsv, anchors/flank_correspondence.tsv,
and the J05 read_summary/read_alternatives tables. No assembly rerun is needed.
