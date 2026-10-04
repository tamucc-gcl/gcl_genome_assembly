# Historical comparison before junction decisions

Order agreed: (1) compare current assessment assemblies with both historical
final haplotypes, (2) investigate the two YaHS joins around the 113.9 kb component,
(3) refine five internal raw-contig transitions, (4) record per-junction decisions.

This is a standalone investigation, not a new production pipeline branch.
Historical haplotype labels and chromosome names are not assumed to correspond.
All four haplotype combinations are aligned sequentially. Historical final
assemblies may include finishing edits; disagreement is not automatically an error.
The supplied old name maps are interpretation aids, not coordinate liftovers.

The comparison retains whole-assembly PAFs with secondary alignments, per-pair
union coverage tables, base-resolved candidate-window PAFs and exact outer-anchor
placements. Current assessment FASTA checksums must match the candidate tables.
No breakpoint is projected from chromosome names or approximate PAF endpoints.
Whole-assembly mapping is assembly alignment, not an independent read check.
Single reported matches do not certify uniqueness. Haplotype correspondence can
change by region; a best whole-genome match must not force local correspondence.

Run from the cluster project root after syncing the new scripts:

```bash
export PYTHONDONTWRITEBYTECODE=1
bash -n gcl_genome_assembly/scripts/comparisons/run_chimera_history.sbatch
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

If checks pass:

```bash
old="$PWD/tst/genome_assembly_store/assembly/final"
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_history.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  --results comparisons/chimera-foundation-01 \
  --assembly Sde-CTlk_104_hap1 --assembly Sde-CTlk_104_hap2 \
  --old "old_hap1=$old/Sde-CTlk_104_hap1.fasta" \
  --old "old_hap2=$old/Sde-CTlk_104_hap2.fasta"
```

Return the printed comparisons/chimera-history-TIMESTAMP-JOBID.tar.gz and
logs/chimera-history-JOBID.out. FASTA sequences are excluded from the archive.
Keep outputs and caches in place; no Nextflow rerun is needed.

Next decision gates:
- Historical comparison: establish region-level correspondence, contiguity and
  alternative placements; do not treat a cleaner historical naming scheme as truth.
- YaHS component: confirm its chromosome assignment and whether the historical
  sequence is contiguous with either current neighbour; then assess library-specific
  Hi-C alternatives and local HiFi mapping with actual read opportunity.
- Raw contigs: refine assignment switches using unique sequence, test read support
  across narrow junctions competitively against both haplotypes, and inspect retained
  graph/read placements. Long unspanned uncertainty intervals are uninformative.
- Decisions: separate origin, adjacency evidence, ambiguity and cut eligibility.
  Retain unresolved events; do not automatically cut or declare biological fusions.
