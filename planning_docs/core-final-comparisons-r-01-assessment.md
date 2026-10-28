# R comparison pilot acceptance — 2026-10-28

The CMat native pilot passes execution, final-input reconciliation and visual
inspection. Expected chromosome-scale synteny has recovered. The remaining
comparison gate is the complete five-pair cohort; no new code fix is indicated
by this packet.

## Execution and provenance

Exit 0; five completed tasks, no failures, retries or ignored tasks. Metadata
selection produced exactly CMat's two finalized haplotypes and one pair.
Minimap2's recorded command uses asm5 without `-c` or `--eqx` and one index.
Native plotting versions are R 4.5.3, ggplot2 4.0.1, pafr 0.0.2 and gggenomes 1.2.0.

Both packaged FAI fingerprints, base totals and source FASTA paths reconcile
independently with the manifest. The FAI contents are byte-identical to the
earlier sparse comparison run. Whole FASTAs are not packaged and have not been
independently rehashed in this review.

Mapping took 4m45s and 33 GB peak RSS, versus 23m32s and 42.1 GB in the earlier
base-level CMat comparison. R plotting took 5.5 seconds. First-launch resume and
unmatched process-selector warnings are harmless in this bounded entry.

## Recovered mapping evidence

| Filtering | Records | Reference chain-span coverage | Query chain-span coverage |
|---|---:|---:|---:|
| Unfiltered | 6,028 | 96.80% | 99.72% |
| MAPQ >=5 | 4,266 | 96.74% | 99.39% |
| MAPQ >=5, length >=10 kb | 1,040 | 96.63% | 98.99% |
| MAPQ >=5, length >=50 kb | 427 | 96.32% | 97.93% |
| MAPQ >=5, length >=150 kb | 185 | 95.88% | 96.44% |

Fractions use interval unions divided by the full assembly lengths. The table's
fractions, monotonic filtering and agreement with the QC summary were independently
checked. This is mapped-chain span coverage, not base-by-base alignment through
every interval. Overlapping mappings are counted once.

The final chain block-length sum is 1,156,291,347 bp, versus 126,430,477 bp in the
earlier base-level comparison: about 9.15 times larger. These block sums are not
union coverage and the old packet cannot supply an exact old coverage fraction.
There are fewer retained records now (185 versus 321), because the recovered
chains span much longer regions. Record count alone was not a completeness metric.

Both R figures draw 180 retained mappings on 15 chromosome-size sequences per
haplotype. The five other retained mappings are excluded by the 1 Mb sequence
display threshold, rather than being lost from the evidence summary. The dotplot
shows the expected near-continuous chromosome correspondence; the riparian plot
again shows broad chromosome ribbons. Both have white backgrounds and readable
chromosome identities. Minor subtitle/label crowding is presentation work, not
a new measurement blocker.

The approximate match fraction of 0.574 is calculated from minimap2's approximate
mapping-chain counts. It must not be interpreted as 57.4% nucleotide identity,
nor compared directly to the earlier base-level 0.893 value. The report explicitly
labels its basis. No identity or biological completeness conclusion is drawn.

## Next acceptance run

Run the same production comparison entry without the diagnostic sample selector:

```bash
mkdir -p logs comparisons
sbatch gcl_genome_assembly/scripts/comparisons/run_final_comparison_checks.sbatch \
  comparisons/core-operational-01/results/assembly_eligibility.json \
  comparisons/core-final-comparisons-r-all-01
```

Return `comparisons/core-final-comparisons-r-all-01.review.tar.gz`. It runs the
five pairs concurrently, without rebuilding assemblies or the pangenome. Review
cohort coverage and figures before closing the optional comparison subgate.
Category 7 remains open for its other operational/dependency/integration work;
R figure migration and unified reporting remain category 6.
