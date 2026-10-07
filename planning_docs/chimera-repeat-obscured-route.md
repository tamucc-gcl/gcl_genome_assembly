# Repeat-obscured gap adjudication

The forward automatic action policy is gap-v2. There is one action schema; old automatic policy stamps are rejected. Manual source-bound actions retain their format.

Two routes can authorize literal-gap unjoins. The original direct route requires informative immediate HiFi flanks, zero qualified bridges, independent chromosome discordance, calibrated support loss in both libraries, checked alternatives, and noncontradictory graph context. The new repeat-obscured route does not require immediate flanks to be informative. It requires all of:

- A verified literal N gap with at least 5 Mb of assembly sequence on either side.
- Unique 50 kb peer anchors sampled 100, 250, and 500 kb beyond the gap. At least two distances must assign each side consistently to different qualified chromosome identities in each supporting peer.
- At least three independent qualified individuals agreeing on the same chromosome pair. Usable haplotype conflicts and usable same-chromosome peer context prevent agreement. Missing peer evidence is not agreement. Exactly one literal gap must lie between the nearest qualified block anchors; ambiguous gaps are unresolved.
- Informative farther HiFi flank coverage, and matched Hi-C controls at the same offsets, independently for every library. Each informative distance needs at least five continuous-sequence and five peer-continuous gap controls, within-flank counts within a factor of two, and at least 100 within-flank pairs per side. Two distances must show support loss in both libraries. Any informative trial lacking the required support loss vetoes this route.
- Zero qualified immediate seam bridges, no qualified positive native continuity contradiction, and checked alternative placements. Unknown graph context prevents authorization. Same primary-contig membership alone is not a molecular continuity veto.
- No contradictory usable block assignment in the individual's other haplotype. Missing counterpart block evidence remains explicitly missing in the packet, rather than claimed homozygosity.

A nearby gap may now become automatically eligible when these independent measurements uniquely localize it; proximity to the original transition alone remains insufficient. Adjacent accepted cuts must also satisfy the selected-cut-set minimum-piece guard. Cuts preserve all sequence and gap Ns and rerun cohort chromosome assignment.

Block evidence samples chromosome context out to 550 kb from each side. This is not whole-chromosome basewise alignment or proof of a fusion. The policy uses conservative heuristic thresholds, not calibrated probabilities. The three-count allowance and tenfold contact-loss threshold are retained. Do not treat the synthetic regression tests as biological sensitivity/specificity estimates.

Implementation expands diagnostic mapping windows to 850 kb either side, while keeping the depth tracks limited to 100 kb either side. Both immediate and farther measurements reuse the same source-bound HiFi BAM. All Hi-C windows are tallied in one pass per assembly, with no pooling of libraries. Contact geometry is matched by offset; literal gap controls were already selected by gap length. Human packets include chromosome_blocks.tsv, per-distance HiFi measurements, raw contact counts, each matched control ID, and farther_contacts plots in addition to the existing plots and IGV sessions.

Local validation includes a deliberately constructed cross-chromosome query, ambiguity from a second intervening gap, continuous peer disagreement, positive HiFi/native vetoes, per-library disagreement at one farther distance, and automatic source-bound action emission from a formerly review-only hypothesis. 98 chimera tests and 10 pre-finishing tests pass. Python 3.11 syntax and whitespace checks pass. Actual tools, plots, thresholds, and biological outcomes require the Crest run.

After syncing, run from /work/birdlab/GCL/spratelloides_delicatulus_genome:

```bash
sbatch --array=0 gcl_genome_assembly/scripts/comparisons/run_chimera_real_test.sbatch ctlk_analysis
```

No new assemblies or read mappings are required. Wider peer context must be mapped because previous packets contain only narrow windows. Return auto.review.tar.gz. This is an isolated workflow test and may make automatic cuts only in its own output. H01 is not hard-coded, and qualification of H01 is not predetermined.
