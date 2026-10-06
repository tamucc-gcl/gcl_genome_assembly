# CTlk round-1 prevention experiments

Submit from the Crest project root after syncing the repository:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_yahs_round1_prevention_batch.sbatch ctlk_yahs comparisons/chimera-foundation-01
```

Six parallel tasks, no array throttle: baseline/cap10Mb/cap10Mb_min100kb for both haplotypes. Both EC stages remain enabled, as they already are in the current round-1 configuration. All runs reuse the exact decontaminated, unscaffolded round-1 contig FASTA and its matching contig-stage filtered Hi-C BAM. BAM reference names and lengths are checked before YaHS starts. They do not reuse the scaffold-stage BAM and do not use raw unfiltered hifiasm contigs.

The default enzyme list is GATC,CTNAG,TTAA,GANTC, matching the prior experiment. Pass a third argument of `none` if the retained baseline round-1 command omitted enzyme correction. The freshly rerun baseline is the controlled comparison for the two treatments; comparing solely to older outputs could conflate tool-version or command differences. YaHS version is recorded.

Baseline: minimum contig 20 kb, resolutions 20 kb through 200 Mb. Treatment 1 changes only the maximum resolution to 10 Mb. Treatment 2 also raises the minimum participating contig length to 100 kb. This excludes smaller contigs from scaffolding participation; it must not be scored as a repair simply because more sequence remains unplaced.

Evaluate intermediate/final AGPs by exact input component identity and interval, ignoring output rank and whole-scaffold reversal. Track all suspect baseline gap adjacencies, rounds at which they appear, retained original contig intervals at the five internal transitions, and chromosome-sized pieces left separate. Compare all components to the retained reference assignments; check new cross-chromosome joins and fragmentation of previously supported chromosomes. An early-stop before >10 Mb makes the cap comparison inactive, as happened in the round-2 batch.

Raw hifiasm transitions already internal to an input contig cannot be prevented by withholding a YaHS scaffold adjacency. YaHS can only correct them if its contig EC cuts the sequence. Their assembly-stage prevention is addressed by the separate hifiasm library/post-join experiments already running. YaHS persistence alone is not sufficient support for retention.

Review archives are written next to each task directory under comparisons/ctlk-yahs-round1-JOBID. Native FASTA/BIN remain on Crest; review archives retain AGPs, logs, commands, version, exact input paths and reference validation.
