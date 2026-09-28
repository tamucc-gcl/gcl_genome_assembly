# Older-join review checkpoint and rerun handoff

2026-09-28. Source reviewed statically; no pipeline, analysis code, or tests executed locally.

## Changes ready for remote validation

- Recover first-round scaffold joins directly against the current pre-finishing FASTA with minimap2 2.28. This avoids assuming Inspector/decontamination preserved coordinates.
- Require full-length exact alignments of both flanks, a single reported placement per flank, MAPQ >=30, consistent orientation, and a retained N-gap of the original length. Defaults: 10 kb flanks, at least 1 kb usable sequence.
- Audit every inventoried first-round gap between components. Changed, missing, duplicated, split, or ambiguous flanks remain unresolved. Heuristic alignment uniqueness is not a proof of global uniqueness.
- Assess reference-chromosome support separately on each side of each recovered junction. Overlapping query alignment spans are unioned within each chromosome. These are descriptive alignment-span summaries, not base-level variant calls or a new cross-individual vote.
- Collect existing Hi-C/telomere evidence for recovered chromosome transitions on scaffolds meeting the existing chimera_min_span threshold (20 Mb). All smaller joins remain in the audit; their profile exclusion is explicit.
- Recovered calls are REVIEW only and excluded from the automatic-cut channel. Native last-round calls remain separate. Keep chimera_break=false.
- Composite chromosome names produce a reference warning rather than excluding an otherwise eligible candidate. Cactus accepts their plus signs.
- Repair coordinate-summary ID matching and normalize evidence lookup IDs.

The recovery branch is enabled by chimera_recover_older when round-two scaffolding and harmonization/chimera detection are enabled. It requests 8 CPUs, 32 GB and 8 hours, with at most two recovery jobs concurrently. Evidence profiling has its existing resources. The direct alignment uses minimap2's documented asm5, CIGAR, extended-CIGAR and secondary-alignment options:
https://github.com/lh3/minimap2/blob/v2.28/minimap2.1

## Remote launch sequence

Sync all modified AND newly added files through GitHub. Retain the current sample sheet without the short-read sample. From the project root on the cluster:

~~~bash
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv validate
~~~

Expect 25 tests in this checkout. Stop and inspect any failure. After successful input validation:

~~~bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv assembly
~~~

The launcher retains Nextflow 23.10.1, -resume, disabled automatic cutting, disabled pangenome/post-assembly analyses and qc_mode=none. Validation has a separate launch directory and cannot replace the assembly resume target. Input validation is not a substitute for executing the new recovery/evidence tasks.

Keep the existing work directory and .nf/assembly cache. Assembly, Inspector, read mapping, scaffolding, finishing and final FASTAs should reuse cache with unchanged inputs/configuration. New recovery tasks, changed evidence tasks, eligibility and reporting are expected to run. Confirm actual reuse in the trace/log; do not assume it from this prediction.

## Return after the run

Send the Nextflow log, Slurm output, assembly_run_summary.md, assembly_eligibility.tsv, chimera_coordinate_status.md, and:
- assembly/chimeras/older_joins/*.older_join_summary.json
- assembly/chimeras/older_joins/*.older_join_audit.tsv
- assembly/chimeras/older_joins/*.review_joins.tsv
- Updated CTlk hap1/hap2 evidence tables, figures and Hi-C pairs audits from assembly/chimeras/evidence.

Check that all ten expected haplotypes have recovery summaries, no assembly was silently dropped, and unresolved counts are visible. Examine CTlk chromosome transitions individually, including the extra hap1 transition. An empty review table does not clear unresolved joins. Coordinates describe pre-finishing assemblies, not the renamed, gap-filled final FASTAs.

## Remaining decisions and later work

- Cross-individual, per-junction biological adjudication and library-specific Hi-C support remain future work. Sister haplotypes share evidence and are not independent confirmation. No cuts are authorized by this checkpoint.
- Multiplicity reference gate remains unchanged: CTlk hap1 deviation 0.369 exceeds 0.30; CTlk hap2 ranks among candidates. Review whether that gate should become a penalty after checking the score tables. Composite-name eligibility changes do not remove this upstream scoring gate.
- Compare older baseline assemblies at /work/birdlab/GCL/spratelloides_delicatulus_genome/tst/genome_assembly_store by sequence correspondence, not chromosome numbers across different selected references.
- Synthetic fixtures cover coordinate movement, reverse placement, duplicate/split/missing or changed flanks, source-gap disagreement, overlap counting, retained-gap recovery and gap replacement. Real minimap2 behavior and Nextflow wiring require this remote checkpoint.
- Short-read regression testing, PSMC integration and the remaining pangenome revamp follow separately. Annotation stays outside this pipeline.
