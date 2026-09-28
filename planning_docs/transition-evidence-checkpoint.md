# Transition-evidence checkpoint

2026-09-28. Prepared after reviewing assembly job 1500517. Static source review only; no pipeline, scientific analysis, compilation or tests executed locally.

## Findings from the completed run

The run succeeded with 267 cached tasks, 23 completed tasks, and no failures or retries. All ten expected haplotypes finalized. Older-join recovery produced all ten audits, including 1,895/1,903 recovered CTlk hap1 joins and 2,004/2,014 recovered CTlk hap2 joins. However, neither CTlk review table contained a selected junction. The local two-flank chromosome test did not cover transitions separated by unassigned sequence.

The absent evidence figures had a concrete source bug: a metric-reporting loop reused the filename-prefix variable, causing PNGs to use a generic filename outside the process output pattern. Candidate verdicts also used the wrong input column in one report field.

## Implemented changes

- Preserve exact-flank recovery and its unresolved statuses.
- Independently scan reference-alignment support in windows of chimera_older_window_bp (500 kb by default), using the existing support and dominance thresholds. Consecutive informative windows assigned to different chromosomes bracket a coarse transition interval; ambiguous/unassigned intervening windows remain part of that uncertainty.
- Publish *.transition_intervals.tsv with interval bounds, chromosome labels, recovered join IDs, selected diagnostic position, profile-selection status and assembly checksum. Coordinates are zero-based, half-open.
- Keep transitions separated by more than chimera_transition_max_bridge (default 5 Mb of intervening windows) in the table without requesting a profile. Retain the existing scaffold size threshold. These thresholds are diagnostic heuristics, not validated biological classification rules.
- Select a transition midpoint for descriptive evidence even when there is no recovered gap. Mark its input row callable=no, evidence_only=yes, REVIEW. Do not snap it to a gap. Report cut_bp=NA and evidence_position_bp separately in the evidence output. Internal input cut_bp is retained solely as the legacy evidence-position transport field.
- Keep recovered join profiles and automatic-cut channels separate. No transition profile authorizes cutting, and chimera_break remains false in the launcher.
- Correct figure filenames and candidate verdict propagation. Shade diagnostic uncertainty intervals on figures. An expected missing figure raises an error; absence of a usable Hi-C profile receives an explicit figure_status.
- Fail if the requested evidence position is absent instead of silently selecting the first candidate row.
- Update Markdown reporting to distinguish review evidence from applied cuts and expose transition/profile request counts.

This does not prove global flank uniqueness, identify base-resolution breakpoints, infer inversion/duplication/translocation calls, or implement independent biological adjudication. Assembly/reference alignment coverage can be absent, especially for the reference itself; an empty transition table is not a clean structural assessment.

## Remote validation and rerun

Sync all modified files and the new tests/test_chimera_evidence_reporting.py. Keep the current testing input sheet without short reads.

From the cluster project root:

~~~bash
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv validate
~~~

Expect 35 tests. The new evidence reporting tests mock numerical and plotting libraries so base Python need not have those packages installed. These test control flow and filenames, not scientific accuracy or actual PNG rendering.

After both checks pass:

~~~bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv assembly
~~~

Keep work and .nf/assembly intact. Recovery, evidence and reports should rerun. Assembly, mapping, scaffolding, gap filling, telomere extension and finalization should stay cached with unchanged inputs/configuration. The full R report is changed for compatible diagnostic wording but is not exercised when the checkpoint disables it.

## Return for review

Collect the same logs/reports as for job 1500517, adding:
- assembly/chimeras/older_joins/*.transition_intervals.tsv
- Current *.chimera_evidence.png and *.chimera_evidence.tsv
- Updated chimera_coordinate_status.md and *.older_join_summary.json

Check ten recovery summaries, visible CTlk transition intervals, explicit profile exclusions, correct figure filenames, retained REVIEW verdicts and no automatic cuts. A transition can remain uncertain or have no uniquely recovered join; that is a reported limitation, not a reason to silently clear it.

Biological cutting decisions, library-specific evidence, reference multiplicity policy, short-read validation and the pangenome/PSMC revamp remain separate subsequent checkpoints.
