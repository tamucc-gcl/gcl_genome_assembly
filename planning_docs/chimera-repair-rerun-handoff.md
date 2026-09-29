# Chimera repair: detection-only rerun handoff

## Status
Changes are prepared for remote validation. No pipeline commands, Python analysis, compilation or tests were run locally. Static source review and Git whitespace checks are the local checks.

This checkpoint tests the repaired detector on the same single-library inputs. Keep chimera_break=false, pangenome off and post-assembly analyses off. Do not restore the second Hi-C library yet. A successful run is evidence for reviewing decisions, not blanket approval of automatic cutting.

## Changes
- Exact inferred chromosome membership is exported as chromosome_member in harmonization name maps and candidates. The existing chromosome-set report retains selection-method/quality context. Detection and both cutting modes use membership, not an absolute 20 Mb cutoff or voter eligibility.
- The old absolute span/member gates are deprecated and ignored. Relative composite/arm-coverage screening remains.
- Alignment query spans define unambiguous chromosome anchors. Same-chromosome overlapping records count once; conflicting overlaps remain unassigned. Whole AGP components are no longer assigned one majority label or filled from neighbours.
- A native gap is callable only if uniquely inside the transition interval, the interval is sufficiently narrow, and no retained alignment spans the proposed position. Anchor bounds, support threshold and compatible-gap count are written in the call table.
- Minimum anchor length uses chimera_component_min_bp, capped at 1% of scaffold length. chimera_max_join_distance now limits half the maximum native interval width; it does not allow snapping. These confidence parameters are separate from chromosome membership.
- Recovery still audits older joins by exact source-flank alignment and current N-gap checks. It enriches existing candidate intervals only; it no longer generates all-scaffold sliding-window candidates or duplicate profiles. Multiple current/recovered gaps remain unresolved. Newly recovered locations remain REVIEW, never automatically promoted.
- Reports, evidence and automatic-cut routing consume the recovered/adjudicated table when recovery is enabled. Supplied-file cutting checks membership, callable status, diagnostic status, interval/gap agreement, current FASTA SHA256 and an N-gap at the cut.
- Evidence evaluates the exact location without moving it. tidk is the only telomere source; missing rows are explicitly unavailable, and tidk logs are published.
- Comparison helpers live in scripts/comparisons/; collected artifacts live in the cluster project's comparisons/.

## Limits to review before auto
The existing scaffold-level concordance heuristic is retained only for a single supported transition whose two chromosome labels match the two-member candidate. Multiple transitions cannot inherit BREAK_CANDIDATE. This checkpoint does not implement new per-junction, distinct-individual concordance inference. Haplotype votes alone cannot distinguish a private heterozygous rearrangement from an assembly error.

An absent safe gap is REVIEW, not permission to cut an aligned component. The previous 30.20 Mb gap should no longer be selected just because it borders an AGP component. Do not force a historical breakpoint or replacement gap into the new results.

Reference self-assessment remains unavailable without an independent alignment frame. Query-span anchors are conservative alignment evidence, not base-resolved proof of a breakpoint.

## 1. Sync and validate
Push all modified AND new files, including py_scripts/chimera_intervals.py, both new test files, and scripts/comparisons/. Pull on the cluster. From the cluster project root:

~~~bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
~~~

Stop if tests fail. The suite currently contains 45 tests (not executed locally). The tests are synthetic; they do not run assemblers or aligners. Remote Nextflow validation is also required because Python tests do not compile DSL2.

For the main sheet shown in this conversation, the original Hi-C pairs are already present inline. Supply an empty additional-readsets table to avoid inadvertently reintroducing the second library:

~~~bash
printf 'sample_id,library_id,readset_id,hic_r1,hic_r2\n' > data/hic_readsets.none.csv
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
    data/assembly_samplesheet.csv data/hic_readsets.none.csv validate
~~~

Use this only with the current single-library sheet containing the inline Hi-C paths, and with the short-read sample still excluded. Confirm the validation report says one Hi-C library/readset per accepted individual and retains the expected Hi-C-only skip. Do not proceed on a failed or unexpected validation.

## 2. Preserve the completed output, keep cache, then rerun
Do not delete the completed baseline. The current run publishes hard links (publish_dir_mode=link), so moving its output directory preserves those files. The pipeline should be stopped before this move.

The following deliberately refuses to overwrite an existing baseline:

~~~bash
mkdir -p comparisons/baselines
if [ -d genome_assembly ] && [ ! -e comparisons/baselines/single_hic_before_chimera_repair ]; then
    mv genome_assembly comparisons/baselines/single_hic_before_chimera_repair
else
    echo 'Baseline already exists or genome_assembly is missing; no move performed.'
fi
~~~

After checking the move succeeded and validation passed:

~~~bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
    data/assembly_samplesheet.csv data/hic_readsets.none.csv assembly
~~~

This launcher uses Nextflow 23.10.1, -resume, the existing work directory, and chimera_break=false.

Keep work/, .nf/assembly/ including .nextflow/cache, shared databases and the historical tst/genome_assembly_store. Fresh published output does not require deleting any cache.

Unchanged assembly, mapping and scaffolding should reuse cache if inputs, resources and launch context match. Harmonization may rerun, including alignment work bundled into its task, because its script/output changed. Chimera detection, recovery, evidence and reports will rerun. Downstream tasks may rerun when the modified name map is an input even though cutting is off. Check the trace rather than assuming every non-chimera task must be cached.

## 3. Collect the completed results
After SUCCESS, compare the preserved immediately previous run with the repair:

~~~bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_comparison.sh \
    comparisons/baselines/single_hic_before_chimera_repair \
    genome_assembly
~~~

For the original historical comparison:

~~~bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_comparison.sh \
    tst/genome_assembly_store \
    genome_assembly
~~~

Download the two timestamped comparisons/chimera-comparison-*-reports.tar.gz files. Each contains sources.tsv identifying its pair of result directories. Retain the corresponding alignments archives for follow-up. The collector follows published file links and packages existing files; it does not run analyses.

## Acceptance checks
- SUCCESS and no unexplained upstream assembly/scaffolding reruns.
- Same input libraries and reference-selection context recorded; inferred membership visible per scaffold.
- Primary evidence only for in-scope BREAK_CANDIDATE/REVIEW composite intervals. Compact older-gap audits may cover all scaffolds.
- No manufactured transition at the old upstream component boundary; absence of a compatible gap remains explicit.
- Differences between native and recovered call tables are explained; the adjudicated decisions match reports/cutting inputs.
- Multiple plausible gaps, unsupported/mixed alignments and diagnostic midpoints remain non-callable.
- tidk output or an explicit unavailable status, never a hidden sequence-count fallback.
- No cuts in this checkpoint. Evaluate individual calls before enabling auto or providing selected validated REVIEW rows.
