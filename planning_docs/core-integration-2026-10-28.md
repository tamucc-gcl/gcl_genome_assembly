# Operational provenance and combined core check — 2026-10-28

The finalized comparison subgate is accepted. This is the next category 7 batch:
run the accepted assembly/QC/graph and R comparison branches together, with an
always-on software provenance collector. Category 7 remains pending native
integration review. Category 6 is the subsequent unified Markdown reporting pass.

## Fixes and scope

1. Software provenance previously ran only inside REPORTING, which itself needs
   both QC and post-assembly enablement. RUN_PROVENANCE now executes independently
   at the terminal end of the normal main workflow, even with the full report off.
2. The old collector concatenated two-column tool/version files with three-column
   process/tool/version files and embedded headers. It now normalizes these two
   currently produced schemas, retains source identity and reports unlabeled,
   malformed, empty and missing version records explicitly. It does not invent
   a tool/version from an unlabeled string or hide it as a successful measurement.
3. All observed tool/version combinations are retained. The current report no
   longer silently chooses one version per tool. The full report redesign is
   still category 6.
4. The collector records source version-file hashes and available task command/
   wrapper hashes. Resolved conda package name/version/build/channel/URL/hash
   records and local Singularity image checksums are captured when accessible.
   Nested version outputs are associated with their original task wrapper.
   Cached producers keep their original task context, rather than being labeled
   as newly executed software.
5. A small declaration catalog fingerprints pipeline/config/module/workflow and
   environment definitions, exposing unconstrained dependencies and mutable image
   references. This catalog is **not a globally locked environment**. Runtime
   package snapshots describe what actually ran; declarations describe defaults
   for future tasks. The permanent collector covers version-emitting tasks; the
   combined check separately inventories every task work directory in its trace.
6. MitoHiFi and FCS image defaults are now exposed through
   `--mitohifi_container` and `--fcs_gx_container`. Their values remain identical
   to the accepted defaults so this audit does not upgrade tools or invalidate
   accepted assembly/graph tasks. Immutable references or audited local SIFs can
   be supplied to freeze a future fresh run. The existing master/latest defaults
   and unconstrained environments are explicitly not described as immutable pins.
7. README's stale Cactus v3.1.4 labels now match the configured v3.2.1 image.
   Existing separate Cactus-package/Toil verification remains authoritative; the
   misleading original construction `cactus --version` output is not relabeled
   as a measured Cactus release.
8. The comparison review launcher now packages provenance, R comparison artifacts
   and acceptance results, excluding raw/filtered PAFs and other large mappings.
   Symlink publishing is preserved. No cache reconstruction was added.

Production artifacts are `provenance/software_versions.tsv`,
`software_version_issues.tsv`, `run_provenance.json`,
`runtime_environments.json`, and `dependency_catalog.tsv`. The normalized version
table is also published at `reports/software_versions.tsv` for the current report.
Image checksums are computed once per distinct local image within each inventory;
whole assemblies, graphs and reads are not rehashed by this collector.

## Operational decisions retained

The general retry policy remains bounded at one retry; validated deterministic
failure routes retain their existing terminate rules. There is no global ignore
policy. The optional PANGENOME_2D_VIZ module has a local best-effort ignore rule;
it remains disabled in this core check and any future enabled failures must be
shown explicitly in the report. Core assembly, mapping, graph and evidence failures
must not be converted into successful no-result states.

Observed resource use is evidence, not a universal ceiling: accepted graph
statistics peaked near 88 GB and coverage/similarity near 66 GB, supporting their
128/96 GB requests. Pairwise mapping used 15–34 GB in the accepted fish cohort,
but earlier larger genomes justify retaining its existing higher memory request.
No broad resource reductions or package upgrades were made from one dataset.

## Local verification

Five provenance tests cover mixed current schemas, multiple versions, explicit
unavailable records, empty inputs, nested task context, actual package metadata
and local image fingerprints. Three acceptance-audit tests cover successful
cached integration with an existing REVIEW disposition, ignored failures/graph
recomputation, and changed final inputs/cohort admission. All eight passed.
Actual Nextflow production-collector cases passed with nonempty and empty version
channels. The full main.nf input-only entry compiled/executed. Shell syntax and
whitespace checks pass. The added trace fields are recognized by Nextflow 23.10.1.
The local resource-instrumentation smoke check hits the existing WSL /proc I/O
sampling limitation (`Unexpected: unbound variable` in the generated wrapper);
native SLURM trace/resource acceptance remains part of the combined run, not a
claimed local success. Resource tracing already worked in the accepted Crest runs.

## Native combined check

From the Crest project directory after syncing:

```bash
mkdir -p logs comparisons
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-operational-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/core-integration-01 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

The launcher uses normal main.nf and its existing `.nf/assembly` / `work` resume
state. It enables all-stage QC, accepted manual cuts, the core pangenome, optional
post-assembly R comparisons within samples, and the current report together.
The sample sheet, manual sheet and baseline paths are project-specific launch
arguments, not hard-coded production rules. The BAM manifest is the existing
test reuse optimization; a fresh production run does not require it.

The final audit compares the current finalized membership/paths/name maps and
reference coordinate identity to the accepted baseline. Cactus construction must
remain cached. Graph accounting must pass; existing source-review dispositions
must remain unchanged. Pair membership and figure/status completeness are checked.
Final software/report artifacts must exist. The trace must not contain failed,
ignored or incomplete tasks. All resolved task environments are inventoried
read-only from the trace; unavailable context stays explicit. Required-input
failures yield a structured acceptance failure rather than disappearing output.

Return `comparisons/core-integration-01.review.tar.gz`. Inspect
`core_integration_check.md/json`, `execution_environments.json`, trace and
`results/provenance/`. Software metadata issues are surfaced for review, not
silently counted as measured versions. This check does not independently hash
the whole graph or establish biological completeness. No new cuts, purging
decisions, graph tool changes or automated biology decisions are introduced.

After native acceptance, close the applicable category 7 integration gate and
move into the agreed run-scoped Markdown report. Reporting will consume the
explicit current-run artifact channels, keep the interactive QC ranking HTML,
omit disabled sections and migrate active Python renderers to R where practical.
PSMC remains the first extension after those core gates.
