# Core interface freeze and source cleanup — 2026-11-05

The core now declares `1.0.0-core`; its supported interface is
[core-interface.md](core-interface.md). This is a core contract, not a claim that
the planned optional extensions have been implemented.

## Changes

- Replaced the obsolete 1,200-line README with current entry points and supported
  routes; documented input identity, manual cuts, coordinates, outputs, cache
  behavior, failures and unsupported extensions.
- Captured 300 configured defaults in a checked inventory. Removed 46 unused
  pangenome settings and 22 unused resource-label blocks.
- Removed 398 obsolete tracked files, plus ignored backup/preliminary-copy
  directories. Retired sources are available in Git history; there is no second
  executable deprecated tree. Scientific outputs, project data and caches were
  not deleted. The active vendored correction, licence and rationale remain.
- Preserved the canonical 13-category checklist, agreed extension requirements,
  PSMC design and relevant native acceptance evidence. Historical links resolve
  at the pre-cleanup revision. PSMC remains first after the core gate.
- Generalized both launchers: no embedded taxid, enzymes, Inspector thresholds,
  organism-specific tuning or fixed optional Hi-C-table requirement. The direct
  launcher delegates to the checkpoint launcher. Explicit caller flags survive.
- The removable integration helper captures supported settings from its supplied
  baseline provenance into `baseline-params.json`, rather than silently relying on
  former project-specific launcher constants. Removed knobs and invocation-specific
  input/output/cut/BAM-override paths are not restored; the caller supplies those.
- Cooler and tidk now emit labelled, single-line version observations even when
  native version output spans several lines. The collector remains strict; no
  legacy-parser or bare-number guessing route was added.
- Added a maintained-source/document-link and default-drift checker, plus launcher,
  integration-parameter and exact producer-statement checks. R sources now have
  explicit LF checkout semantics alongside the existing executable source rules.

Twenty-seven optional-analysis source files are retained unchanged for later
reuse, covering variants/rearrangements, private/accessory evidence, progressive
graphs and 2D layouts. Their production wiring/defaults are not enabled by this
restoration; see the inventory for exact files and remaining validation.

No biological threshold, alignment policy, assembly correction, plotting calculation
or graph preprocessing rule changed. Some previously ignored prototype toggles now
fail explicitly instead of suggesting that a requested analysis ran.

The pipeline sketches in `assets/` remain for a separate later update. See the
[retired-component inventory](retired-components.md) for replacement versus deferred
functionality; removing a prototype does not cancel its roadmap requirement.

## Checks completed locally

- 321 general Python tests pass after removal and launcher/configuration changes.
- 176 maintained Nextflow includes, literal staged-source references, 30 current
  local documentation links and the 300-default snapshot pass repository checks.
- Nextflow 23.10.1 main input/compile, full-cohort chimera wiring, two-species graph
  construction/role/sequence/source-audit contracts and multiple-species main-report
  wiring pass. Controlled tools test contracts, not biological accuracy.
- The graph contract compiles the retained vendored correction with a real C++
  compiler; the actual accepted upstream differences were not removed or replaced.
- Exact production version statements injected from the maintained modules run
  through Nextflow and the production provenance collector with zero issues for
  single-line Cooler and multiline tidk output.
- Launch tests confirm no implicit project biology, optional Hi-C handling,
  normalized output paths and explicit taxid/ploidy forwarding. Integration-parameter
  tests include a manifest block before params and reject retired/run-specific keys.
- `git diff --check` and shell syntax checks pass. Rendering calculations were not
  edited in this cleanup; their native acceptance remains detailed-report-01.

## Native confirmation accepted

[Core-freeze-01 assessment](../planning_docs/core-freeze-01-assessment.md): native
SUCCESS, core audit PASS, 875 cached/64 completed tasks, graph construction cached,
543 report links valid and zero version-record issues. All 84 selected baseline
evidence files were byte-identical. The command below records the accepted test;
it is not an instruction to repeat it before starting PSMC.

## Accepted cached Crest confirmation command

Sync deletions as well as modified/new files. Keep the project-level reviewed cuts,
accepted results, `.nf/assembly` session and work directory intact. Do not copy an
obsolete source tree back over the cleaned repository.

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
mkdir -p logs

CORE_INTEGRATION_ALLOW_REBUILD=false \
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/detailed-report-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/core-freeze-01 \
  --run_report true --report_chromosome_painting true \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

The baseline helper preserves this project's recorded supported settings while
the general launcher stays organism-neutral. The BAM manifest remains an explicit
test-reuse override, not a requirement for a fresh production run.

Expect tidk/contact-map producers and dependent figures/QC/report/provenance tasks
to rerun because version-emission script hashes changed. Assembly, scaffolding,
read mappings, manual cut results and graph construction should remain cached.
The post-run core audit checks this after execution; it is not a scheduler veto.

`CORE_INTEGRATION_ALLOW_REBUILD` is a test-wrapper environment setting, not a
Nextflow pipeline parameter. Its default is `false`, so the explicit prefix can
be omitted unless overriding an inherited value. `false` requires graph construction
to have stayed cached and baseline identities/locations to match at the final audit.
`true` explicitly permits a planned rebuild while retaining content and acceptance
checks. Neither value changes Nextflow scheduling or cache behavior. The shell
prefix passes the setting through `sbatch` to the nested test wrappers; the trailing
backslash only continues the command onto the next line. Production runs do not
need this setting.

Accept launcher exit 0, core audit PASS, report-validation PASS and no new malformed
version observations. Old cached files are not rewritten. Return
`comparisons/core-freeze-01.review.tar.gz`. This confirms the cleanup on Crest; it
does not reopen the accepted chimera biology or ask for another graph rebuild.
