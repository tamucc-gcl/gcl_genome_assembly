# Maintained repository and retirement policy

## Maintained sources

| Location | Purpose |
| --- | --- |
| `main.nf`, `nextflow.config` | Production entry point and configuration |
| `functions/`, `workflows/`, `modules/` | Shared contracts, active workflow composition and processes |
| `py_scripts/`, `r_scripts/` | Maintained helpers and renderers referenced by production |
| `environments/` | Task dependency specifications |
| `vendor/cactus-gfa-tools-regional/` | Required pinned regional-assignment correction, source manifest and licence |
| `docs/` | Current interface, code ownership/retirement, verification and tool rationale |
| `planning_docs/` | Canonical checklist, retained acceptance evidence and agreed extension designs |
| `tests/` | Current contract/unit/known-answer tests and bounded integration fixtures |
| `scripts/comparisons/` | Removable acceptance launchers and controlled replay/input helpers, not production modules |

Scientific outputs under `comparisons/`, scratch contents and project handoffs are
ignored run artifacts, not maintained pipeline source. This cleanup does not delete
those outputs or external project data, caches, reviewed cut sheets or assemblies.

## Removed in the core freeze

398 obsolete tracked files were retired: superseded chimera detection/reporting
implementations, duplicate/uninvoked workflows, uncalled pangenome prototype modules,
their helpers/tests, old launchers, backup sources, obsolete plans and rendered report
mockups. Twenty-two unused resource-label blocks and 46 unused pangenome defaults
were removed. Ignored backup/preliminary-copy directories were removed after checking
their location and contents. The active vendored correction is **not** an obsolete patch.

The last pre-cleanup revision is
[`341193c18878886ec18e26c3c843a44507e6162a`](https://github.com/tamucc-gcl/gcl_genome_assembly/tree/341193c18878886ec18e26c3c843a44507e6162a).
Git history preserves retired sources and investigation narratives; do not keep another
executable `deprecated/` tree. Retained acceptance records link historical evidence at
that revision where necessary. The generic SV known-answer fixture generator remains
for the planned method validation; it does not enable unsupported graph SV analysis.

The restored pipeline sketches in `assets/` are retained for a separate later
update. The [retired-component inventory](retired-components.md) records original
purposes, actual replacements and deferred functionality for all 79 reviewed
module/helper/workflow sources: 52 retired and 27 retained prototypes.

## Retained optional-analysis prototypes

Twenty-seven sources remain in their original module/helper locations at the
user's request for later reuse. They cover variant classification, inversion
rescue/rearrangements, private/accessory evidence, progressive graphs and 2D
layouts. They are not called by the production workflow and are not accepted
extensions. The inventory identifies each file. Keep them until those roadmap
categories are reviewed; do not mistake their presence for a working interface.

## Rules for future changes

- Keep specimen IDs, reviewed cut sheets, hardcoded cluster work hashes and project
  biological settings outside general production implementations and launchers.
- Put temporary investigations under the project's ignored `comparisons/` directory.
  Promote only general, tested evidence calculations into the pipeline.
- Remove replaced implementations and their obsolete tests once the replacement is
  accepted. Preserve biological rationale in maintained documents and Git history.
- Keep optional extensions explicit in the roadmap. Uncalled code is not proof an
  extension works; rebuild and validate it when that extension is implemented.
- Keep the source-bound manual-cut and final-coordinate output contracts stable.
  Evidence generation and plotting remain separate. Plotting changes preserve
  retained measurements, and revised figures use R/ggplot2 where practical.
- Run `python3 scripts/check_repository.py` for maintained includes, staged sources,
  documentation links and the frozen parameter-default inventory. A deliberate
  parameter change requires interface review and `--write-defaults`.

The generated default inventory is a review aid, not an automatic semantic validator
for every advanced parameter. Scientific regression and native execution checks remain
necessary when an implementation or environment changes.
