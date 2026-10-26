# Core operational run 01 assessment

Reviewed 2026-10-25 from `core-operational-01.review.tar.gz`.

## Outcome

Accepted as a successful resume and native taxonomy-guard check. Exit status 0;
857 tasks: 855 cached, two completed, no failures, ignored failures or retries.
The completed tasks were `RESOLVE_TAXONOMY (373251)` (18 seconds including task
overhead) and `ASSEMBLY_RUN_SUMMARY` (2.4 seconds). The log contains no WARN/ERROR
entries. All ten expected haplotypes finalized; the explicitly unsupported
Hi-C-only sample remained a documented skip.

All five HiFi organelle assemblies remained cached. Updated required-context
routing delivered their outputs downstream without invalidating accepted work.
This is not a fresh native organelle assembly test and does not exercise the
short-read GetOrganelle/MITOS branch.

## Accepted results preserved

- Cactus construction remained cached at `fc/fa231d` (original duration 16h16m).
- All 111 packaged pangenome files are byte-for-byte identical to production 02.
  The existing five low-retention source flags remain REVIEW; they are unchanged,
  not new operational failures. The sequence-graph, clipping, variant and regional
  audits remain PASS.
- All 160 common packaged reference, QC and final name-map files are byte-for-byte
  identical to production 02. The final linear reference remains CMat hap1 with
  SHA256 `f9d2e52fcbd832e3af3832603bfe8343e67bad00dae2e8217b978aa9ade04260`.
- Manual-cut verification confirms exact parent sequence reconstruction for both
  CTlk haplotypes. Hap1 retains cuts at scaffold_1:63051325,
  scaffold_9:42831000 and scaffold_5:35426000; hap2 retains scaffold_7:39704000.
  No cuts or biological interpretations were changed by this review.

## Review packaging correction

The archive omitted `results/organelle`, although those files are published by
the pipeline. Consequently, this archive validates the organelle tasks' cache
status and downstream completion, but does not allow direct inspection of their
recovery status tables, annotations or figures. The comparison launcher now adds
compact organelle TSV/JSON/Markdown, execution logs, GenBank and image files below
20 MB each. Reads, FASTAs and assembly graphs are excluded. This changes review
packaging only, not task commands, production publishing or cache identity.

Validation of this correction: Bash syntax check and `git diff --check` passed.
No graph rerun is required for the packaging correction.

## Remaining acceptance work

Categories 2/3/7 remain open for the native short-read pilot, remaining nuclear
short-read routing/finishing checks and operational retry/dependency audit. The
current MITOS adapter explicitly rejects multiple FASTA records; implementing
multi-record annotation remains separate work. Full reporting remains category 6.
Do not claim these untested branches passed because this all-HiFi resume passed.
