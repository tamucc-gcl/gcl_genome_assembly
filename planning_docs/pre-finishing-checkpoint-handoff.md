# Pre-finishing assessment checkpoint

Prepared 2026-09-24. Source edits and static review only: no Nextflow, assembly,
analysis scripts, or tests were executed locally. Nextflow remains 23.10.1.

## Implemented order

1. Contig assembly/refinement.
2. First Hi-C scaffolding round.
3. Existing Inspector correction and optional scaffold decontamination.
4. Existing Hi-C remapping and optional second scaffolding round.
5. Harmonization assessment; exact last-round join detection/evidence; optional cuts.
6. Gap filling of HiFi/Hi-C assemblies.
7. Optional telomere extension of HiFi assemblies.
8. Final chromosome naming/orientation and configured QC.

ASSEMBLY_FINISHING is the new permanent workflow boundary. Short-read assemblies
pass through it unchanged; HiFi-only assemblies bypass gap filling and can receive
telomere extension. Finishing retains the post-cut name assignments. Finalization
checks sequence-ID correspondence, rejects duplicated names, and refreshes lengths
from the finished FASTA. The combined report name map now uses final maps.

## Coordinate contract

- Round two: use Hi-C pairs already remapped to corrected/decontaminated scaffolds,
  then lift through the round-two AGP only.
- Unmodified round one: use contig-stage pairs and the round-one AGP only.
- Correction/decontamination with round two explicitly disabled: report exact join
  evidence unavailable; do not apply the old AGP to the altered assembly.
- No new read-mapping task was added.
- Liftover retains split-component source intervals and reverses strands correctly.
  Ambiguous duplicated placements are discarded and counted rather than resolved by
  whichever placement is encountered first. Zero supporting pairs is reported.
- Validate AGP object lengths, sequence IDs, and gap intervals against the exact
  pre-finishing FASTA. Called tables carry assembly_sha256, coordinate_stage and
  join_scope. Reviewed break tables must match that FASTA; old tables must be regenerated.
- Assemblies without a usable name map or called table pass through automatic cutting.
- Do not interpret these evidence positions as final-assembly coordinates.

**Coverage limitation:** first-round joins inside corrected round-two components
are not reconstructed. They remain unresolved for exact cutting. Harmonization can
still flag the current scaffold as a candidate. The new Markdown status section
makes this limitation explicit; absence of a callable join is not a clean result.
The reference itself has no self-PAF and is not tested for reference-relative joins.
A complete older-join provenance solution remains future work.

## Other harmonization repairs already included

Failed chromosome-scale voters are no longer promoted because too few passed or
because an assembly was selected as reference. With no candidate passing the
existing criteria, two-pass harmonization emits no reference and assemblies pass
through. A low voter count is reported. This does not introduce new biological
thresholds or make the size heuristic a guarantee of chromosome completeness.

## Expected resume behavior

Keep the existing work directory, Nextflow history, input sheets, and launcher.
Hifiasm/SPAdes, read preparation, correction, both scaffolding rounds, and their
existing read mappings retain their process boundaries and input contracts.
Those tasks should resume when their prior outputs and cache records are intact;
verify this in the new log rather than assuming every task will be cached.

Expect harmonization, chimera handling, gap filling, telomere extension,
finalization, and affected downstream QC/reporting to rerun. Finishing moved to
a new workflow name, which changes its task identity even where input sequence
happens to be unchanged. No automatic output-directory relocation was added.

Wait for the current remote run to finish before syncing. Keep the next checkpoint
at run_pangenome=false, chimera_break=false and the current QC settings.

## Remote validation

The following commands are for the user to run after syncing. They have NOT been
executed by Codex.

From the repository, with Python 3 available, run the synthetic checks:

~~~bash
python3 -m unittest discover -s tests -p 'test_*.py'
~~~

The tests cover split source components, reverse-strand projection, ambiguous
duplicate placement, same-length gap alteration, stale reviewed break tables,
and an empty supporting-pair set. They need no biological data or third-party
Python libraries.

Then, from the project root, retain the existing validation and assembly commands:

~~~bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv validate
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch data/assembly_samplesheet.csv data/hic_readsets.csv assembly
~~~

Run validation first and inspect its result before submitting the assembly resume.
The validation entry checks parsing/startup; it does not exercise all assembly
workflow branches. Remote runtime validation is still required.

Confirm:
- Five HiFi assemblies resume from hifiasm and yield ten finalized haplotypes;
  the short-read sample yields one final assembly; the Hi-C-only sample stays skipped.
- The two libraries for Sde-CTlk_104 remain on the existing mapping route.
- ASSEMBLY_FINISHING runs downstream of harmonization/chimera routing.
- assembly_run_summary.md includes the collapsed coordinate-status section.
- Per-assembly coordinate_audit.tsv files describe pre_finishing / round2_only.
- No original-contig-to-round-two chain is used for chimera evidence.

Before enabling automatic cuts in a production run, validate the changed workflow
with a small cut/no-cut case. Existing reviewed post-finishing tables are deliberately
rejected and cannot be reused blindly.

## Remaining revamp work

The assembly-only checkpoint includes the reordered finishing workflow and an
explicit final-assembly eligibility audit. The CLIP/GREF/FULL construction boundary
has also been drafted, with Cactus pinned at 3.2.1 (the user-confirmed image from
the successful GREF run), but the pangenome rewrite is NOT ready to enable.

The audit records declared ploidy and assembly representation, distinct biological
individuals, missing final members, chromosome-scale suitability, and planned
graph membership per species. Fragmented phased assemblies can be members without
being voters. The short-read primary assembly remains an assembly output but,
under the current declared diploid input, is not a haploid graph member.

Failure isolation is NOT implemented in this checkpoint. The existing retry then
terminate policy remains. Missing-member checks protect cohort selection when an
audit is reached; they do not make Nextflow continue after a required task failure.
A failed run may not produce a fresh summary. Check its run name and final log status,
and never treat an older report or existing FASTA as evidence of completion.

Before pangenome enablement, finish primary Markdown graph-report integration,
published-tool sharing/ordination migration, Cactus 3.2.1 artifact verification,
and safe per-sample failure continuation with an explicit incomplete-run outcome.

Later work also includes the published-tool sharing/ordination refocus and targeted
INV/DUP/translocation tool comparisons. Annotation, BlobTools rebuilding, TellSeq
and ONT remain outside this checkpoint.
## Review the current remote run before syncing

Keep its final Nextflow log and Slurm output. Record its success/failure status,
run name, and revision. Preserve the work directory and .nf/assembly history.

Review the current assembly_run_summary.md, final FASTA/FAI/name-map inventory,
and harmonization/chimera reports. Expected completion for this sheet is ten HiFi
haplotypes plus Sde-CMat_061_primary, with Sde-CTlk_101 explicitly skipped.
Each FASTA should have its matching index and map; no accepted samples should be
missing. QC mode none does not establish biological correctness.

Save copies of the small current reports/logs before resuming for comparison.
Do not delete, move, or clean the work/cache or shared databases. The pipeline
continues publishing directly to the user-selected output directory.

After syncing, first run the synthetic tests and input-validation checkpoint.
If either fails, send the complete error and do not launch the assembly resume.
If they pass, run the existing assembly checkpoint command. It keeps pangenome,
post-assembly extras, and automatic cuts off. This is the first runtime validation
of the changed assembly workflow, not a claim of prior successful execution.

## New-run acceptance checks

- Confirm the log's final status and the run name in assembly_run_summary.md agree.
- Confirm all 11 expected final assemblies are listed for that invocation.
- Inspect assembly_eligibility.tsv and assembly_eligibility.json: five eligible HiFi
  individuals / ten haplotypes, one excluded diploid short-read primary. A graph-ready
  cohort additionally requires an actual passing chromosome-scale reference.
- With pangenome disabled, a qualifying cohort reads disabled_but_eligible; otherwise
  it records why it is withheld. No graph jobs should execute.
- Inspect pre-finishing coordinate audits, candidate/joins tables, and final maps.
  Finishing can legitimately change lengths; final maps use finished lengths.
- Check the log for cached upstream tasks. Investigate unexpected hifiasm/SPAdes
  cache misses before allowing expensive replacements to proceed.
- Compare scaffold identities across stages using their corresponding name maps.
- Leave automatic cutting disabled until the separate small cut/no-cut validation passes.

## Static verification performed locally

Read the changed workflow boundaries and process input/output contracts; checked all
Nextflow include paths resolve, including implicit .nf extensions; searched for
references to removed finishing outputs and the nonexistent statistics versions
channel; ran git diff --check. No pipeline commands, Python analysis, compilation,
or synthetic tests were executed. Runtime compatibility remains to be established
on Nextflow 23.10.1.

Nine additional synthetic eligibility tests cover individual counting, disabled but
eligible cohorts, fragmented non-voting members, collapsed versus declared haploid
short-read assemblies, missing members, species separation, reversible graph names,
numeric taxonomy keys, and duplicate identities.
## Deferred CTlk comparison (2026-09-28)

User-confirmed previous single-library assembly baseline:
/work/birdlab/GCL/spratelloides_delicatulus_genome/tst/genome_assembly_store

The previous reference_id.txt identifies Sde-CMat_203_hap2; the completed current
run uses Sde-CMat_203_hap1. Do not compare chromosome numbers as stable identities
between these runs without reconciling their references.

The user supplied the old chimera_candidates, chromosome_sets, reference_id and
both CTlk chimeric_joins tables. Copies are retained in the task outputs under
previous-single-library-baseline. Defer investigation until the updated assembly
checkpoint completes, as requested. Keep chimera_break=false.

Pending code moves assessment ahead of finishing and replaces the two-round
coordinate chain. It does not yet implement per-junction verdicts or evidence for
REVIEW cases. Those remain explicit follow-up items, along with older joins inside
corrected components that the last-round-only assessment cannot resolve. Missing
calls after the update must not be interpreted as resolution of the CTlk candidates.
### Deferred reference-selection audit

Explicit user request: investigate the change from CMat_203_hap2 to CMat_203_hap1
rather than assuming it is benign. Compare the old and updated candidate lists,
per-candidate score tables, aggregated reference_scores.tsv, chosen_reference.txt,
reference_id.txt and actual downstream reference inputs. Verify exclusions,
missing/failed score tasks, voter composition (CPla must remain a passenger),
score direction, tie-breaking and deterministic ordering. Separate changes due to
the extra Hi-C library and pre-finishing assessment from accidental fallbacks,
stale published tables or reference propagation errors. Confirm the chosen ID is
consistent across harmonization and the future shared pangenome/PSMC reference
manifest, or document an explicit eligibility-driven replacement. A changed winner
alone is neither evidence of failure nor proof that selection is working correctly.

This investigation remains deferred until the new assembly checkpoint results are
available. It does not enable pangenome, PSMC or automatic chimera cutting.