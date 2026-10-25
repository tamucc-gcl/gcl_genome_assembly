# Enabled organelle and taxonomy failure contracts

This pass continues categories 2/3/7 after core pangenome acceptance. It does not
change assembly cuts, contig/scaffolding parameters or the accepted graph policy.

## Issues corrected

1. **GetOrganelle errors were masked.** Both CLI attempts used `|| echo`, allowing
   execution failure to become a successful task with status `failed` and no assembly.
   The pinned GetOrganelle 1.7.7.1 source also catches exceptions/SystemExit and can
   return zero after logging an error. The new command wrapper propagates nonzero
   exits, rejects logged errors/tracebacks and requires the normal completion marker.
   Logs stream to the task output and a sample/organelle-specific command log.
2. **Retry graph discovery assumed plants.** The retry selected only graph names
   containing `extend-embplant`, excluding ordinary animal/fungal SPAdes graphs.
   It now considers the generated FASTG graphs independent of organelle naming.
3. **Circular classification used a substring.** `*complete*` also matches
   `incomplete`. Selection now requires the `.complete.` filename token. A usable
   non-circular path remains linear; verified completion without a path is
   `no_recovery`, which does not imply biological absence.
4. **MITOS failures became empty GenBank files.** A failed executable or missing/
   empty BED now fails annotation. Successful conversion must produce a nonempty
   GenBank. Multi-record inputs are explicitly rejected by this single-record MITOS
   route rather than silently annotating only the longest sequence. The organelle
   assembly is retained, but multi-record annotation is not yet a supported success.
   This is the pipeline adapter's limitation: upstream MITOS 2.1.10 does support
   per-record output directories. Its executable name/options were checked from the
   published wheel, and the previously unpinned MITOS environment now pins 2.1.10.
5. **Missing species context could silently drop samples.** A shared keyed-context
   helper now broadcasts each taxid's reference, genetic code or target list to all
   requesting samples. Missing keys and duplicate context records fail explicitly;
   multiple samples from the same species remain separate. Disabled organelle
   branches still release their nuclear samples without requiring those contexts.
6. **Zero-exit unresolved taxonomy could appear valid.** RESOLVE_TAXONOMY now checks
   for eight rank fields and a nonempty, resolved lineage before emitting its TSV.

These changes distinguish execution failure, unsupported input, verified no-recovery,
linear recovery and circular recovery. A no-recovery result releases an empty bait
set; it does not erase the nuclear sample. Plant organelles retain explicit manual
annotation notes rather than claiming an automated annotation exists.

The ORGANELLE workflow additionally emits `recovery_status` tuples
`(meta, organelle, status_tsv)` for later report integration. Existing bait and
annotation interfaces remain in use. Short-read statuses and full command logs are
published for successful tasks in the organelle output directory. A failed task's
command log and `.command.out` remain in its task directory rather than being
published as a successful result. The consolidated report remains category 6.

## Source verification

The behavior above was checked against the pinned upstream
[reads command](https://github.com/Kinggerm/GetOrganelle/blob/1.7.7.1/get_organelle_from_reads.py)
and [assembly command](https://github.com/Kinggerm/GetOrganelle/blob/1.7.7.1/get_organelle_from_assembly.py).
Their final completion log and exception handling explain why checking return codes
alone is insufficient. The wrapper is conservative: a logged ERROR is not reclassified
as a harmless biological outcome. Warnings and ordinary completed extraction failure
without an error are allowed. Future tool upgrades must revalidate this log contract.

MITOS executable/options and its multi-record output behavior were inspected in the
published [MITOS 2.1.10 wheel](https://pypi.org/project/mitos/2.1.10/). The remaining
multi-record restriction belongs to this pipeline adapter, not to upstream MITOS.

## Tests and practical limits

`tests/run_organelle_contract_checks.py` exercises the actual ORGANELLE,
SHORTREAD_ORGANELLE, ORGANELLE_ANNOTATION, database modules and taxonomy process.
Tests use two samples sharing one taxid to check broadcasting and one output per
sample. The GenBank converter runs with real Biopython/NumPy. Native assemblers,
database downloads, MITOS and circular-plot rendering are deterministic boundary
doubles: this is workflow/format/failure validation, not a biological assembler pilot.

The matrix covers circular and linear recovery, an `incomplete` filename, no-recovery,
animal/fungal retry, both plant organelles, nonzero native failure, logged zero-exit
error, incomplete command log, GetOrganelle/MITOS database failure, annotation failure,
missing BED, multi-record annotation rejection, missing/duplicate targets, missing
genetic code/reference, successful taxonomy, unknown taxid, empty lineage and native
taxonomy failure. `tests/test_organelle_command.py` tests the error/completion rules.
The two-species case also checks that organelle targets and annotation genetic codes
stay associated with the correct sample, rather than leaking across species.

Completed local validation: all 25 matrix cases passed across the batch and targeted
follow-up runs, as did four command-log unit tests. All 20 general core integration
scenarios were rerun successfully, and the accepted production tool/graph/source
audits reused ordinary Nextflow caching. Logs are retained under
`comparisons/organelle-contract-validation/`. No native assembler/annotation accuracy
claim follows from the external-tool doubles used in these checks.

Run the fixture matrix with Python plus Biopython/NumPy and the existing local
Nextflow test runtime. No production reads or reference databases are downloaded
by these fixtures. The local isolated wheel directory is a development convenience,
not a required production input or a new pipeline dependency mechanism.

Categories 2/3/7 remain open pending native pilot evidence, the remaining nuclear
short-read routing/finishing checks and the broader operational retry/dependency
audit. In particular, multi-record animal/fungal annotation requires a separately
designed route before it can be declared supported. Do not describe the current
single-record checks as acceptance of every fragmented organelle input.

## Next Crest check

Resume the normal all-stage QC pipeline into `comparisons/core-operational-01`, with
the same reviewed cuts, graph settings and cache. This validates the real taxonomy
records and existing HiFi organelle context against the updated guards; the Cactus
construction command and graph helpers are unchanged and should remain cached.
The existing five-individual dataset does **not** exercise a native short-read
organelle assembly, so that pilot remains separate and must use actual shotgun reads.

After syncing, run from the Crest project directory:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-operational-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_regional_chunk_bp 25000
```

Return `comparisons/core-operational-01.review.tar.gz`. Inspect the trace for normal
cache reuse and the new taxonomy task's successful resolution. This is a real
pipeline invocation, not the local external-tool-double fixture matrix.
