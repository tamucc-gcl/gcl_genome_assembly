# Chimera context integration checkpoint — 30 September 2026

## Implemented

- `CHIMERA_SEQUENCE_CONTEXT` is a separate module called from the chimera workflow. `chimera_sequence_context=true` by default; set false to disable. Only chromosome-member=yes REVIEW/BREAK_CANDIDATE intervals marked callable or evidence-only receive sequence diagnostics. Assemblies without eligible intervals receive a short status report without mapping.
- Each candidate assembly's intervals are compared directly with other same-taxid pre-finishing assemblies, including its sister when available. Peer identity comes from metadata, not sample-name conventions. The table reports individual relationships and whether a single gapped PAF record brackets the full interval with at least 1 kb on each side. This is descriptive evidence, not a new cutting threshold or a vote.
- Published outputs include raw PAFs, interval sequences, 100-bp tidk windows, coordinate offsets, assessment checksum, peer paths and tool versions, plus Markdown explanations under `assembly/chimeras/sequence_context/<assembly>.sequence_context/`.
- `chimera_hifi_context=false` by default. Enabling it maps the individual's complete HiFi input against the complete assessment assembly once per eligible assembly, retaining sorted BAM/index, regional SAM and MAPQ>=20 primary depth. This new optional path does not import historical Inspector BAMs or standalone comparison BAMs. Its output is cached by Nextflow normally; changing cohort inputs may invalidate this combined diagnostic task, so keep HiFi off for the first integration checkpoint. Future mapping/report separation can improve reuse if this becomes routine.
- Existing chimera evidence now publishes raw tidk windows as declared outputs, and records peak offset instead of implying a distant window peak demonstrates a fusion.
- Added five small tests for chromosome scope, non-callable review candidates, midpoint-only versus interval-wide peer coverage, and required flanks. No tests or Nextflow commands were run locally. Static whitespace review passed.

The diagnostic module does not alter candidate tables or authorize/veto automatic cuts. Existing conservative coordinate/gap checks remain in force. A validated integration of evidence into cutting policy requires a positive error case as well as this no-cut example; do not treat the new reports as an implemented automatic adjudicator.

## Cluster validation after syncing

Keep the original single-library input for this checkpoint. From the project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
```

After successful validation:

```bash
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv assembly
```

Do not delete work or restart assembly from scratch. The new diagnostic tasks and modified evidence tasks need execution; upstream assembly/scaffolding should remain reusable. Verify actual cache behavior from the run log. HiFi context is off, so this checkpoint should not repeat the standalone mapping.

## Return outputs

From the project root after completion:

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
tar -czf "comparisons/chimera-context-${stamp}.tar.gz" \
  genome_assembly/assembly/chimeras
ls -ltrh "comparisons/chimera-context-${stamp}.tar.gz" logs/*assembly* logs/*validate*
```

Download that archive and the relevant validation/assembly logs. If optional HiFi is enabled later, exclude `*.bam`, `*.bai` and `*.sam` from a first transfer and retain them on the cluster. The default peer-only archive contains no new read BAMs.

## Pinned for after broader improvements and the two-library rerun

Revisit CTlk's chromosome arrangement after the other pipeline changes are complete and both Hi-C libraries have been restored. Specifically locate chr12-homologous sequence across every relevant scaffold, distinguish chromosome naming from sequence absence, compare the same intervals between single- and two-library assemblies, and reassess reference/peer chromosome assignments. Preserve the existing comparison archives and mapping BAM as the baseline.

The current local transition has HiFi support and sister interval continuity, and is not dominated by canonical telomere repeats. Do not hard-code that conclusion into executable sample-specific rules. The biological interpretation remains deferred; it need not block the broader revamp.
