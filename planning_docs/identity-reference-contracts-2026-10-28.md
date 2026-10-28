# Identity and reference contracts — 2026-10-28

## Changes

The parser previously accepted scientific-name-only samples but `main.nf` only
scheduled taxonomy resolution for samples already carrying a taxid. Samples could
therefore lose required downstream species context. `INPUT_TAXONOMY` now runs
before the declared assembly/trait side channels and taxonomy consumers.

- Existing positive-integer taxid metadata and read handles pass through unchanged.
- Name-only requests are resolved together by one streaming pass through the local
  NCBI `names.dmp` and a node-existence check. Matching is exact scientific-name
  matching after case/whitespace normalization; common names are not guesses.
- Unknown names, multiple matching taxids, missing nodes, invalid taxid syntax and
  absent identities fail explicitly. No result is fabricated or silently dropped.
- Name-only metadata receives the canonical name and taxid; the original query is
  preserved in `taxonomy/resolved_names.tsv`. A supplied effective taxid remains
  authoritative, preserving the existing explicit name/search override policy.
- The numeric-only route schedules no name-resolution task. Required-context joins
  validate completeness/uniqueness for the name route. No read copying or new
  metadata fields are added to accepted numeric-ID assembly tasks.

`FIND_MITO_REFERENCE` no longer takes the first output from potentially empty or
ambiguous globs. It requires exactly one nonempty FASTA/GenBank pair with matching
filename stems, one FASTA record, the configured minimum sequence length and an
organism annotation. These are output integrity checks, not independent sequence
identity verification or a guarantee of the best biological reference.
The [upstream findMitoReference source](https://github.com/marcelauliano/MitoHiFi/blob/master/src/findMitoReference.py)
uses matching accession stems for FASTA/GenBank outputs and defaults to one
reference, consistent with these checks. This source inspection does not pin the
container version installed on Crest; the native check remains necessary.

Reference retrieval now requests one CPU, 4 GB and two hours, while retaining the
same installed MitoHiFi reference tool. Mitochondrial assembly resources are unchanged.
`run_post_assembly=false` is now initialized explicitly; an absent parameter no
longer warns or implicitly enables optional comparisons. Existing launchers already
supplied false, so this preserves their effective setting.

## Retry boundary

| Work | Policy | Reason |
|---|---|---|
| Local name resolution | terminate | Same snapshot/request cannot fix an unknown/ambiguous name by retrying |
| Offline taxonkit resolution | terminate | Invalid taxid, missing dump or local lookup failure needs correction |
| Missing/invalid downloaded reference outputs | exit 65, terminate | No usable unambiguous pair exists |
| Unavailable/nonexecutable reference command | exits 126/127, terminate | Dependency failure needs correction |
| Other external reference-tool failures | one bounded retry | Tool exit alone does not distinguish transient network failure from other external errors |

This is a targeted lookup/reference policy. The pipeline-wide default remains one
retry; this pass does **not** claim every assembly/QC process now distinguishes all
deterministic/transient failures. Database downloaders also retain their existing
internal network retries. Further dependency/resource/retry audit remains category 7.
The current MitoHiFi container uses a mutable `master` tag; pinning/acceptance is
still part of that dependency audit, not silently claimed complete here.

## Local evidence

Seventeen workflow scenarios passed using actual identity/reference modules and a
controlled remote reference executable. They cover mixed/known-only records,
unknown/ambiguous names, missing nodes/identity, invalid taxids, a valid pair,
missing/empty/mismatched/multiple/short/header-only references, absent organism
annotation, unavailable executable and external failure. The last case proves one
retry; deterministic reference failures prove no retry.

The mixed-identity case was rechecked after adding canonical-name propagation.
The entire main pipeline compiled and its input-only entry completed successfully.
The native-check workflow also completed locally with controlled external tools.
These checks do not establish that Crest's installed tools or remote NCBI access
pass; the small native check below is the remaining evidence for this chunk.

The local full-main check disables Linux trace instrumentation, as other local
fixtures do; one attempted run with instrumentation enabled failed inside the
Nextflow wrapper, not inside the input-summary command. Successful native Crest
traces remain the evidence for production trace behavior. Independent Java runs
must not concurrently unpack the same Capsule runtime; the final compilation was
run serially after the matrix.

## Crest native acceptance

Sync the repository, then from the project directory:

```bash
mkdir -p logs comparisons
sbatch gcl_genome_assembly/scripts/comparisons/run_identity_reference_checks.sbatch \
  comparisons/core-identity-reference-01 \
  373251 'Spratelloides delicatulus'
```

Return `comparisons/core-identity-reference-01.review.tar.gz`.

This removable diagnostic imports the production identity, taxonomy and reference
modules. It verifies two synthetic sample identities (one name, one taxid) against
the actual local taxonomy database, runs installed taxonkit and obtains a real
reference through findMitoReference. It uses no reads and launches no assembly or
graph work. Only the diagnostic entry is copied temporarily to the repository root
to preserve normal project-relative helper paths; its launch/cache/output paths
are separate from the accepted pipeline run. Nextflow handles ordinary resume.

Acceptance requires exit 0, `identity_reference_acceptance.md`, both identity
records, agreement on taxid, a complete lineage record and a validated reference
pair. This is not another native shotgun assembly pilot. The general-purpose
workflow change is already wired into `main.nf`; the diagnostic can be removed
after acceptance without disabling production name handling.
