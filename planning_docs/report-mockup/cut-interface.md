# Manual cut file interface

Use one TSV for all assessed assemblies: [chimera_review.tsv](chimera_review.tsv). This first-pass file has eight candidate rows across two assemblies, all selected=NO. All ten assemblies are listed in the assessment registry. Tabs separate fields; blank cells are intentional.

## Approve or leave a detected candidate

Copy the generated file outside the output directory before editing. For an approved cut, set selected=YES and supply reviewer and reason. Keep the identity, assembly, coordinate stage and assessment checksum unchanged. Leave retained or unresolved candidates selected=NO; optionally record why in reason.

| Field | How to use it |
|---|---|
| selected | YES applies the cut; NO does not |
| id | Unique identifier for the decision |
| assembly | Exact assessed assembly ID, including haplotype |
| coordinate_stage | pre_finishing |
| assessment_sha256 | SHA-256 of the exact assessed FASTA file, binding the cut to its source |
| scaffold | Original scaffold ID in that assessed FASTA, before final chromosome renaming |
| action | UNJOIN_UNSUPPORTED for a verified all-N gap; BREAK_PROBABLE_MISJOIN for a reviewed internal sequence cut; UNRESOLVED only on unselected rows |
| cut_bp | Integer position between bases; divides [0,cut_bp) and [cut_bp,length) |
| gap_start / gap_end | Required all-N interval for a gap cut; blank for an internal cut |
| decision_source | review |
| evidence_packet_id | Reference identifying the evidence supporting this decision |
| reviewer / reason | Required for every selected row |
| localization_status | proposed_gap for a verified proposed gap; localized for an exact internal cut |
| evidence_summary / report_path | Human context and a link to supporting evidence; retained in the audit |

Coordinates are zero-based and intervals are half-open. A gap [63051225,63051325) contains 100 bases. Cutting at 63051325 retains those Ns on the left piece. The applicator preserves all sequence; a cut is not trimming.

### Select a gap cut

Candidate C02 in CTlk hap1 has cut_bp=63051325 and the verified all-N gap [63051225,63051325) on Sde-CTlk_104_hap1 / scaffold_1. Review its evidence, change NO to YES, and add reviewer and reason if you approve that exact cut.

## Add a cut the detector did not find

Append a row to the same TSV. Detection is not required: the current applicator validates selected rows against the original FASTA rather than requiring a matching candidate call.

1. Copy a row for the correct assembly to retain its exact assembly ID, pre_finishing stage and assessment checksum. If there is no candidate row for that assembly, obtain those values from its coordinate audit/provenance and the original assessment FASTA; do not borrow another assembly's checksum. The [assessment registry](assembly-registry.tsv) supplies assembly identities, including assemblies with no candidates.
2. Give the row a new unique id, for example manual-additional-001. Enter the exact original scaffold ID.
3. Enter an exact cut_bp. For a verified literal gap, use UNJOIN_UNSUPPORTED and its gap bounds. For an internal sequence cut, use BREAK_PROBABLE_MISJOIN, leave gap bounds blank and set localization_status=localized.
4. Record reviewer, reason, decision_source=review, and evidence_packet_id. Link your supporting IGV observations, plots or review note in report_path. A user-added row does not automatically create new candidate-specific evidence assays in the current workflow.
5. Set selected=YES only when the coordinate and evidence have been reviewed. Keep unresolved entries at NO.

[manual-addition-template.tsv](manual-addition-template.tsv) contains a valid unselected example row with a real CTlk hap1 assessment identity and blank cut coordinate. Append its row, not its header. It supplies no cut location and grants no approval.

Multiple cuts on one scaffold use separate rows with unique IDs, all in original coordinates. Do not convert later cuts into coordinates of pieces created by earlier cuts.

## Apply and verify

Repeat the original pipeline command with the same source assemblies and add:

```bash
--chimera_break /absolute/path/reviewed-cuts.tsv -resume
```

The pipeline rejects stale checksums, unknown selected assemblies or scaffolds, invalid coordinates, false gap descriptions, duplicate cuts and pieces below the configured minimum length. It records applied decisions, checks exact parent reconstruction, writes a coordinate lift, and repeats chromosome assignment for affected cohorts. Original evidence BAMs belong to the original assessed FASTA; use the lift when comparing corrected outputs.


