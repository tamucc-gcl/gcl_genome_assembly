# Source-bound correction interface

The pipeline's reviewed-file route now accepts only the source-bound action format. This route preserves every base, including gap Ns, and assigns neutral identities to split pieces. It does not infer chromosome names from the parent composite. Split pieces are marked `chromosome_assignment_pending` and carried as unplaced until reassignment is implemented.

Use `--chimera_break /absolute/path/to/selected-actions.tsv` with the existing pipeline. Only rows for each assembly are applied. Every selected row must pass validation; malformed rows fail the task rather than being silently skipped.

Required TSV columns:

`id assembly coordinate_stage assessment_sha256 scaffold action cut_bp gap_start gap_end decision_source evidence_packet_id reason`

These are tab-separated columns. `coordinate_stage` is `pre_finishing`; `assessment_sha256` is the checksum of the literal assessed FASTA file. Coordinates are zero-based. `cut_bp` divides the sequence into `[0,cut_bp)` and `[cut_bp,length)`. `UNJOIN_UNSUPPORTED` requires a literal all-N interval `[gap_start,gap_end)` and permits either gap edge or an internal gap cut. Cutting at `gap_end` retains all gap Ns on the left piece.

For an explicitly reviewed internal cut, use `BREAK_PROBABLE_MISJOIN`, add `reviewer`, and set `localization_status=localized`. Rationale and packet identity remain mandatory. Automatic internal cutting is rejected.

Each application writes a coordinate lift and verification JSON alongside the cut audit. The lift records each piece's exact original interval and sequence checksum. Verification reconstructs every parent and checks the written output. Original evidence coordinates remain pre-correction coordinates; existing BAMs must not be attached to renamed outputs without remapping or a validated coordinate transformation.

The automatic applicator guard accepts only decision-stage actions carrying `auto_eligible=yes` and `policy_version=gap-v1`; those fields must never be hand-added to candidate votes. The workflow now supports unattended `--chimera_break auto`: it collects evidence, publishes decisions, and applies only eligible actions. Current collectors supply positive local HiFi support but not the complete matched-control measurements required for automatic gap eligibility. Accordingly, this initial version leaves unsupported/unknown cases unresolved and cannot automatically approve a gap unjoin yet. The policy's accepted-action branch is unit-tested, not biologically calibrated. Full automatic correction and post-cut chromosome reassignment remain implementation work.

