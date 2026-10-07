# Evidence-first manual chimera review

Automated cutting is deferred as of 2026-10-07. Both the pipeline and standalone applicator reject auto mode. Prototype automatic decision code remains experimental reference material and is not called by the production workflow. See [deferred automated cutting](chimera-automated-cutting-deferred.md).

Run 1 uses `--chimera_break false`, with chimera detection, sequence context, and evidence enabled. HiFi context is enabled by default; an exact retained BAM manifest avoids remapping. Explicitly disabling HiFi context leaves that evidence unavailable and must be visible in review. The pipeline publishes:

- `assembly/chimeras/review/README.md`: review entrypoint.
- `assembly/chimeras/review/chimera_review.tsv`: one merged editable file.
- Per-assembly reports with peer chromosome tracks, independent-individual summaries, per-library contacts, controls, HiFi measurements, graph context, IGV links, and measurement audits.

Every row starts with `selected=NO`, blank reviewer and blank reason. Copy the published TSV to an input location before editing. Preserve unselected rows and use them to record retained or unresolved decisions. A proposed coordinate is evidence for inspection, never permission to cut. Non-gap transitions remain unlocalized with blank cut coordinates.

Required columns in the forward review file are:

`selected id assembly coordinate_stage assessment_sha256 scaffold action cut_bp gap_start gap_end decision_source evidence_packet_id reviewer reason localization_status`

Additional generated columns `evidence_summary` and `report_path` help review and are carried in the file. Selection must be exactly YES or NO. To approve a literal-gap cut, set selected=YES, keep action=UNJOIN_UNSUPPORTED, set reviewer and reason, and verify coordinates. To approve an internal cut, set action=BREAK_PROBABLE_MISJOIN, supply an exact cut_bp, and set localization_status=localized. Keep decision_source=review. Edit the cut rather than accepting an interval midpoint by default.

Run 2 repeats the original pipeline launch with the same input assemblies, parameters, launch directory, and work directory, adding `--chimera_break /absolute/path/to/edited-review.tsv -resume`. Leave evidence settings unchanged so cached evidence tasks can be reused. Do not rerun from the finalized renamed FASTA. The review checksum binds to the original pre-finishing assessment FASTA.

The selected file is a staged input to cutting, so edits invalidate the cutting stage and its dependent tasks. Only assemblies with selected rows enter cutting. Other assemblies pass through. A species cohort containing actual cuts reruns chromosome harmonization using the same reference. The pipeline's downstream dependencies rerun as needed; -resume does not guarantee zero upstream reruns if inputs, tools, configuration, or scripts have changed.

Validation rejects malformed selection values, selected unknown assembly IDs, stale assessment checksums, invalid/out-of-range coordinates, false gaps, duplicate/conflicting selected cuts, and cuts creating undersized pieces. Selected rows require reviewer and rationale. Unselected unresolved rows can contain blank coordinates and cannot cause a cut. Only the current review-file interface is supported; no legacy-file conversion is maintained.

All coordinates are zero-based, half-open on the original assessed FASTA. A cut divides [0,cut_bp) and [cut_bp,length). For an all-N gap [gap_start,gap_end), cutting at gap_end preserves all gap Ns on the left. Every base is preserved; no sequence or gap trimming is inferred. Pieces initially get neutral names, then chromosome assignment is computed from corrected sequences. Evidence stays in original coordinates, with a coordinate lift and reconstruction audit for each applied cut. Do not attach original BAMs directly to corrected renamed FASTAs.

The manual application and reassignment previously passed real-data H01 testing. The new report/index/manual-selector wiring is locally verified but requires its own Nextflow execution on Crest. The existing real test launcher now runs evidence-only and reviewed-manual cases; it no longer runs automatic cutting.

The evidence pass now publishes Markdown per-assembly report.md packets, a cohort README.md, assembly_registry.tsv including zero-candidate assemblies, static chromosome-track SVGs and cut-instructions.md. The main report integration is deferred to its separate reporting update. User-added cuts need not match detected candidates; use the registry identity and original coordinates.
