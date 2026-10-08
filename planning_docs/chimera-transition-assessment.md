# Transition assessment and manual review

Candidate detection, independent confirmation, localization and cut authorization are separate results. Detected chromosome identities survive weak local confirmation. Below-threshold peer identities are reported as observations, not qualified votes.

Each review row now records detected_transition, assessment_status, transition_id, related_candidate, preferred_candidate and bridge_status. retention_supported rows use action=RETAIN and remain unselected. They are reported separately from review_required boundaries. supporting_measurement rows refer to a localized transition rather than creating another cut decision. Original measurement IDs, raw evidence and all rows remain available.

A nearby verified gap is associated with a source transition only through the existing source_candidates provenance, the same scaffold and matching ordered chromosome pair. Multiple source transitions or intervening candidates prevent consolidation. Proximity alone never establishes precedence.

Sequence context now preserves exact chromosome-assigned alignment segments within 250 kb of each transition, in addition to the coarse plot bins. Bridge assessment projects segments into original scaffold coordinates, records coverage and chromosome order per independent peer, and looks for reversals or other chromosomes between the transition interval and gap. A peer must cover at least half of the bridge without discordant segments to establish adequate observability. Discordant assignments totaling at least 10 kb in a peer block consolidation; smaller observations remain recorded. These are explicit provisional grouping criteria, not calibrated criteria for automatic cutting. Sparse evidence prevents precedence, even if no substantial reversal is detected.

When the bridge is adequately observed without a reversal, the verified gap becomes the preferred candidate and the source interval becomes supporting_measurement in the same transition group. The gap still requires human review. All selected fields remain NO. Automatic cutting remains deferred.

Validation includes supported consolidation, reversal rejection, sparse bridge rejection and retention classification, plus the existing chimera test suite. The archived run can regenerate reports locally; exact bridge segments can be reconstructed from its published peer PAFs and chromosome labels without mapping reads again. New pipeline runs collect these segments directly.
