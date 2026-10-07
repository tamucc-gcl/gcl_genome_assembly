# Automated cutting: deferred

Status: deferred on 2026-10-07 at the user's request. No return date is scheduled.

Current priority is a robust evidence-first manual workflow: detect candidates, create understandable evidence plots/tables/reports and an editable unselected review file, apply only explicitly reviewed cuts on a second run, preserve every base, and recompute chromosome assignment.

Automatic cutting is disabled in the production workflow and standalone applicator. Prototype algorithms and historical experiments are reference material, not validated production authorization. Do not spend further effort adjusting automatic thresholds until this work is explicitly resumed.

Revisit later after the manual workflow is dependable and there is a reviewed calibration set covering supported joins, unsupported gap joins, repeat-obscured boundaries, internal misassemblies, same-chromosome fragments, and contradictory haplotype/library evidence. Assess localization, sensitivity, false cuts, and evidence quality on examples beyond H01. Require human-readable explanations and reproducible source-bound decisions before enabling unattended cutting.
