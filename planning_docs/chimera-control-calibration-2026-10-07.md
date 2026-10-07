# Control calibration and peer anchor audit

The control collector now selects up to 48 matched-length literal scaffold gaps and 12 gap-free sequence pseudo-junctions per candidate. These populations are explicitly tagged and displayed separately. Pseudo-junctions are spaced by 500 kb, avoid all suspect intervals and literal N runs within 250 kb, and require qualified HiFi flanks and at least two narrow HiFi bridges. Gap controls require continuous context from three independent qualified peer individuals and do not require HiFi bridges across artificial gaps.

Contact comparisons match within-flank counts within a factor of two, require 100 within-flank pairs per side and at least five qualified controls from EACH population in EACH library. A conservative support-loss result requires the candidate's three-count allowance ratio below 10% of the lowest control ratio. These are calibration heuristics; automatic eligibility still separately requires candidate HiFi informativeness, localization, peer discordance, and graph conditions. Control selection itself is never approval.

Peer anchors now test lengths 5, 20, and 50 kb at offsets 0, 10, and 50 kb from each interval edge. MAPQ, identity, aligned-base coverage, alternative placements, and target-size filters are retained. Inconsistent informative targets across trials invalidate the combined comparison. Farther anchors measure context, not basewise continuity across the original junction or precise cut localization. A peer_anchor_trials.tsv table records every accepted/rejected flank and interval.

Local audit of the 1511311 packet reassessed 297 peer/interval comparisons and 2,673 anchor trials using the existing PAFs. It recovered 27 continuous-context comparisons, one different-chromosome comparison, and two with unknown chromosome identity; 267 remained uninformative. Rejection reasons across both flanks: insufficient anchor coverage 3,583; low identity 655; accepted 545; no alignment 221; short target 215; low MAPQ 74; ambiguity/CIGAR aligned-base insufficiency 53. All H01 peer comparisons remained uninformative. Only one archived gap control reaches three independent continuous peers. No decisions or cuts were changed by this retrospective diagnostic audit.

Audit outputs: comparisons/chimera-peer-audit-1511311/{anchor_trials.tsv,peer_summary.tsv}. The reproducible audit entrypoint is scripts/comparisons/audit_chimera_peer_packet.py.

Next Crest run measures the new sequence controls and broader gap-control pool using existing BAMs and pairs. After syncing, run from the project root:

```bash
sbatch --array=0 gcl_genome_assembly/scripts/comparisons/run_chimera_real_test.sbatch ctlk_analysis
```

Return auto.review.tar.gz. Inspect control_qualification and per-library control_populations, peer_anchor_trials.tsv, raw contact counts, and source-bound decisions. H01 remains review-only; this calibration change does not automatically authorize the previously reviewed manual cut. No assembly rerun or repeat of the passed manual application test is needed.
