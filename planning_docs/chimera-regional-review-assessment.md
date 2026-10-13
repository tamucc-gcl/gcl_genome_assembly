# Assessment of the Crest regional review

Reviewed packet: `chimera-regional-review.review.tar.gz`; job 1532402.
This document records findings, not new cut decisions or implemented fixes.

## Execution and consistency

- Exit status 0; Nextflow SUCCESS. 33 tasks completed, 313 cached, no failures or retries.
- Wall clock approximately 17 minutes. CTlk assessment took 10m57s / 10m02s.
- All ten assembly source digests match the previous production-final packet.
- Existing whole-assembly pairwise alignments and finalization were cached.
- Both CTlk HiFi BAMs were reused; updated discovery and assessment ran.
- Eight assemblies have usable chromosome-scale comparisons. CPla haplotypes remain
  outside the resolved chromosome-scale scope, not assessed negatives.
- Six transitions remain: four CTlk hap1, two CTlk hap2.
- 38 review rows: 27 alternative regions, eight gap options, three original-transition
  rows. Every selection is NO. No cuts applied.
- Actual bridge counts reconcile with the bridging-read tables. Eleven figures have
  opaque white backgrounds; focal chromosome colors are distinct and consistent.
- Chimera report links resolve. Main-summary links to final FASTAs and eligibility
  are not available in this intentionally restricted review archive.

## Candidate findings

| Assembly / event | Finding | Interpretation |
| --- | --- | --- |
| CTlk hap1 C01, chr9/7 | Broad bracket 60.8–64.4 Mb. Gap 63,051,225–63,051,325; zero direct gap-spanning reads; both neighboring read neighborhoods observable. Regional probe 63.051 Mb also has zero bridges at every placement tier. Narrow chromosome tests qualify one individual; broader tests qualify three. Both Hi-C libraries have zero matched controls at the gap. | Still the strongest manual-review gap. The previous SUSPECT classification was lost to the broad-localization veto, rather than changed read evidence. It is not an automatically qualified cut. |
| CTlk hap1 C02, chr7/12 | Broad bracket 94.7–100.7 Mb. Each of six independent-haplotype comparisons has 12.85–14.62 Mb aligned chr7 sequence. Regional loss 95.783–95.788 Mb has zero bridges at all tiers, expected dominant chromosome order in the flanks, incomplete right-flank purity, incomplete left read-neighborhood observability, and no demonstrated Hi-C loss. | The chromosome composite is supported; the old narrow chr12-dominated view was misleading. The 95.783–95.788 Mb loss warrants contextual review, but is not an established cut. |
| CTlk hap1 C03, chr4/14 | Original interval 35.568–35.649 Mb has minimum six confident spanning molecules. Merged weak region 35.217–35.496 Mb; probe 35.426 Mb has zero confident bridges but 12 low-confidence primary bridges, all with competing placements. Hi-C within-left counts 4 / 0. | Upstream loss is measured, but placement ambiguity and poor contact observability prevent treating it as a straightforward unsupported join. |
| CTlk hap1 C04, chr2/5 | Alternative 42.831–42.832 Mb has one confident bridge and three independent separating individuals. Neither library demonstrates contact loss at this probe. Original 42.855–42.856 Mb has minimum nine bridges and contact loss in both libraries. | Nearby read and Hi-C evidence point to different locations. A real weak point was found, but it does not yet justify selecting a cut. |
| CTlk hap2 C01, chr12/7 | Original 43.069–43.077 Mb has minimum 23 confident bridges. Ex3 demonstrates contact loss; Ex2 does not. Only nominated alternative 43.489–43.490 Mb is not chromosome-localized. | No persuasive new cut site. Keep the library disagreement visible. |
| CTlk hap2 C02, chr14/4 | Loss 39.864–39.894 Mb; probe 39.865 Mb has zero confident bridges but 311 low-confidence primary bridges, all with competing placements. All-primary minimum across the region is nine. Ex2 has 507 / 1673 within-side pairs but only two matching controls; Ex3 has 87 / 819 and no matching controls. | The trough is real for confident placements. It does not mean the region has no reads. Contact calibration remains insufficient, and the loss lies within a mostly chr4-assigned region. |

## Remaining decision and presentation defects

1. Broad-localization warnings override all site-specific recommendations. This
   demotes hap1 C01 and expands its plot despite useful exact-gap evidence. Keep
   event-level chromosome-composite confidence separate from site-level precision.
   A failed narrow flank test must not erase existing evidence for a manual candidate.
2. `gap_association` accepts a corridor assigned entirely to either focal chromosome
   as consistent. For hap2 C02R01, this yields "confirmed" even though no individual
   has qualified different-chromosome flanks at that location. Compatible corridor
   order is not a localized chromosome transition; those conclusions need separate
   labels and separate use in recommendations.
3. Alternative ranking includes known gaps and same-chromosome weak sites. Hap1
   C02's first-ranked alternative is a gap near 100.07 Mb without qualifying
   separating flanks, while the loss near 95.78 Mb is closer to the expected switch.
   Discriminating chromosome context must outweigh generic minimum-read counts.
4. The same physical gap is represented by both a G option and an R assay row.
   Consolidate physical review options while retaining their distinct measurements.
5. A broad original interval compares contacts across up to six Mb; an alternative
   uses point flanks. Those are different biological-distance assays and their
   calibration outcomes should not be treated as equivalent site measurements.
6. All six overview recommendations are REVIEW, hiding the differences above.
   Provide practical manual-review prioritization and reasons without implying
   automatic cutting or treating inadequate controls as support for retention.
7. Report size is excessive: cohort README 976 lines, hap1 report 898 lines.
   Retain a concise decision table and put routine excluded locations/control
   diagnostics in expandable detail. Some figure peer labels are clipped and
   alternative labels overlap the chromosome-panel title; tighten layout validation.

The next work should address these bounded classification/ranking/presentation
defects using this packet. Existing read mappings and pairwise alignments remain
useful; another assembly or mapping rerun is not required to diagnose them.
