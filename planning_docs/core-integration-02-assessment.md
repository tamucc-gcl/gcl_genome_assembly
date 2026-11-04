# Core integration 02 — accepted 2026-11-04

Source: `core-integration-02.review.tar.gz`. Pipeline and audit exit status 0;
explicit planned-rebuild acceptance mode.

## Verified results

- 597 completed and 326 cached tasks, with no failed, ignored or incomplete tasks.
- All five hifiasm assemblies remained cached. Cactus legitimately rebuilt after
  the YaHS fix: 14h55m50s execution, 167 GB peak RSS.
- The native audit verified all ten final FASTA contents and FAI fingerprints
  against the accepted operational baseline. Final name maps, cohort admission
  and reference identity match. FASTA hashes were checked on Crest; the review
  archive excludes large FASTAs and does not support an independent local rehash.
- CLim hap1 is restored to its accepted final assembly. The previous increase
  of 642 records is resolved; the bounded YaHS replay restored its accepted AGP.
- CTlk hap1 has three selected cuts and hap2 one. Both verification files confirm
  exact parent reconstruction in pre-finishing coordinates.
- Graph clipping and sequence/coordinate accounting pass. Source audit remains
  REVIEW with 6,182 inputs, 1,330 names absent from both representations and five
  chromosome-scale retained-path review flags, matching baseline dispositions.
  These flags do not by themselves establish biological loss.
- All five intended within-sample pairs have dotplot and ribbon outputs marked
  plotted. CLim's rendered dotplot was inspected and shows chromosome-wide diagonal
  correspondence. The combined Markdown renderer completes in six seconds.
- All task wrapper contexts are available. Thirty explicit version metadata issues
  remain, chiefly unlabeled producer lines; environment snapshots retain package
  identities. Their labeling is carried into the reporting pass. Windows extraction
  needed a repair for an archive hardlink, not a Crest execution failure.

## Category 6 reporting work

1. All ten embedded graph image references use paths relative to their source
   reports and do not resolve from assembly_report.md. Images exist in the graph
   output directories. Rebase embedded paths or link standalone reports.
2. Imported full-document headings disrupt the main report hierarchy.
3. The chimera summary calls original pre-cut candidates pending review despite
   four verified applied cuts. Show original evidence, recorded CUT/KEEP decisions,
   applied actions and final assembly outcomes distinctly.
4. Pipeline release/commit appears as unknown; use recorded run provenance.
5. Resolve version labels from retained producer/environment context while keeping
   genuinely unavailable measurements explicit.
6. Preserve prior requirements: alternate HiFi placements, R/ggplot2 where practical,
   readable evidence, GitHub-relative links and omission of disabled sections.
   Presentation changes must preserve accepted evidence calculations.

Category 7 operational integration is accepted. Category 6 remains open. Develop
and inspect the report against this accepted packet before a cached native report
check; report development alone does not require another graph rebuild. PSMC is
category 8 and remains the first extension after the core reporting gate.
