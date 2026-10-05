# Chimera junction tracker

As of 2026-10-05. Sample **Sde-CTlk_104**. All ranges are **0-based, half-open bp in the current pre-finishing assessment assembly**. Seven pipeline transitions (J) plus four historical-switch leads (H), on four scaffolds. No cuts authorized.

| ID | Haplotype | Current scaffold | Current chromosome name | Position range (bp) | Origin | Historical haplotype / chromosome | Evidence summary |
|---|---|---|---|---|---|---|---|
| CTLK-J01 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 62,618,168–62,941,917 | Present in raw hifiasm contig; error mechanism unresolved | old hap2 / chr7_2 (-) | 523,749 bp full window exact. Same local sequence existed historically. Does not establish biological fusion. |
| CTLK-J02 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 95,877,814–96,047,526 | Present in raw hifiasm contig; error mechanism unresolved | old hap2 / chr12_1 (-) | 369,712 bp full window exact. Same local sequence existed historically. Sister-haplotype recurrence is not independent validation. |
| CTLK-J03 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 115,992,478–116,065,811; gap 116,057,447–116,057,547 | YaHS round 1 component join | old hap2 / chr12_1 upstream; unplaced_122 matches parts of inserted component | No continuous historical placement through both neighbours established. First edge of a 113,941 bp component. Physical gap length unknown despite 100 placeholder Ns. |
| CTLK-J04 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 116,171,469–116,221,388; gap 116,171,488–116,171,588 | YaHS round 1 component join | old hap1 / chr12_1 downstream; old hap2 / unplaced_122 partial matches | Downstream 100,694 matching bp plus 1 bp deletion. Second edge of the same component as CTLK-J03. Different local haplotype affinity is unresolved. |
| CTLK-J05 | hap1 | scaffold_5 | chr4_1+chr14_1 | 35,443,804–35,581,476 | Present in raw hifiasm contig; error mechanism unresolved | old hap1 / chr15_2 (-) | 303,862 bp exact, includes full transition and both anchors. Repeat ambiguity is substantial at the left outer anchor. |
| CTLK-J06 | hap2 | scaffold_5 | chr7_1+chr12_1 | 43,055,156–43,233,354 | Present in raw hifiasm contig; error mechanism unresolved | old hap1 / chr12_1 (+) | 378,198 bp full window exact. Same local sequence existed historically. Old hap2 also has a related, non-identical alignment. |
| CTLK-J07 | hap2 | scaffold_7 | chr4_1+chr14_1 | 39,549,503–39,688,365 | Present in raw hifiasm contig; error mechanism unresolved | old hap2 / chr15_1 (+) | 338,862 bp full window exact. Right outer anchor has alternative placement within the current haplotype. |
| CTLK-H01 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 63,051,203–63,051,334 | Unresolved; intersect with AGP before attributing to YaHS | old hap2 / chr7_2 to chr6_3 | Adjacent long-block endpoints, 131 bp apart.. Related scaffold has CTLK-J01, but these are different positions. Historical fragmentation can explain a target switch. |
| CTLK-H02 | hap1 | scaffold_1 | chr7_1+chr9_1+chr12_1 | 94,151,483–95,221,878 | Unresolved; intersect with AGP before attributing to YaHS | old hap2 / chr6_3 to chr12_1 | Broad 1,070,395 bp bracket between long blocks. Smaller alignments not yet reconciled.. Related scaffold has CTLK-J02, but these are different positions. Historical fragmentation can explain a target switch. |
| CTLK-H03 | hap2 | scaffold_5 | chr7_1+chr12_1 | 54,123,484–54,123,676 | Unresolved; intersect with AGP before attributing to YaHS | old hap1 / chr12_1 to chr6_2 | Adjacent long-block endpoints, 192 bp apart.. Related scaffold has CTLK-J06, but these are different positions. Historical fragmentation can explain a target switch. |
| CTLK-H04 | hap2 | scaffold_7 | chr4_1+chr14_1 | 46,160,462–46,160,583 | Unresolved; intersect with AGP before attributing to YaHS | old hap2 / chr15_1 to chr4_2 | Adjacent long-block endpoints, 121 bp apart.. Related scaffold has CTLK-J07, but these are different positions. Historical fragmentation can explain a target switch. |

## Evidence and next checks

| ID | Source components | HiFi | Anchor placements | Local Hi-C | Next check | Status | New evidence / decision rationale |
|---|---|---|---|---|---|---|---|
| CTLK-J01 | h1tg000004l | 0 span entire interval; no retained read is long enough. Longest 32,674 bp. | 1 / 1 | 62.900 Mb; ratio 0.159; inside interval | After YaHS review, refine assignment switch and test narrow junction with competitive HiFi mapping. | REVIEW; No cut authorized | |
| CTLK-J02 | h1tg000133l | 0 span entire interval; no retained read is long enough. Longest 33,258 bp. | 1 / 1 | 95.900 Mb; ratio 0.439; inside interval | After YaHS review, refine assignment switch and test narrow junction with competitive HiFi mapping. | REVIEW; No cut authorized | |
| CTLK-J03 | h1tg000549l_1 (+) to h1tg002235l_1 (+) | 0 primary molecules bracket both 1 kb gap flanks. Absence not calibrated. | 1 / 1 at gap | 115.200 Mb; ratio 0.639; outside interval | Validate component assignment and alternative neighbours; local HiFi opportunity and library-specific Hi-C. | REVIEW; No cut authorized | |
| CTLK-J04 | h1tg002235l_1 (+) to h1tg001566l_1 (-) | 0 primary molecules bracket both 1 kb gap flanks. Absence not calibrated. | 1 / 1 at gap; 1 / 2 outer interval | 115.200 Mb; ratio 0.639; outside interval | Validate component assignment and alternative neighbours; local HiFi opportunity and library-specific Hi-C. | REVIEW; No cut authorized | |
| CTLK-J05 | h1tg000153l | 0 span entire interval; no retained read is long enough. Longest 33,027 bp. | 50 / 1; possible saturation | 35.000 Mb; ratio 0.242; outside interval | After YaHS review, refine assignment switch and test narrow junction with competitive HiFi mapping. | REVIEW; No cut authorized | |
| CTLK-J06 | h2tg000028l | 0 span entire interval; no retained read is long enough. Longest 33,258 bp. | 1 / 1 | 43.200 Mb; ratio 0.486; inside interval | After YaHS review, refine assignment switch and test narrow junction with competitive HiFi mapping. | REVIEW; No cut authorized | |
| CTLK-J07 | h2tg000298l | 0 span entire interval; no retained read is long enough. Longest 31,142 bp. | 1 / 2 | 39.700 Mb; ratio 0.243; outside interval | After YaHS review, refine assignment switch and test narrow junction with competitive HiFi mapping. | REVIEW; No cut authorized | |
| CTLK-H01 | Not established | Not assessed at this position | Not assessed at this position | Not assessed at this position | Locate verified AGP component boundary and reconcile smaller/alternative alignments; then assess read/contact evidence. | LEAD; No cut authorized | |
| CTLK-H02 | Not established | Not assessed at this position | Not assessed at this position | Not assessed at this position | Locate verified AGP component boundary and reconcile smaller/alternative alignments; then assess read/contact evidence. | LEAD; No cut authorized | |
| CTLK-H03 | Not established | Not assessed at this position | Not assessed at this position | Not assessed at this position | Locate verified AGP component boundary and reconcile smaller/alternative alignments; then assess read/contact evidence. | LEAD; No cut authorized | |
| CTLK-H04 | Not established | Not assessed at this position | Not assessed at this position | Not assessed at this position | Locate verified AGP component boundary and reconcile smaller/alternative alignments; then assess read/contact evidence. | LEAD; No cut authorized | |

## Definitions and provenance

- **Tracking scope:** As of 2026-10-05. Sample Sde-CTlk_104. Seven pipeline transitions and four explicitly identified historical-switch leads on four scaffolds. Entries are not eleven independent chimeric scaffolds.
- **Coordinates:** All coordinates are 0-based, half-open bp in the current pre-finishing assessment FASTA, not the final renamed/polished FASTA. Start is included; end is excluded.
- **Range meaning:** J entries: original chromosome-assignment uncertainty interval. H entries: bracket between selected historical alignment endpoints, not a verified gap or proposed cut. Separate gap columns record verified 100-N spans for J03/J04. Blank gap fields mean not established.
- **Chromosome names:** Current names are copied from the candidate tables. Historical names refer only to that historical assembly. Their numbers are not assumed homologous between runs.
- **Origin:** Present in raw contig establishes sequence history, not the biological or computational cause. YaHS joins connect components; 100 Ns do not measure the true physical separation.
- **Historical matches:** Exact extents are from minimap2 CIGARs. Historical and current assemblies may share reads, contigs and errors. Local haplotype correspondence varies.
- **HiFi:** Whole uncertainty intervals exceed retained read lengths. Zero full-interval bridges is uninformative. Gap-level zero support also requires local coverage, mapping ambiguity and physical-gap opportunity checks.
- **Anchors:** Reported placements are a limited search within one current haplotype. One reported placement does not certify uniqueness across both haplotypes.
- **Hi-C:** Local minima and ratios summarize the prior evidence report. They do not independently validate joins made with the same Hi-C. Library-specific alternative-adjacency tests remain outstanding.
- **Maintain this tracker:** Keep stable IDs. Add new evidence, update review status and decision rationale in amber columns. New rows may be appended to the Excel table. Retain original coordinates and hashes if a new assembly is compared; create a new version rather than silently changing the coordinate frame.
- **Scope limitations:** Historical-switch leads are not an exhaustive switch scan. Other historical fragments may warrant entries later. CMat hap1 reference remained unassessed in the candidate run; CPla unresolved chromosome inference is not clearance.
- **Source archives:** chimera-review-20261004-061953-1506189.tar.gz; chimera-origins-20261004-064841-1506290.tar.gz; chimera-history-20261004-083159-1506291.tar.gz. Source field identifies the relevant file families.

Assessment SHA256:

- hap1: 1fd4707ae4958ce31f7e95b5732e01d70686755c6a62a414b080fac24e4fd96a
- hap2: 121ce840edf743fb0a1c9bcc39ef210c241a585bd9e2f6dac3758176bc12ff19
