> **Presentation example:** detailed CTlk measurements are retained evidence; the full-report cohort overview contains example values.

# Chimera Detection and Review

## Assessment summary

2 assemblies assessed; 8 candidate boundaries in 2 assemblies. Evidence generation applies no cuts. All review selections start at NO.

[Editable review TSV](chimera_review.tsv) · [All assessment identities](assembly-registry.tsv) · [Cut instructions](cut-interface.md)

| Assembly | Candidate boundaries | Evidence |
| --- | ---: | --- |
| Sde-CTlk_104_hap1 | 6 | [Report](Sde-CTlk_104_hap1.review/report.md) |
| Sde-CTlk_104_hap2 | 2 | [Report](Sde-CTlk_104_hap2.review/report.md) |

## Candidate boundaries

| Candidate | Assembly | Scaffold | Chromosomes left → right | Region to review, bp | Exact cut, bp | Review priority | For cutting | Against cutting |
| --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| [C01](Sde-CTlk_104_hap1.review/report.md#candidate-c01) | Sde-CTlk_104_hap1 | scaffold_1 | chr9 → chr7 | 62618168–62941917 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CMat_203: chr9 → chr7 | No opposing chromosome evidence observed |
| [C02](Sde-CTlk_104_hap1.review/report.md#candidate-c02) | Sde-CTlk_104_hap1 | scaffold_1 | chr9 → chr7 | 63051225–63051325 | 63051325 | Prioritize gap-cut review | Separate chromosomes: Sde-CBau_104: chr9 → chr7; Sde-CMat_203: chr9 → chr7 | No opposing chromosome evidence observed |
| [C03](Sde-CTlk_104_hap1.review/report.md#candidate-c03) | Sde-CTlk_104_hap1 | scaffold_1 | Unresolved from qualified local alignments | 95877814–96047526 | Not assigned | Insufficient evidence to propose a break | No informative independent chromosome evidence supporting a break | No opposing chromosome evidence observed |
| [C04](Sde-CTlk_104_hap1.review/report.md#candidate-c04) | Sde-CTlk_104_hap1 | scaffold_1 | chr12 → chr12 | 116057447–116057547 | 116057547 | Evidence favors retaining this sampled boundary | No informative independent chromosome evidence supporting a break | Same chromosome: Sde-CBau_104: chr12 → chr12; Sde-CMat_203: chr12 → chr12 |
| [C05](Sde-CTlk_104_hap1.review/report.md#candidate-c05) | Sde-CTlk_104_hap1 | scaffold_1 | chr12 → chr12 | 116171488–116171588 | 116171588 | Evidence favors retaining this sampled boundary | No informative independent chromosome evidence supporting a break | Same chromosome: Sde-CBau_104: chr12 → chr12; Sde-CLim_110: chr12 → chr12; Sde-CMat_203: chr12 → chr12 |
| [C06](Sde-CTlk_104_hap1.review/report.md#candidate-c06) | Sde-CTlk_104_hap1 | scaffold_5 | chr4 → chr14 | 35443804–35581476 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CLim_110: chr4 → chr14; Sde-CMat_203: chr4 → chr14 | No opposing chromosome evidence observed |
| [C01](Sde-CTlk_104_hap2.review/report.md#candidate-c01) | Sde-CTlk_104_hap2 | scaffold_5 | Unresolved from qualified local alignments | 43055156–43233354 | Not assigned | Insufficient evidence to propose a break | No informative independent chromosome evidence supporting a break | No opposing chromosome evidence observed |
| [C02](Sde-CTlk_104_hap2.review/report.md#candidate-c02) | Sde-CTlk_104_hap2 | scaffold_7 | chr14 → chr4 | 39549503–39688365 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CBau_104: chr14 → chr4; Sde-CLim_110: chr14 → chr4; Sde-CMat_203: chr14 → chr4 | No opposing chromosome evidence observed |

IDs C01, C02, etc. are scoped to an assembly. Separate-chromosome assignments support reviewing a fusion; same-chromosome assignments oppose a fusion at that sampled boundary. No informative assignment means the measurement cannot decide; discordant within-individual assignments require investigation. Prioritize gap-cut review means a verified gap and at least two informative independent individuals with separate-chromosome assignments, without opposing or discordant chromosome assignments. This is a review aid, not a calibrated cut classifier.

No detected candidate is not a guarantee of structural correctness. Review ranges are measured chromosome-transition intervals, not permission to cut at their midpoint. An exact coordinate is needed only when requesting a cut. No opposing evidence observed is not equivalent to positive evidence for cutting. Detailed reports distinguish independent chromosome context, local support and measurement limits.
