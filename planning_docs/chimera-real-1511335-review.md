# CTlk automatic calibration review: 1511335

Exit status 0, all 21 tasks succeeded. Both CTlk assemblies remained unchanged and passed exact sequence verification. H01 was discovered and emitted as a review-only source-bound proposal.

Continuous-sequence controls now qualify: five to nine per original candidate, except hap1 scaffold5 with four. H01 has six qualified continuous controls and five independent-peer-qualified gap controls. However, neither library has any controls matching H01 within-flank coverage within a factor of two on BOTH sides. Candidate Hi-C flank counts are 157/4819 (Ex2) and 163/5565 (Ex3). Immediate HiFi flank measurements are zero/four qualified molecules with median depths 2/9. Thus support absence is confounded by poor left-flank mappability; increasing the same control pool does not establish automatic eligibility.

Hap1 candidate2 has seven contact-matched continuous controls in each library; hap2 candidate1 has nine Ex2 and six Ex3. Ex3 candidate/control ratios show strong contact depletion, but Ex2 does not reach the tenfold support-loss threshold. Neither broad interval is a uniquely localized gap and neither has usable independent chromosome-discordance votes. These are not supported automatic cuts. A weak library alone cannot override contrary/stronger support in another library.

The full automated result remains eight unresolved hypotheses, zero actions, and no declaration of confirmed biological fusion. Previously reviewed manual H01 correction and post-cut reassignment already passed separately.

The next work should use retained measurements locally to inspect farther-flank mappability/contact profiles and existing localized read evidence, report which eligibility criteria fail, and reconcile interval-level output with previously reviewed narrow seam assays. Avoid another identical whole-workflow rerun or threshold relaxation just to obtain a cut. The implemented joint requirement for both control populations must be interpreted as conservative calibration for literal-gap unjoins, not as a universal assay of every internal transition; non-gap intervals remain ineligible for automatic internal cutting independently.
