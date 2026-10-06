# Finalizing chromosome-scale chimera adjudication

Date: 6 October 2026. Updated after user clarification and restoration of comparisons. Status: proposed investigation and implementation plan; no cuts authorized or executed.

## Recommendation and limits of this review

Finish a bounded adjudication of chromosome-scale suspect structures before expanding optional analyses. Use cross-assembly chromosome-block correspondence to discover candidates, exact boundary evidence to adjudicate them, and a separate validated action manifest to apply corrections. Add an origin-directed prevention experiment alongside adjudication. Do not make a reference-label switch, zero HiFi bridges, or haplotype vote an automatic cut rule.

The working assumption is that an upstream algorithm had evidence or a heuristic basis for each connection. Reconstruct and evaluate that basis; it is not a presumption that the connection is correct. Consistent separation into substantial chromosomes in six better assemblies is positive comparative evidence against the CTlk composite representation, subject to verifying homology and block order. Do not require biological impossibility before removing an unreliable scaffold connection. A suspect structure remains suspect until explicitly adjudicated; presence of some contacts/alignments does not by itself clear it.

The operational target is a defensible assembly action, not proof of evolutionary history. A poorly supported scaffold connection can be removed as UNJOIN_UNSUPPORTED even when a biological adjacency remains possible. An unresolved internal-contig transition can remain intact with a reported limitation. Neither outcome should trigger an indefinite investigation.

Reviewed initially: the current handoff, junction assessment/focus/history reports, adjudication checkpoint, CTlk registry, and current discovery/review/evidence/cutting source. Comparisons was initially inaccessible but is now restored. Follow-up inspection includes archived name maps/candidate tables, foundation transition tables/configuration, origin ledgers/provenance, assessment gap links, focused status, and current mapping/scaffolding modules. This verifies specific recorded findings below, not every underlying alignment or molecule. No biological analysis/test runs were executed. The unrelated existing .gitignore modification is outside this plan.

## Updated CTlk scope and findings

User confirms job 1511195 was cancelled and did not complete. Inspect and salvage its surviving outputs without another status-confirmation request; no local archive for that job was found in the current top-level comparisons listing.

The restored maps confirm:

| Haplotype/scaffold | Length | Assigned composition | Required associated pieces |
|---|---:|---|---|
| hap1 scaffold_1 | 138,446,979 bp | chr7_1 + chr9_1 + chr12_1 | scaffold_11, chr7_2, 38,363,695 bp; all other chr7/9/12 sequence |
| hap1 scaffold_5 | 74,765,583 bp | chr4_1 + chr14_1 | scaffold_13, chr4_2, 35,299,899 bp; all other chr4/14 sequence |
| hap2 scaffold_5 | 83,430,767 bp | chr7_1 + chr12_1 | scaffold_14, chr7_2, 41,302,352 bp; chr9 scaffold_9 and scaffold_17 |
| hap2 scaffold_7 | 82,879,935 bp | chr4_1 + chr14_1 | scaffold_15, chr4_3, 35,652,342 bp, and scaffold_18, chr4_2, 5,742,418 bp |

Thus 4+14 is shared at chromosome-pair level, and 7+12 is also shared; the additional chr9 incorporation distinguishes hap1 scaffold_1. Both haplotypes have separate large chr7-labelled pieces. Names describe inferred correspondence, not proven whole chromosomes or proven breakpoints. The question is whether the same *subregions and oriented boundaries* recur, not whether the same two labels appear anywhere on a scaffold. The existing 2f/6s pair vote cannot answer that question.

User clarification: chr4_2 occurs in both haplotypes; chr4_3 occurs only in hap2. Whether the numbered pieces match across haplotypes is unknown. Include every such piece in the content/copy audit, identifying correspondence by aligned sequence rather than part suffix. Determine whether all pieces cover complementary chromosome intervals, duplicate the same intervals, or include phase/copy mixtures. In particular, matching 4+14 labels in both haplotypes does not clear the fusion if the rest of chr4 is partitioned differently.

The six better comparator assemblies are the two haplotypes each of CBau, CLim, and CMat in the restored cohort; CPla is highly fragmented and should provide local sequence evidence rather than chromosome-order votes. These are three other individuals, not six independent biological replicates. The current CTlk name maps show broad overlapping reference spans for separately named pieces. Determine complementary versus duplicated coverage from base-resolved alignment unions; `_1`, `_2`, and endpoint ranges alone cannot distinguish them.

**A crucial origin distinction:** the apparent major reference transitions are already present within raw hifiasm contigs:

- hap1 scaffold_5, chr4 -> chr14 uncertainty interval 35,443,804–35,581,476: the origin ledger maps the surrounding window exactly to h1tg000153l, 619,979–957,651.
- hap2 scaffold_7, chr14 -> chr4 interval 39,549,503–39,688,365: surrounding window maps exactly to h2tg000298l, 538,671–877,533, reverse.
- hap2 scaffold_5, chr12 -> chr7 interval 43,055,156–43,233,354: surrounding window maps exactly to h2tg000028l, 560,967–939,165, reverse.
- hap1 scaffold_1's major chr9 -> chr7 and chr7 -> chr12 label transitions likewise map to raw contigs h1tg000004l and h1tg000133l.

These are uncertainty windows, not executable cuts. Their survival from raw contigs does not validate chromosome identity or physical correctness. It establishes that changing YaHS joining parameters alone cannot be assumed to eliminate these exact sequence transitions. Conversely, H01 and other verified round1 gaps may attach large blocks incorrectly even when the nearest reference-label transition lies inside a raw contig. The correspondence transition and the erroneous physical adjacency may differ. We must compare block layouts across peers to locate the actual problem rather than moving the cut to the nearest gap.

**Terminology correction after inspecting the supplied old name maps:** continuity across a current reference-assignment transition is not evidence that the old assembly contained the current chromosome-scale fusion. The old hap1 map records scaffold_1 as two subintervals, 0–37,144,322 and 37,144,322–111,642,300, named chr5_1 and chr9_1. Old hap2 records scaffold_3 split at 43,694,678 into chr12_1 and chr6_3. This is evidence that those scaffold objects were divided before final naming, although the map alone does not establish when, why, or whether the original connections came from hifiasm or YaHS. A name map describes final records and does not enumerate every adjacency within a contig or earlier scaffold. Old hap2 also retains other composite names, chr2_1+chr3_2 and chr7_3+chr11_2; neither is the current CTlk composite identity.

For example, the current hap1 chr4/14-labelled transition window aligns through the uncertainty interval inside old chr15_2. That demonstrates local sequence continuity under the historical alignment, not an old chromosome explicitly named chr4+chr14. Different naming frames, incomplete or incorrect assignment, and wrong connections elsewhere on a large scaffold can explain this distinction. The causal investigation must locate *actual structural adjacencies*, rather than treating a changing reference label inside continuous sequence as a demonstrated chromosome-fusion breakpoint. Historical recurrence alone does not clear the current composite or establish its cause.

The current disposition table below assesses individual boundaries. Retaining H02's final gap or finding reads across J05 does **not** clear the surrounding composite scaffold. The scaffold-level concern stays open until its block structure has been explained.

Sources: assignment-review-1505984 harmonized name maps/candidate table; local-review-1506189 foundation chimeric_joins.tsv; local-review-1506290 origin_ledger.tsv/provenance.json; local-review-1506302 selected_gaps.tsv/seed_gap_links.tsv. Historical chromosome names belong to a different frame and must not be equated directly to current chromosome numbers.

## 1. Establish the reusable evidence baseline

Before launching biological work, obtain the surviving outputs/accounting of the cancelled, incomplete 1511195. Preserve its directory. Inspect output validity per assembly and stage; a log printed before mapping is not completion evidence. Do not assume cache reuse.

Build a small evidence manifest for all ten assessment assemblies: individual, haplotype, FASTA digest, coordinate stage, AGP digest and lineage, read mappings and their exact source FASTA, library identities, available assembly comparisons, and missing artifacts. Import completed sister mappings and earlier packets only after identity, parameters, and coordinate checks. Final-coordinate BAMs cannot substitute for pre-finishing mappings without an explicit verified translation.

**Exit product:** one inventory identifying reusable, invalid, missing, and unfinished evidence. No new assembly or whole-genome HiFi remapping is needed for this step.

## 2. Discover the complete, bounded candidate set

Reconcile production candidates, historical leads, and the CTlk investigation registry. The registry is not batch-wide discovery. Current production discovery in chimera_joins.py depends on chromosome-member composite names and reference assignments; this can miss a fusion in the reference itself or a connection whose name does not expose the conflict.

Use existing assembly alignments first to build a reciprocal block correspondence map across sisters and suitable peers. Assign chromosome-group identities by distributed, distinctive sequence correspondence, retaining substantial competing placements. Names and hap1/hap2 labels are not identity evidence. Measure aligned bases by interval unions, coverage of the corresponding peer chromosome/group, order, orientation, and unassigned or duplicated sequence.

Candidate triggers should include:

- A scaffold containing two or more substantial chromosome groups that are separate in informative sister/peer assemblies.
- A sister pair sharing one chromosome group but giving it different large partners, such as A+B versus A+C.
- A major correspondence transition within a contig, or near a verified component boundary, even without a composite scaffold name.
- Existing flagged historical transitions and complex out-and-back arrangements requiring reconciliation.

Separate whole/near-whole chromosome combinations, partial-arm combinations, local insertions/reorderings, and assignment ambiguity. Define "substantial" relative to empirically inferred chromosome-group sizes and distinctive coverage; do not use a universal fish size threshold. The reported 15 chromosomes and genuine fusion history are context, not a forced scaffold count.

Audit the selected graph reference with the same method and repeat the candidate screen with that reference excluded or replaced. Cross-assembly consensus can share an error; identify shared reads, contigs, libraries, and reference influence rather than counting correlated assemblies as independent votes.

Cheap scaffold/AGP inventory can be complete, but expensive boundary analyses must receive an explicit selected list: suspect regions plus a fixed, matched control panel. Within a complex region include its relevant component boundaries; do not analyze every gap on every scaffold. Record every excluded lead and why it was excluded. Merge overlapping registry/gap records into one physical-boundary identity with multiple discovery reasons.

**Exit product:** batch-wide candidate table, chromosome-block maps, suspect boundary groups, and fixed controls. If missing current alignments require mapping, map only the selected informative blocks to sisters and quality-selected peers before considering whole-assembly mapping.

## 3. Test competing arrangements, with a declared decision purpose

For each candidate specify the current adjacency, the split hypothesis, and any sequence-supported alternative order/orientation/partner. Include phase-switch or repeat/copy ambiguity when relevant. A+B versus A+C requires review of both boundaries and the full diploid content, not selection of a winner by majority.

### Test A: reciprocal chromosome and component correspondence

Compare multiple distinctive anchors along both large blocks, not just two nearby flanks. For each boundary, retain several anchor distances and explicitly mark whether an anchor crosses another join. Confirm reciprocal coverage, orientation, copy identity, and whether the sister really contains the relevant sequence.

Different target scaffolds support contradiction only if they are informative chromosome groups; sister fragmentation, scaffold ends, missing sequence, or uncertain phase must be distinguished from a confidently incompatible structure. Ordered outer flanks establish regional correspondence, not the intervening component's exact adjacency.

**Resolves:** chromosome-scale fusion versus local anomaly/assignment artifact; credible comparative contradiction versus missing comparative information.

### Test B: boundary Hi-C structure and alternatives

Reuse retained deduplicated contacts with verified coordinate projection. Produce a whole-suspect-scaffold map and targeted boundary maps, separately by library. Examine continuity of distance decay, cross-boundary enrichment/depletion, and whether the two blocks behave like separate contact domains. Show competing partner regions with the same plotting and counting rules.

Keep the existing 25/100/250/500 kb scales for comparability, adding a chromosome-block view where needed. Mark every other join included in a window. Normalize interpretation for usable sequence, flank contact opportunity, repeat/mappability, coverage, and position; raw counts and per-million totals alone are insufficient. Inspect raw counts alongside normalization so poor coverage is not hidden. Low-information flanks must return uninformative, not failed.

Fit empirical expectations from supported intra-chromosomal regions matched for opportunity and location, and between established chromosome groups. Select controls independently of the score being evaluated where practical. Nearby gaps are comparison joins, not truth. Nested windows and Ex2/Ex3 are not independent votes. Hi-C already used for scaffolding is consistency evidence; do not label it independent validation. A retrospective read split does not create a holdout if all reads influenced the existing assembly.

**Resolves:** coherent chromosome continuity versus a supported boundary discontinuity; whether alternative partners improve the explanation; whether weak counts are interpretable at all.

### Test C: targeted molecule continuity and phase

First determine whether available HiFi lengths and local coverage could span the nearest distinctive anchors. Scaffold placeholder Ns do not estimate physical separation. If opportunity cannot be estimated, absence cannot become negative evidence.

Inspect existing retained reads and CIGARs first. Where a decision depends on ambiguous placements, extract the relevant reads and compare them competitively against both haplotypes and relevant repeat/alternative copies, retaining alternatives. Distinguish unrelated repeat copies from homologous alleles. Count independent molecules with distinctive anchors on both sides, sufficient aligned sequence, acceptable error/indel structure, and compatible local phase where informative variants exist. Absence of an emitted secondary alignment is not proof of uniqueness.

For internal-contig transitions localize a candidate breakpoint using read/path discontinuity; then test reads across that exact site. A tiled support track does not establish one coherent phased path across a long repeat. Graph evidence requires actual connected paths/overlaps and their molecule support, not isolated S/A records. Attempt graph reconstruction only if the retained graph can answer the named unresolved question.

**Resolves:** direct support for the claimed boundary, a localized unsupported discontinuity, or inability of current molecules to distinguish alternatives.

### Test D: diploid consistency and phase errors

Summarize both homologues' chromosome-group content, copies, partners, and local phase transitions. Test whether apparent A+B versus A+C can be explained by fragmented representations, swapped homologues, copy ambiguity, or phase switching. Investigate reciprocal products where expected under the proposed model; do not require them for every biological mechanism.

Treat an apparent haplotype-specific whole-chromosome fusion as an exceptional claim requiring stronger affirmative boundary and phase evidence. No numeric biological prior is established for this species. Disagreement is a scrutiny trigger, not automatic proof of misassembly. Conversely, shared fusion structure in sisters can reflect shared assembly error and still needs boundary support.

**Resolves:** whether the assembly difference supports a credible diploid arrangement or remains an assembly/phase ambiguity. Without decisive phase information, retain biological uncertainty even if unjoining is the defensible assembly action.

## 4. Current CTlk priorities and terminal outcomes

These are starting dispositions from reported evidence, not accepted cuts. The table is not the complete batch candidate list. All locations refer to assessment assemblies; gap starts are zero-based.

| Case | Existing evidence | One bounded next question | Terminal rule |
|---|---|---|---|
| H01, hap1 scaffold_1, 63,051,225 | Weak multiscale contacts; informative displaced anchors separate in history and one current peer; depleted immediate left flank | Do block-scale contacts and reciprocal sister/peer correspondence support separate chromosome groups after accounting for opportunity? | BREAK_PROBABLE_MISJOIN if informative contradiction and exact gap are established; UNJOIN_UNSUPPORTED if the connection remains unreliable without positive validation; otherwise RETAIN or UNRESOLVED according to evidence |
| H02 final gap, 95,221,756 | Ordered correspondence in current sister and old hap1; contact recovery at wider scales | Is there a specific boundary contradiction beyond historical separation? | Default to retaining with documented support unless a new contradictory observation changes the case; no rerun merely to reconfirm history |
| H02 internal gap, 94,575,157 | Depressed contacts and local distance/order discrepancy | What are the component copies/order at the clean immediate anchors? | Resolve as local structural issue; do not call a whole-chromosome fusion from this discrepancy alone |
| J03/J04, 116,057,447 and 116,171,488 | Two boundaries around a 113,941 bp component; outer continuity masks individual boundary support | Does this component belong here, in this orientation/copy, and which of its two boundaries is supported? | Evaluate jointly; preserve the component as sequence even if detached. Existing minimum-piece cutter may reject this action; use an explicitly supported component-preserving route |
| H03/H04 | Contacts broadly comparable to available comparison joins | Does current batch-wide sister/block analysis identify credible contradiction? | Retain provisionally unless the new screen supplies a specific contradiction; historical splits alone do not justify cuts |
| J01/J02/J06/J07 | Continuous screened local read support; historical recurrence | Does the block/phase screen reveal a specific discontinuity requiring localization? | No internal cut on current evidence; close as no cut justified, with any biological uncertainty recorded |
| J05 | Low-MAPQ pocket; one bracketing primary molecule; no demonstrated repeat competition | Can inspection of the retained spanning read and its alternatives change the action? | No cut justified currently; at most one targeted ambiguity test if decision-relevant; no midpoint cut |

## 5. Decision policy and automation

Record two separate fields: assembly action and biological interpretation. RETAIN does not certify an evolutionary fusion; UNJOIN_UNSUPPORTED does not disprove adjacency.

| Action | Required basis |
|---|---|
| RETAIN | Affirmative boundary support or coherent informative comparative/contact evidence sufficient for the assembly claim, with serious contradictions resolved. Haplotype-specific chromosome fusions require stronger direct/phase support and remain manual initially. |
| BREAK_PROBABLE_MISJOIN | Informative structural contradiction plus local evidence against the connection and a validated boundary. Comparative votes alone are insufficient. |
| UNJOIN_UNSUPPORTED | Verified scaffold connection lacking adequate affirmative support after bounded assessment; conservative loss of asserted adjacency is preferable to retaining an unreliable chromosome-scale connection. Manual initially. |
| UNRESOLVED | Evidence opportunity inadequate, competing structures remain indistinguishable, or safe internal cut cannot be localized. Specify the remaining limitation and its effect on release. |

Start with automatic candidate discovery, evidence packets, triage, and validated execution of reviewed decisions. Automatic unjoining/classification is a later gated capability: require exact gap eligibility, informative and consistent evidence, no credible direct support for the existing join, no unresolved phase/copy conflict, and calibrated performance on held-out tests. Until calibration passes, route those cases to manual review. Keep internal sequence cuts and exceptional biological-fusion claims manual. Do not invent an "extraordinarily high" count threshold now.

The decision table must include assembly digest/stage, boundary and component IDs, chromosome-block identities, anchors and ambiguity, comparative evidence and dependencies, contact opportunity/results by library, molecule/phase evidence, alternatives, action, rationale, reviewer/policy version, and a specific missing decisive observation or explicit statement that assessment is closed.

## 6. Validation that can justify automation

Use a bounded benchmark covering:

1. Supported native chromosome boundaries and genuine supported fusion structures where available; do not label nearby gaps as correct by assumption.
2. Artificial chromosome-scale wrong joins made from confidently separate blocks, including A+B in one haplotype and A+C in its sister. Include repeat-rich/poor-coverage endpoints and local insertions as separate classes. Transform read coordinates correctly or remap only relevant reads; injecting Ns into a FASTA does not by itself generate valid contact evidence.
3. Fragmented sisters, missing peer sequence, phase switches, homolog ambiguity, both-haplotype shared errors, a chimeric reference, and a true heterozygous rearrangement challenge case. If empirical true-fusion examples are unavailable, report that generalization limit; simulation is not biological validation.
4. Internal-contig errors with known breakpoints, including repeats longer than reads, to verify that unsafe cuts abstain.

Split calibration from evaluation by individual/chromosome group, avoiding nearby correlated windows in both sets. Select thresholds for low false-cut risk and report detection sensitivity, false-cut rate, abstention rate, breakpoint accuracy, and results per evidence-opportunity class. State sample-size uncertainty; zero false cuts in a small benchmark is not established safety. Confirm decisions are stable to sensible anchor/window/filter changes; unstable cases become manual/unresolved.

Required implementation checks include duplicated candidate identities, checksum/stage mismatches, coordinate reversal, ambiguous lineage, missing evidence, compound boundaries, sequence conservation, and retaining all fragments. Run tests and biological pilots on the user's cluster after source synchronization; this review did not execute them.

## 7. Bounded implementation and release sequence

### Prevention experiments: change the stage that created the problem

First reconstruct the effective command, installed version, exact input digests, libraries, and intermediate AGPs/graph path for each suspect adjacency. Configuration defaults are not proof of the parameters used by cached upstream tasks. Identify the earliest YaHS resolution/round at which a gap connection appeared or returned after correction. For an internal-contig transition inspect the native graph/path and hifiasm inconsistency output where retained. If exact internal support/scoring was not saved, report that limitation rather than claiming to reproduce the algorithm's rationale.

Current-source and archived-default audit identifies these concrete tests:

| Experiment | Why test it | Scope and interpretation |
|---|---|---|
| YaHS round2 scaffold error correction enabled | Current/archived defaults set yahs_round2_scaffold_ec=false; round1 has it enabled | One change using unchanged round2 inputs. May reject/rebreak weak joins, but cannot be assumed to remove within-contig transitions; verify actual stage commands |
| Effective mapping-quality contract | YaHS receives a coordinate-sorted filtered BAM; its `-q` is suppressed for that sort order | Upstream pairtools already defaults to hic_min_mapq=30. Audit retained alignments and actual run settings first. A tested pair-level filter or separate name-sorted YaHS input is needed for an effective threshold change; preserve indexed coordinate BAMs for downstream consumers. Sorting itself changes how YaHS places links, so isolate that effect |
| Conservative YaHS stopping resolution/iterations | Joining at coarse resolution may attach blocks whose local boundary evidence is poor | Identify the first offending intermediate AGP, then stop before that resolution in a diagnostic variant. Avoid an arbitrary genome-wide maximum; compare preservation of legitimate long-range joins |
| YaHS telomere-end protection | Current wrapper does not expose `--telo-motif`; terminal telomeres can supply an evidence-based prohibition on end joins | Check installed support and taxon motif. Test only where input ends actually have terminal telomere evidence. Internal motif occurrences cannot veto genuine fusion sequence; modern YaHS may already use a motif database |
| Hifiasm post-join disabled | Current/archived hifiasm_post_join=1 permits a contig post-join step; major suspect transitions occur in raw contigs | If graph/origin evidence implicates that step, run one CTlk variant with -u 0, using validated compatible saved correction/overlap files in an isolated experiment. Compare exact contig paths, phase and retained sequence. It may have no effect if the connection formed earlier |
| Hifiasm coverage/overlap or misjoin settings | --b-cov is 0 and --h-cov is -1, so coverage-triggered breaking is disabled; --l-msjoin is already enabled at 500 kb | Only test after local coverage/overlap evidence identifies a mechanism. Choose coverage settings from sample distributions and matched controls; no arbitrary threshold sweep. --lowQ creates an inconsistency annotation, not a cutting rule. --l-msjoin is a minimum eligible unitig size, not a generic confidence threshold |

Dual scaffolding is already disabled, and hifiasm already receives a taxon telomere motif in the wrapper; do not propose disabling/enabling those as if they were unexamined remedies. The recorded foundation launch already supplies the Arima four-enzyme motif list to both YaHS rounds; confirm it matches library chemistry rather than treating missing enzyme settings in defaults as the cause.

Run prevention variants on CTlk plus one better individual as a control, maintaining both haplotypes. Start with the baseline and the relevant single-stage/single-factor change; promote combinations only after a single-factor result warrants them. No full-cohort assembly restart and no large parameter grid. Compare suspect block adjacencies, valid chromosome continuity, copy/content retention, phase consistency, contact explanation, fragmentation, and resource cost. N50 or attaining exactly 15 scaffolds cannot be the acceptance objective. Prefer a setting only if it removes independently adjudicated weak connections without unacceptable loss of supported ones. If prevention resolves an unsupported gap, still record why the old join was rejected and audit the new structure.

Official references, to be checked against installed versions:

- YaHS input sort/MAPQ behavior, scaffold error correction, intermediate AGPs, telomere protection: https://github.com/c-zhou/yahs .
- Hifiasm post-join, inconsistency BED, coverage break and unitig misjoin options: https://github.com/chhylp123/hifiasm/blob/master/CommandLines.cpp .
- Hifiasm retained corrected reads/overlaps and self-scaffolding behavior: https://github.com/chhylp123/hifiasm .

### Recommended next bounded packet

#### Prioritize the one-library versus two-library causal comparison

Historical outputs are optional diagnostic context for this particular regression. They are not required production evidence, a prerequisite for general-purpose detection, or a reason to postpone current-cohort adjudication. The general algorithm must work from current raw graphs/reads, stage lineage, sisters, and current quality-selected peers. Prefer current evidence to resolve structure; use the old example only where it can cheaply isolate what changed. Controlled library ablation can use current inputs even if historical assemblies are unavailable.

The user reports that the chromosome-scale composites were absent with one Hi-C library and appeared with two, while HiFi is unchanged. This warrants direct investigation of hifiasm as well as YaHS. It does not yet establish that the local raw-contig connections were newly created: the historical comparison packet reports sequence continuity through all five major internal-transition windows in older assemblies, with exact CIGAR matches across the uncertainty intervals and anchors. Therefore compare old/new *adjacencies and paths*, not merely old/new chromosome names or final scaffold compositions. Local continuity can predate the added library while chromosome-scale scaffolding changes around it.

HiFi supplies corrected sequence/overlap structure. Hi-C contributes long-range phasing decisions that can change which paths are emitted from ambiguous structure. It is not a simple low-weight vote appended to a fixed contig FASTA. If there is a uniquely supported HiFi connection, changing Hi-C should not casually overturn it; if there is repeat/haplotype ambiguity, a new contact distribution can change partitioning, path selection, and subsequent contig outcomes. These are hypotheses for this case, not a demonstrated cause.

The wrapper supplies comma-separated `--h1` and `--h2` lists. The recorded new manifest contains Ex2 and Ex3. It provides no explicit per-library weights to hifiasm. Do not interpret `--n-weight` as a library-balance control: it sets rounds of link reweighting. Confirm which library/read set was used historically rather than assuming the single-library baseline is Ex2. Verify paired list order, same individual, preprocessing/HiFi digests, actual versions, coverage estimates, seeds, and effective commands. Library quality, depth and repeat-specific contact differences are candidate explanations; none has been measured here.

**First use retained old/new outputs:** obtain full `r_utg.gfa`, `p_utg.gfa`, relevant phased/primary contig graphs, inconsistency BEDs, hifiasm logs and compatible corrected-read/overlap caches. Map suspect windows and distinctive neighboring sequence into both old/new unitig graphs. Reconstruct the actual oriented path and each supporting overlap across the connection; identify competing paths and local phase/copy ambiguity. The previously extracted S/A contig snippets cannot supply this reconstruction. A contig S record is already an assembled result, and its A records alone do not prove the correctness of its internal path. If intermediate path evidence is absent, record it and regenerate only the necessary assembly stage with supported diagnostic outputs.

For each implicated path record: present in old/new processed unitigs; present in old/new final contigs; selected neighbor/orientation; distinctive versus repeated anchors; read-overlap counts/error/coverage; phase markers where available; Hi-C link support for selected and alternative branches separately by library where recoverable; post-join involvement; and the first downstream scaffold connection that changes its chromosome-scale composition. Lack of exported internal scores must be explicit.

**Then, if the comparison does not resolve causality, run a controlled CTlk ablation:**

| Treatment | Held fixed | Question |
|---|---|---|
| Historical single library only (A) | Current pinned hifiasm version, identical HiFi/overlap inputs and settings | Can the one-library result be reproduced without historical software/parameter confounding? |
| Added library only (B) | Same | Does B alone favor the suspect path or assignment? |
| A+B | Same | Does combining libraries alter that path relative to either alone? |
| HiFi-only | Same compatible graph-generation settings | Is the connection already present without Hi-C? Treat haplotype/phase differences as expected; compare local sequence/path structure rather than hap1 labels |

Inspect hifiasm outputs before running YaHS. Reuse the existing valid A+B treatment where possible. Share only validated compatible HiFi correction/overlap intermediates in isolated experimental directories. Hi-C caches must match the exact treatment: do not copy stale `*hic*.bin` across A/B/A+B, rely on unchanged graph identity to validate changed Hi-C reads, or delete production caches. A changed read list requires newly valid Hi-C mappings. If A and A+B show a reproducible difference, add one depth-matched paired subsampling comparison to separate extra contact depth from library composition. Seed sensitivity is a targeted follow-up only if the branch decision appears unstable, not a large repeat matrix.

If unitig structure is unchanged but final haplotype paths change, prioritize phasing/purge/path-selection evidence. If the same internal path exists in every treatment, prioritize incorrect reference assignment or a shared HiFi ambiguity, and identify the new YaHS connections creating the larger composites. If graph-generation structure changes unexpectedly, audit hidden input/version/settings/cache differences before blaming contact weighting. HiFi-only contiguity does not establish phase correctness or biological validity.

**Settings follow-up:** disable post-joining in one discriminating treatment if the implicated path forms during that step. The installed binary's help/source is authoritative: older docs describe `-u` as a switch, while current source accepts 0/1 and this wrapper supplies `-u 1`. Confirm syntax and behavior before using a `-u 0` variant. Examine inferred homozygous coverage and duplicate purging if copies/phase are implicated. Do not change error-correction/overlap parameters first simply because a final haplotype path changed; that would destroy the controlled comparison. `--l-msjoin` is already enabled and its size cutoff is not a general Hi-C trust setting.

Official Hi-C integration/cache behavior and parameter distinctions: https://github.com/chhylp123/hifiasm/blob/master/docs/source/hic-assembly.rst and https://github.com/chhylp123/hifiasm/blob/master/docs/source/faq.rst . Verify details against the installed release; online documentation and current source have version-dependent differences.

For the four CTlk composites and their associated chr4/7/9/12/14 pieces, reuse the existing reference/peer PAFs to build an ordered, copy-aware chromosome-block correspondence map against the six better assemblies. Add only missing targeted comparisons. Overlay raw-contig identity, all informative physical boundaries, the existing H/J tests, and terminal/phase evidence. This packet answers whether the problem is misleading assignment, raw-contig continuity, wrong scaffolding of real blocks, or a mixture; it selects the precise next boundary tests and the relevant prevention experiment. It should not repeat a full contact scan or HiFi remap before selecting those tests.

The immediate deliverable is a four-scaffold hypothesis table with supporting and contradictory observations, exact candidate boundaries or localization limits, the original assembly rationale where recoverable, and a named decisive test for each remaining ambiguity. H01 remains a priority, but neither H01 alone nor the current CTlk registry exhausts the four composite-scaffold problems.

**First change: scope, provenance, and runtime.** Add explicit selected-candidate input distinct from inventory. Save progress/checkpoints by assembly and evidence stage, with validated cache keys including source digests, tool/options, and relevant candidate/anchor definitions. Index contact partner counts by library/junction/side instead of rescanning the entire partner dictionary for each side. Stream contacts once per assembly for all selected windows/scales. Log pair counts, elapsed stage times, cache hits, and evidence availability. Profile a small pilot before assigning batch resources; the prior bottleneck is suspected, not measured.

**Second change: comparison and decision packets.** Produce the batch-wide block screen, bounded candidate/control list, and the four named tests. Render compact chromosome/block and boundary evidence panels with ggplot2/tidyverse and patchwork, separate from analysis. Each test declares which decision it can change. Run a pilot including H01, a comparatively supported join, and an ambiguous/internal case; then process remaining selected candidates only if the pilot validates scope and outputs.

**Third change: action export.** Current junction_review decisions and break_chimeras input contracts differ. Add a validated exporter that selects only BREAK_PROBABLE_MISJOIN/UNJOIN_UNSUPPORTED rows, binds them to the current FASTA/AGP, and refuses unresolved/noneligible boundaries. The cutter currently rejects auto mode, treats supplied rows as instructions, requires chromosome-member validated gaps, and defaults to 1 Mb minimum pieces. Do not pass the whole decision table directly. Add a separate reviewed, sequence-preserving component-detachment route if J03/J04 need small fragments; never evade the guard by fabricating gap eligibility or silently dropping fragments.

**Fourth change: correction and core acceptance.** Dry-run each action to show exact boundaries, output fragments, sequence conservation, and naming/provenance changes. Apply accepted actions before finishing, preserve rejected component adjacencies so later scaffolding cannot silently recreate them, and audit lineage after any downstream rerun. Verify intended joins are absent, untouched regions are conserved, all sequence/copies remain represented, coordinate maps remain valid, and graph filtering retains corrected fragments. Choose the graph reference after applying the same review to it. Rerun only affected descendants and complete core QC/pangenome/report acceptance before optional analyses.

Budget the biological investigation as one evidence-reuse audit, one small pilot, one bounded batch, and one targeted follow-up per unresolved decision group. The budget limits the investigation, not the need to fix software failures. After the follow-up, assign a terminal action/limitation; new testing requires a named observation likely to change it. Do not continually add evidence classes.

**Core release gate:** every candidate has a closed action and rationale; all cuts have validated coordinates and conservation audits; residual chromosome-scale uncertainties are prominently reported and excluded from unsupported adjacency-dependent claims; no silent fragment loss or recreated rejected joins; current-run QC, graph outputs, and consolidated report pass their existing acceptance checklist. Automatic biological certainty is not a prerequisite for a working core pipeline.

## Scientific context

The proposed cross-assembly evidence layer supplements, rather than repeats, native assembly/scaffolding support. YaHS already uses Hi-C for scaffold construction and optional error correction, so reuse of those contacts must be interpreted accordingly: https://pmc.ncbi.nlm.nih.gov/articles/PMC9848053/ . Hifiasm's phased-graph approach explains why local sequence continuity and haplotype correctness need separate assessment: https://www.nature.com/articles/s41592-020-01056-5 . Neither source establishes a species-specific probability for heterozygous fusion; the stricter review policy above is an assembly decision policy, not a measured biological prior.
