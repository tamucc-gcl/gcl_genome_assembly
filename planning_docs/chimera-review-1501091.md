# Scope-fix checkpoint: run 1501091

Static review of chimera-comparison-20260929-171356-reports.tar.gz.

- mighty_cantor completed SUCCESS: 55 succeeded, 234 cached, zero failures/retries.
- Correct input paths were included: assembly_samplesheet_original.csv and hic_readsets.none.csv. Five accepted individuals each had one Hi-C readset; ten assemblies finalized; CTlk_101 remained skipped.
- Submitted tasks were harmonization, chimera processing, finalization, eligibility and summaries. Upstream assembly/scaffolding reused cache.
- Reference remains Sde-CMat_203_hap1.
- CPla hap1/hap2 retain 73/51 compact composite candidate records, all NOT_A_CANDIDATE with unresolved chromosome scope. Their adjudicated join tables contain no calls and they produce no primary evidence tables.
- Exactly one primary evidence table remains: CTlk hap2 scaffold_3. The adjudicated result is REVIEW, callable=no, interval 33,660,739..33,830,451, diagnostic position 33,745,595, no compatible gap. Evidence cut_bp=NA, tidk source, figure generated.
- The CTlk hap2 assessment checksum remains 2cbed25444ba9e50e49aea7751d6068a22da49c89602413f4c5976ecf21acf23, identical to the earlier checkpoints.
- No adjudicated callable BREAK_CANDIDATE rows remain. The upstream scaffold-level CTlk BREAK_CANDIDATE is correctly downgraded during location assessment.

The scope correction passes this dataset checkpoint: evidence tables fell from 239 to 1 and the 10 CPla callable BREAK_CANDIDATE rows disappeared. This validates the observed negative/scope behavior, not positive automatic-cut sensitivity or a safe breakpoint for CTlk. An auto run on these calls would not resolve its remaining interval. Review that interval and the remaining concordance policy before calling the entire automatic-cut system validated. No further identical scope-only rerun is needed.

Archive extraction must include the baseline and current directories together because tar stores cached hard-linked files as links between those trees. The archive is usable; selective extraction of only current AGPs omits their baseline link targets.