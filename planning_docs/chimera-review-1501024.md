# Review of run 1501024

Static inspection of both 20260929-153558 and 20260929-153610 report archives. No pipeline or analysis program executed.

## Execution
- Run condescending_sinoussi completed SUCCESS: 55 succeeded, 234 cached, zero failures/retries.
- Invocation used assembly_samplesheet_original.csv and hic_readsets.none.csv; every accepted sample had one Hi-C library/readset. Five individuals, ten finalized assemblies. CTlk_101 skipped as expected.
- chimera_break=false, pangenome=false, post-assembly=false. No automatic cutting in this run.
- Submitted work was harmonization, chimera processing, finalization, eligibility and summaries. Assemblers and scaffolding did not appear among submitted tasks.
- Archives compare current output against the preserved immediate baseline and against tst/genome_assembly_store, respectively.

## Location repair succeeded on the motivating case
CTlk hap2 scaffold_3 has one REVIEW interval: 33,660,739..33,830,451 in pre-finishing coordinates. Position 33,745,595 is diagnostic only, cut_bp=NA in the evidence report. No compatible gap was found. Native and recovered tables retain this non-callable result. The old 30,202,188 proposal is not selected.

Assessment FASTA SHA256 is identical to the immediate baseline: 2cbed25444ba9e50e49aea7751d6068a22da49c89602413f4c5976ecf21acf23. This is a detector change on the same pre-finishing sequence, not a reassembly difference.

CTlk hap1 has zero candidate transition profiles. This does not independently establish whether the historical hap1 fusion was biological or how it disappeared; historical/current sequences differ.

## Scope regression: not ready for auto
CPla hap1/hap2 are passengers with no_sharp_dropoff and threshold_fallback chromosome-set methods. The fallback selects 283/255 scaffolds. My membership export incorrectly treated membership in that fallback set as confirmed chromosome-scale eligibility.

Evidence table counts:
- CPla hap1: 145
- CPla hap2: 93
- CTlk hap2: 1
- Other assemblies: 0

The final adjudicated CPla tables include 8 and 2 callable BREAK_CANDIDATE rows, respectively. These are recommendations that could reach automatic cutting, subject to cutter validation/minimum-piece rules; they are not observed cuts. Cutting was disabled.

The fix must distinguish a usable naming fallback from confidently inferred chromosome membership. A failed dropoff inference should remain explicitly unresolved for primary chimera profiles and cuts. Do not replace this with another absolute size floor or exclusion by sample ID. Do not use voter/passenger status as a substitute for chromosome inference quality. Preserve compact audits for excluded/unresolved scaffolds.

Add integration coverage for threshold_fallback, successful dropoff, and successful chromosome inference in a non-voter. The current tests checked the membership flag downstream but missed its incorrect upstream derivation.

## Other observations
All 239 evidence tables report tidk as their telomere source. The detailed REPORTING workflow was skipped because QC and post-assembly analysis were off; the assembly_run_summary and individual evidence outputs are available. This checkpoint therefore does not validate the detailed report renderer.

The collector archived default input filenames, not the actual alternative filenames used by this invocation. The launch log and accepted-input summary verify the chosen paths and library counts, but future collections should include the exact invoked sheets.

## Next checkpoint
Repair inference-quality propagation and its tests, then repeat detection-only with unchanged single-library inputs and retained work/cache. Keep auto off and do not restore the second library yet. No reason to repeat hifiasm or scaffolding to test this correction. The longer-term per-junction, individual-aware concordance policy remains separate from this scope fix.