# CTlk parallel decision batch

This batch replaces the sequential regeneration-first plan. The user establishes that HiFi and all hifiasm settings are unchanged except for addition of the second Hi-C library. That is accepted; no historical parameter audit is a prerequisite. No historical assembly is required by these tests.

## Launch both arrays after source sync

From `/work/birdlab/GCL/spratelloides_delicatulus_genome`:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_hifiasm_fusion_batch.sbatch \
  HIFIASM_ENV genome_assembly comparisons/chimera-foundation-01 ANALYSIS_ENV

sbatch gcl_genome_assembly/scripts/comparisons/run_yahs_fusion_batch.sbatch \
  YAHS_ENV comparisons/chimera-foundation-01
```

Replace environment names with existing working environments. HIFIASM_ENV must contain the original hifiasm 0.25.0-r726 and gfatools; ANALYSIS_ENV needs Python, samtools and minimap2. They can be the same environment. YAHS_ENV needs YaHS and samtools. No package installation or pipeline restart is requested. Logs directory is the existing project logs directory. Input paths follow the current published layout; missing inputs fail promptly instead of being rebuilt silently.

Hifiasm runs up to three exclusive-node tasks concurrently (48 CPUs, 96G each, 8h limit). YaHS runs up to three exclusive-node tasks concurrently (8 CPUs, 96G each, 4h limit). Arrays are independent and can overlap; Crest availability controls actual concurrency. Native files remain under comparisons, never production. Each treatment starts with a new directory and independent caches. Repeated large input hashing is omitted; manifests record exact resolved input paths, sizes and mtimes. Original inputs are not modified.

## Nine treatments with common evaluation

| Array/task | Treatment | Decision question |
|---|---|---|
| hifiasm/0 | Ex2+Ex3, original settings | Reproduces baseline and retains missing graph context, concurrently with all variants |
| hifiasm/1 | Ex2 only | Does Ex2 support the same raw contig path? |
| hifiasm/2 | Ex3 only | Does Ex3 support the same raw contig path? |
| hifiasm/3 | HiFi-only | Is the connection present without Hi-C-dependent phasing? Phase completeness is expected to differ |
| hifiasm/4 | Ex2+Ex3, post-joining disabled | Does -u 0 remove the implicated path while preserving its sequence? |
| YaHS/0–1 | Existing round2 inputs, scaffold error correction on, both haplotypes | Does native scaffold correction remove the suspect physical connections? |
| YaHS/2–3 | Same, scaffold correction on and maximum resolution bin size capped at 10 Mb | Do late coarse joining rounds create/recreate the large composite connections? The cap is an experiment, not a chosen default |

A/B refers to current Ex2/Ex3, not an unverified claim about the historical library identity. The hifiasm parameter list copies the confirmed original command except for the named treatment differences. No coverage, purge, overlap or perturbation sweep is added. Full HiFi correction/overlap must run independently because original caches were not retained; the newly retained caches enable future compatible variants. Hifiasm native GFA/BIN/BED outputs stay in each treatment directory. No downstream finishing or pangenome runs.

Each hifiasm treatment maps its two new haplotypes to the same CMat hap1 reference and exports reference-scaffold composition summaries. Secondary placements remain in PAFs; composition is a MAPQ/query-span screen, not calibrated unique aligned bases or an automatic fusion call. Preserve/reference names via the copied CMat name map. Do not equate reference-scaffold composition with physical misassembly without boundary evaluation.

The five implicated original raw contigs, selected in an investigation TSV, are mapped to BOTH new haplotypes for each treatment. This avoids unnecessary whole-original-genome comparisons and does not assume hap1/hap2 labels stay fixed. The mappings localize the previous suspect sequence, changed paths and split representations. Each treatment retains its new graphs for branch/overlap inspection. The broader scaffold-level composition must still be interpreted against the full sister/peer chromosome pieces; an unchanged label transition is not proof of the same physical fusion.

YaHS uses unchanged corrected-scaffold FASTA and retained filtered Hi-C BAMs. The script retains every intermediate AGP and original command. MAPQ20 remains specified for parity with baseline; the coordinate-sorted-BAM suppression caveat remains, and upstream filtering is unchanged. Both variants enable initial contig error correction and scaffold error correction. Gap length and chromosome-count/N50 are not success criteria.

## Return once, evaluate together

Outputs are `comparisons/ctlk-fusion-batch-ARRAYID/` and `comparisons/ctlk-yahs-batch-ARRAYID/`. Every completed task writes a small `*.review.tar.gz` with commands, logs, manifests, PAF/TSV or AGP outputs, excluding large FASTA/GFA/BIN/BAM files. Return those review archives plus failure logs for any unsuccessful tasks. Retain native graphs/caches on Crest for targeted inspection, rather than transferring them all.

Evaluate one matrix across the four composites and associated chr4/7/9/12/14 pieces:

1. Whether the suspect physical adjacency/path survives, changes or splits, and at which raw/processed/phased graph or YaHS round it changes.
2. Whether alternatives preserve all component sequence and improve conformity to the six better current assemblies, resolving complementary versus duplicate chr4/7 pieces.
3. Whether existing uniquely anchored reads/contact evidence favors the surviving connection or a conservative unjoin. Do not count disappearance under a parameter change as proof of error.
4. Retain, break probable misjoin, unjoin unsupported, or close unresolved with a specific localization limitation. Exceptional haplotype-specific chromosome fusion claims require affirmative boundary/phase support.

Choose the least disruptive successful setting/correction from this batch, then run one confirmation through the affected downstream stages and finish the core. If a second batch is needed, it must address a specific contradictory result from this matrix; do not resume generic artifact hunting or a broad parameter grid. A missing graph from the old run is no longer a blocker because the baseline regeneration runs alongside the alternatives.

## Validation

Meaningful treatment-isolation and alignment-union tests are included and executed at the start of each hifiasm array task. They have not been executed locally, consistent with the cluster-validation workflow. Static whitespace checks are separate from runtime validation. Version mismatch, missing input, or a failed subprocess fails the treatment; failures are not biological absence. No production defaults or assembly cuts are changed by this batch.

Local Python AST syntax checks and git whitespace checks passed. Runtime/tool integration and biological evaluation remain cluster checks.
