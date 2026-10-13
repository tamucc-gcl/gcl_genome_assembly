# Observable regional Hi-C scoring (v7)

The regional search retains the existing 150 kb radius, 100 within-side pairs
per flank requirement, shoulder matching and depletion threshold. It now chooses
the minimum ratio among observable bins first. A sparse zero cannot hide an
observable zero or override an observable nondepleted ratio. Equally low
observable bins are ranked by distance to the assay, then coordinate, making
selection independent of input order.

If no observable bin exists, the sparse minimum is retained with observable=false
and cannot provide affirmative contact evidence. Sparse bins, including undefined
ratios, remain in unobservable_bins. Selected cross/within-side counts and the
selection method are stored in evidence.json. The report identifies the selected
regional position, retains sparse-zero limitations, and shows both local calibrated
and regional findings when both exist. No exact cut is inferred from a contact bin.

343 tests passed, including shifted-coordinate and shuffled-order selection,
sparse-only zeros, sparse zeros beside observable nonloss, undefined ratios,
regional corroboration without an exact cut, and joint local/regional reporting.

Replaying all ten v6 packets changed only CTlk hap2 C02 from REVIEW to SUSPECT.
Its priority region remains 39.703–39.708 Mb and broader decision range remains
39.684–39.736 Mb; this is regional advice, not an exact cut. Ex3 now contributes
the observable 39.8 Mb bin (within left=1068, within right=599, cross=0). Sparse
39.6/39.7 Mb bins remain limitations. Ex2's observable 39.8 Mb bin does not meet
the regional depletion threshold. All manual selections remain NO. C01's exact
manual gap recommendation and all other decisions are unchanged.

For a full cached Crest rerun after syncing:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/chimera-observable-contact-v7 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-manual-review-v6/results/assembly/chimeras/mapping_manifest.json"
```

The local replay uses retained measurements; it does not claim a new Crest run
or apply cuts. The full pipeline command reuses existing mappings and applies
no cuts unless an explicit reviewed file is supplied.
