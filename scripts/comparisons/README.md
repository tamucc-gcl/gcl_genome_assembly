# Comparison helpers

Keep temporary review, collection and comparison scripts in this directory. They are separate from production pipeline modules.

Run helpers from the cluster project root (the directory containing `gcl_genome_assembly`, `data`, `work`, and the result directories). Store all generated comparison artifacts under that project's `comparisons/` directory. Use a descriptive, timestamped prefix for each comparison; preserve earlier outputs until manually cleaned up.

## Chimera comparison collector

```bash
bash gcl_genome_assembly/scripts/comparisons/collect_chimera_comparison.sh \
    tst/genome_assembly_store \
    genome_assembly
```

This collects existing reports without running analyses or changing pipeline results. Outputs are:

- `comparisons/chimera-comparison-<timestamp>/`: inventory and provenance notes.
- `comparisons/chimera-comparison-<timestamp>-reports.tar.gz`: upload this first.
- `comparisons/chimera-comparison-<timestamp>-alignments.tar.gz`: existing reference alignments, when available; retain for follow-up.

The result-directory arguments are project-relative. The defaults above are specific to the current comparison, not sample-specific analysis rules. Each future helper should document its inputs and follow the same output-directory convention.

Existing archives created before this organization change remain valid and need not be recreated. Do not place assembly work/cache directories inside `comparisons/`.
