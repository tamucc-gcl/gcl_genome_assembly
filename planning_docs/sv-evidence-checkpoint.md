# SV evidence and diagnostic plotting checkpoint

2026-10-02. Experimental comparison workflow only; production calling unchanged.

## Plot revision

Candidate plots now include local context/detail and two one-axis context panels: reference interval plus 250 kb flanks against the entire tested query chromosome, and the converse. These expose distant partners hidden by the previous rectangular crop. Full event names, native copy status and direct SyRI alignment-child counts are shown. Direct children are outlined in black and exported to TSV. CPG/CPL/TDM may have no direct alignment children; this is not a missing-evidence verdict. Bottom panels use unequal axis scales and cannot demonstrate interchromosomal matches outside the tested pair.

Only static review and git diff checks performed locally. User must render on cluster and return archive for visual verification. Reuse the successful SyRI retry; no alignment/calling rerun needed:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_diagnostics.sbatch \
  "$PWD/comparisons/environments/sv-diagnostics-r" \
  comparisons/sv-syri-retry-20261002-070635-1503156
```

## Evidence policy to implement after review

Preserve native calls and hierarchy. Separate assembly-derived calls, calls passing documented confidence checks, and read-supported calls. Filtering is not biological validation. Record unavailable evidence separately from contradictory evidence; short-read-only inputs must not silently acquire HiFi support requirements.

- Assess alignment identity, aligned fraction, competing matches and both breakpoint neighborhoods; flag gaps/contig boundaries. Use parameters calibrated against fixtures and real examples, not arbitrary universal thresholds.
- For inversions, seek a reversed block with appropriate flanking adjacencies. For translocations/dispersed duplications, inspect source and destination context; local flanks need not align to one another. Duplication calls require copy context, not merely one reverse alignment.
- Use chromosome-scale/full comparison scope before interpreting interchromosomal events. Retain ambiguity in repetitive regions; do not discard all repeat-associated variation.
- For stronger support, use original reads from each respective individual mapped competitively against the complete relevant assembly, checking breakpoint-spanning/split alignments with unique flanking anchors, coverage and competing placements. Diploid reads need haplotype-aware interpretation. Lack of uniquely placeable support can mean unassessable, not false.
- Agreement with another method or other haplotypes is corroboration, not proof. Reads used to build the assembly provide additional evidence but are not independent sequencing. Orthogonal data or targeted experimental confirmation can resolve important ambiguous events.
- No bespoke confidence score or hard filtering added in this plot-only revision. Prefer published tools for read-based calling/genotyping; benchmark event-class support before selecting one.

Sources: https://schneebergerlab.github.io/syri/fileformat.html and https://schneebergerlab.github.io/syri/pipeline.html . SyRI classifies assembly alignments; its authors recommend testing alignment settings for the biological question.
