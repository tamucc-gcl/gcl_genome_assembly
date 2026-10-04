# Targeted join-origin and junction checks

This batch follows the successful foundation run 1506189. It is a standalone
comparison job, so it does not invalidate assembly or harmonization caches or
alter production defaults. It consumes all existing inferred-chromosome candidate
intervals without selecting a particular sample, chromosome or midpoint.

## Implemented

- An ordered stage inventory covering raw contigs, organelle filtering, purge,
  contig correction, round-one input/output, scaffold correction and assessment.
  Missing stages are explicit; an optional JSON manifest supports other layouts.
- Assessment FASTA SHA256 verification against called tables; stage FASTA and
  AGP checksums; retained HiFi provenance must match the assessment assembly/hash.
- Targeted base-resolved minimap2 alignments with secondary placements retained.
  Full-window presence requires direct comparison of sequence letters, including
  reverse complements. Gapped approximate PAF spans are not accepted as exact
  coordinate transformations. Exact CIGAR anchor placements are also retained
  when the entire interval does not match.
- Component ends and complete gaps projected only through literal verified
  full-window matches. AGP layout/length checks and local N-gap checks fail on
  disagreement. Multiple exact placements remain ambiguous.
- A bounded inventory of hifiasm/YaHS task roots named in the run trace, to locate
  retained raw graphs/intermediate AGPs without searching every work directory.
- Contig GFA record counts and links involving exact-match source segments.
  This does not reconstruct unitig paths or read-overlap graphs.
- Junction checks at actual N gaps, verified source boundaries and the entire
  uncertainty interval. No midpoint is invented as a proposed breakpoint.
- Two-sided anchor placement screening within the assessment haplotype. Secondary
  mappings are retained, and possible secondary-output saturation is flagged.
  A single emitted placement is NOT certified genome-wide uniqueness.
- Reuse of retained regional HiFi SAMs: distinct molecule counts, full two-anchor
  bracketing, CIGAR-chain support and alignment-wide NM screening. Large internal
  indels, supplementary/secondary records, duplicate/QC-fail records and unknown
  MAPQ cannot be promoted to the screened molecule count. Missing SAMs remain
  unavailable; there is no new read-mapping job.
- Pooled Hi-C contacts from each candidate's external flanks to other scaffolds,
  including zero-count flank rows. These descriptive counts are not normalized
  support ratios or independent confirmation of a scaffold built with those reads.
- Small report archive and nonzero exit status on failure; FASTA inputs/symlinks
  stay on the cluster. Scope audit retains the reference's unassessed status.

The current read screen uses MAPQ >=20, alignment-wide NM rate <=2%, and chains
without indels larger than 50 bp. These are descriptive sensitivity settings
inherited from the previous support exploration, not calibrated classifier
thresholds. Local read lengths are sampling-biased and do not estimate expected
bridge counts. Physical scaffold-gap lengths are not inferred from placeholder Ns.

## Run on Crest

Sync the code first, then run from the project root. No output moves or full
Nextflow rerun are required. Existing production results are used only to find
the prior HiFi-enabled context; the new comparison supplies the assembly frame.

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'

sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_origins.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  comparisons/chimera-foundation-01 \
  --hifi-results genome_assembly
```

The first argument must be an existing Conda environment containing Python,
samtools and minimap2. The previously used SyRI benchmark environment should
contain these tools; the launcher verifies executables before running. If the
HiFi-enabled results were moved, replace `genome_assembly` with their retained
results root. Wrong assessment hashes are rejected, not silently reused.

This runs stages sequentially with four threads on an exclusive node and a
96 GB request, following the node-sharing issue observed on Crest. It does not
rerun hifiasm, YaHS, harmonization or HiFi mapping. Eight hours is the scheduling
limit, not an estimated runtime. Index construction, stage I/O and streaming the
large contact files are the main costs. The single-index resource budget is
configurable with `--index-bases`; it is not a chromosome-size eligibility rule.

## Return files

The launcher automatically creates and prints:

```text
comparisons/chimera-origins-YYYYMMDD-HHMMSS-JOBID.tar.gz
logs/chimera-origins-JOBID.out
```

Attach both. If manual collection is necessary, use the exact directory printed
by the job (replace RUN_DIRECTORY and JOBID):

```bash
out=comparisons/RUN_DIRECTORY
tar --exclude='*.fa' -czf "${out}.tar.gz" "$out"
ls -lh "${out}.tar.gz" "logs/chimera-origins-JOBID.out"
```

Keep the unarchived comparison directory and shared work/cache. The report
archive retains targeted PAFs with CIGARs, input hashes, small graph metadata,
origin/boundary/anchor/read/contact tables and tool versions. No large genome
FASTA, GFA, BAM or Hi-C pairs upload is required.

## How the result determines the next experiment

1. Exact full-window presence in a raw contig establishes that the represented
   adjacency predates YaHS. Inspect the raw contig's internal repeat/graph/read
   support rather than modifying a later scaffolding gap.
2. Presence only at a later stage is a clue, not proof of creation there. Use
   detailed PAFs and exact flanks to distinguish a new component join from a
   sequence correction, clipping, missing input or ambiguous placement.
3. Retained source boundaries provide concrete junction hypotheses. Missing
   boundaries leave an uncertainty interval; they do not justify its midpoint.
4. Compare current and alternative connections at those hypotheses using
   competitive read placements, repeat-aware anchors and library-specific Hi-C.
   The pooled counts and existing regional reads produced here help choose those
   tests but cannot replace them.

This is the origin-tracing and descriptive evidence layer, not the completed
fusion/misassembly classifier. Independent reference assessment, haplotype-aware
competitive mapping, calibrated read opportunity, matched Hi-C controls and
automatic editing remain separate stages. Every action in this batch is REVIEW.

## Validation and method references

Added synthetic regression coverage for reverse coordinates, indels, malformed
AGPs/tables, direct sequence disagreement, repeated placements, stale provenance,
read-name deduplication, false CIGAR bridging, contact endpoint ordering and a
mocked end-to-end origin run. Python/R/Nextflow tests were not run locally; execute
the cluster suite above before submitting the diagnostic.

Coordinate/format behavior follows the [NCBI AGP specification](https://www.ncbi.nlm.nih.gov/genbank/genome_agp_specification/).
Base alignments and secondary-reporting limitations follow the
[minimap2 documentation](https://github.com/lh3/minimap2) and
[manual](https://github.com/lh3/minimap2/blob/master/minimap2.1).
