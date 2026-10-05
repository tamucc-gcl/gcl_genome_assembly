# All-junction diagnostic test

This standalone test makes no production changes or assembly cuts. It inventories
every scaffold in every supplied assessment assembly, rather than selecting only
the previously suspicious CTlk scaffolds. Sister-haplotype discordance raises
review priority; it is not assigned a numerical fusion prior or treated as proof
of heterozygosity. Both haplotypes could have assembly errors or phase switches.

## Run

From the cluster project root after syncing:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
bash -n gcl_genome_assembly/scripts/comparisons/run_all_junction_review.sbatch
sbatch gcl_genome_assembly/scripts/comparisons/run_all_junction_review.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  --results comparisons/chimera-foundation-01 \
  --registry gcl_genome_assembly/planning_docs/investigations/ctlk-junction-registry.tsv \
  --evidence-packet comparisons/junction-assessment-20261005-084623-1506302 \
  --hic
```

Tests are supplied for cluster execution; they have not been run locally.
The default layout discovers all published round2 haplotype FASTAs. It fails on
unrecognized identities rather than guessing sisters. Check input_manifest.json
against the intended batch; an assembly absent from the directory cannot be
discovered. An explicit --manifest JSON list supports other layouts: each entry
has assembly, sample, haplotype, sister (empty when absent), and stages. Stages
use the origin-tracer schema (stage/fasta/agp, optional graph), with exactly one
assessment stage. Paths are relative to the submission directory.

The registry adds the existing flagged intervals; it is not a discovery filter.
Omit it for another dataset or provide a checksum-bound registry for that batch.
--inventory-only skips mapping/contact streaming. --hic-inputs uses the existing
assessment diagnostic TSV schema (assembly/fasta/sha256) if the layout adapter
cannot locate the exact Hi-C source FASTA.

## What runs and what it costs

- Inventory all FASTA N runs and all recoverable AGP gaps, plus adjacent ungapped
  component boundaries in the assessment AGP. N runs without verified AGP origin
  remain explicit unknowns. The scaffold inventory includes scaffolds with zero gaps.
- Earlier-stage origins currently require exact whole-scaffold correspondence,
  including reverse orientation. Changed earlier scaffolds cannot silently be
  assigned an origin. Embedded inherited N runs remain visible even when origin
  recovery fails. This is not yet complete recursive AGP provenance reconstruction.
- Extract 10 kb flanks at offsets 0, 25 and 100 kb and map them to each explicitly
  assigned sister. Targets run sequentially on one exclusive node. This is roughly
  one small-query mapping per assembly, not another genome all-versus-all alignment.
  Full PAF alternatives are retained. Missing/ambiguous mappings and flanks that
  cross other inventoried boundaries cannot count as clean sister disagreement.
- Reuse the prior assessment packet's HiFi summaries by exact ID, coordinates and
  FASTA checksum. Other joins explicitly have missing HiFi measurements. No new
  whole-genome HiFi remapping runs. These chain counts are descriptive, not a
  repeat-uniqueness or phase certificate; zero is not an absence test.
- With --hic, stream the retained pairs once per assembly and report library-specific
  contacts and alternative bins for all inventory entries using the existing
  verified coordinate adapter. These are the scaffolding reads, not independent
  evidence. Counts are not yet opportunity/distance/mappability calibrated and
  do not automatically determine a decision. This can be the longest part of the run.

## Review output and stopping rule

join_decisions.tsv is the master table. It keeps evidence, automated triage and
human decisions separate. sister_comparisons.tsv provides placement coordinates,
orientation, separation and window status; separation is not interpreted as exact
gap length. Different sister scaffolds can indicate fragmentation, not a different
biological chromosome. Support on one sister scaffold does not prove adjacency
through its intervening gaps or repetitive sequence. Historical assemblies are
not required. Other-individual and historical evidence from existing packets can
be consulted during review; this test does not remap every anchor to every peer.

Review every prioritized gap with the existing contact and molecular evidence.
Gather additional local uniquely anchored molecule/phase evidence or calibrated
contact evidence only when it can change the decision. Record RETAIN,
BREAK_PROBABLE_MISJOIN, UNJOIN_UNSUPPORTED, or UNRESOLVED with reviewer, reason,
evidence_for and evidence_against (use explicit 'none observed' where appropriate).
--reviews accepts a subset of rows from the decision table and validates identity
and checksum. It allows unjoin decisions only for verified AGP gaps. Internal
sequence cuts require a separately localized breakpoint and are not accepted here.
No cut file is exported, even after review. UNRESOLVED is a recorded outcome, not
a requirement to keep collecting evidence indefinitely.

Before production integration, review the batch, verify origins for any selected
cuts, resolve unsupported joins according to the agreed conservative assembly
policy, preserve sequence/coordinate provenance, prevent rejected joins from being
recreated, and audit fragment retention in the pangenome. Do not treat graph
placement as evidence that sample-specific adjacency has been restored.
