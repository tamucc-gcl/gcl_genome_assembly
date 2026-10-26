# Native GetOrganelle retry and fragmented annotation correction

The CMat_061 native pilot completed its reads-based GetOrganelle run in about
818 seconds and recovered six linear mitochondrial scaffolds. The optional
from-assembly retry ended immediately after graph slimming with process exit 0,
but without the closing `Thank you!` marker. Our wrapper rejected it.

Inspection of the locally retained pinned 1.7.7.1 source explains the legitimate
early-return branch: when `slimmed_assembly_graph.fastg` is empty, `main()` returns
the string `finished with no target organelle contigs found!` before the normal
closing logger runs. The script entry point does not print that returned string.
The pasted log alone cannot prove the artifact state; the updated wrapper verifies
it on Crest before accepting the outcome.

The wrapper now accepts this branch only for `get_organelle_from_assembly.py`,
with the completed slimming log, a nonempty processed input graph and an existing
empty slimmed output graph. Logged errors/tracebacks and nonzero exit statuses
still fail, as do missing or nonempty artifacts without normal completion. A
verified empty retry is logged explicitly and the existing reads-based fragments
remain available for selection. It is not a failure of the primary recovery.

The six-record output would next hit the known single-record MITOS adapter limit.
That limit is now removed: MITOS 2.1.10's numbered per-record output directories
are converted in input order, with a required nonempty BED for every fragment.
Every sequence and its own features are preserved in the multi-record GenBank;
the combined BED retains original FASTA IDs. No longest-fragment selection,
invented concatenation or circular label is applied to fragmented output. A plain
`--bed` converter input still requires a single FASTA record, and a multi-record
input cannot be labelled circular. Missing per-fragment annotations remain errors.

The full native MITOS biological result still requires Crest acceptance. Local
workflow fixtures use external-tool doubles and real Biopython conversion; they
establish routing/data-preservation contracts, not native annotation accuracy.

Validation: all 26 workflow scenarios passed, including empty-target retry,
multi-record annotation, missing BED, executable/database/annotation failures,
missing context and taxonomy success/failure. Five completion-check unit tests
and one fragment preservation/matching-annotation regression passed. The latter
checks distinct sequences/features, linear topology, combined BED IDs and missing
per-fragment annotation rejection. `git diff --check` passed.

After syncing, resume only the existing pilot submission; preparation is complete:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_shortread_smoke.sbatch \
  comparisons/core-shortread-native-01 assembly
```

The modified recovery task will rerun under normal Nextflow caching. Unchanged
successful tasks can resume. Return the same `.assembly.review.tar.gz` packet.
