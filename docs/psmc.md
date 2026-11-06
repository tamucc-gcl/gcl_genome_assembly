# Optional PSMC demographic history

Implementation status: integrated into `main.nf` and tested with native tools on
bounded known-input fixtures. Biological Crest acceptance is pending. This does
not change the accepted assembly, chimera or core graph procedures.

## Enable and select inputs

Set `--run_psmc true`. PSMC is independent of `run_pangenome`, `qc_mode`,
`run_post_assembly` and contact-map/visualization controls. Disabled PSMC registers
no analysis tasks. The always-on shared final linear reference supplies the same
coordinate frame whether or not a graph is built; graph VCF/GREF is not used.

`--psmc_mode auto` chooses paired-assembly mode when both same-individual
haplotypes are finalized, otherwise the complete-HiFi route for a single diploid
assembly. `assemblies` or `hifi` explicitly selects a route; `hifi` can also be
used on an individual with two assemblies for a method comparison. Haploid,
incomplete, shotgun-only or otherwise unsupported inputs get explicit skips.
Different individuals or species are never joined into a diploid pair.

By request, scope is **long reference scaffolds only**. The default minimum is
`--psmc_min_contig_bp 1000000`. Shorter sequences are excluded, not concatenated.
No sex-chromosome identification is required or inferred. The retained scope may
therefore include sex-linked sequence. The same minimum applies before optional
BED masks. Full-reference mapping is retained to allow competing placements;
read variant calling and consensus construction are restricted to long scaffolds.

An optional `--psmc_manifest file.tsv` selects individuals and overrides mode,
scope masks or scaling. It is project data, not pipeline source. Columns:

```text
sample	mode	include_bed	exclude_bed	mutation_rate	generation_time
individual_A	assemblies				
individual_B	hifi				
```

Only `sample` is mandatory. Omitted columns/blank values use global defaults.
When a manifest is supplied, only its listed individuals are analyzed. Duplicate
or unknown samples and unknown columns fail explicitly. BED paths are resolved
relative to the original manifest location; coordinates are zero-based half-open
on the exact selected final reference. Wrong names or out-of-range intervals fail.

## Scaling and fitting controls

`psmc_mutation_rate` (per nucleotide per generation) and `psmc_generation_time`
(years) both default to **null**. With both unset, output is explicitly unscaled:
time is `t * theta_per_base`, and size is `lambda * theta_per_base`. Set both to
positive finite values for years/Ne. A fractional generation time is supported.
For mixed species, scaling must be supplied per individual in the manifest;
one global species-specific assumption is not silently broadcast across taxa.

| Control | Default | Purpose |
| --- | ---: | --- |
| `psmc_bin_size` | 100 bp | Upstream fq2psmcfa binning and scaling |
| `psmc_min_good_bases` | 10,000 | fq2psmcfa minimum usable bases per sequence |
| `psmc_min_callable_bp` | 10,000,000 | Execution guard, also checked after sequence retention by fq2psmcfa |
| `psmc_min_callable_fraction` | 0.1 | Fraction of selected non-gap scope |
| `psmc_min_hets` | 1,000 | Minimum heterozygous SNP guard |
| `psmc_iterations` | 25 | Fitting iterations; output must reach the requested iteration |
| `psmc_pattern` | `4+25*2+4+6` | Upstream example, not an organism-specific optimization |
| `psmc_max_time` | 15 | Initial maximum scaled coalescent time |
| `psmc_theta_rho_ratio` | 5 | Initial theta/rho ratio |
| `psmc_plot_iteration` | 20 | Explicit selected complete iteration |
| `psmc_show_last_interval` | false | Omit terminal free-parameter block, as upstream plotting does |
| `psmc_bootstraps` | 10 | Pilot replicates; request 100 explicitly for a full set |
| `psmc_bootstrap_block_bp` | 5,000,000 | Physical block size, converted to PSMC bins for splitfa |
| `psmc_bootstrap_max_forks` | 8 | Concurrent bootstrap tasks |
| `psmc_seed` | 1729 | Derives stable per-individual/per-replicate resampling seeds |

Guard thresholds are operational floors, not guarantees of demographic power.
Insufficient callable data or an input without both homozygous/heterozygous bins
produces a **not fitted** report with retained evidence. A failed requested tool
or malformed/truncated fit fails the pipeline rather than becoming an empty plot.
Per-parameter-group expected recombination diagnostics are exported; groups below
the upstream advisory of about ten are highlighted for review, not automatically
removed. Iteration/pattern choices need biological review after the pilot.

## Genotypes and callable sequence

Paired assemblies use pinned Dipcall 0.3 with `-a` so scaffold naming does not
silently exclude sequence. Its confident BED defines where homozygous reference
can be assumed. Called SNPs supply homozygous alternate or heterozygous IUPAC
bases. Filtered/missing calls and indel/complex reference spans plus five-base
flanks remain unknown. Gaps are always unknown.

The read route maps the individual's **complete HiFi read set** using minimap2
`map-hifi`, then uses bcftools 1.21's general diploid multiallelic caller with
all-site output and the unseen allele retained. No human-trained neural model is
assumed. A site must have a complete diploid genotype, acceptable depth and
genotype confidence; heterozygotes also require allele balance. A site absent
from the VCF stays unknown. Default controls:

- `psmc_min_mapq=20`, `psmc_min_baseq=20`, `psmc_min_gq=20`;
- `psmc_min_depth=0` and `psmc_max_depth=0` derive bounds from mean depth over
  long scaffolds: floor(mean/3), at least one, and ceil(2*mean);
- explicit positive depth values override those derived bounds;
- `psmc_min_allele_balance=0.2`, `psmc_indel_padding=5`;
- `psmc_max_pileup_depth=10000` is the caller's computational depth cap.

**Confidence definition:** use FORMAT/GQ where supplied. At many monomorphic
sites bcftools emits no GQ; use the retained PL likelihood separation between the
called genotype and its best competitor instead. That separation is not posterior
GQ. `psmc_min_gq` is the phred threshold applied to the available measure; the
rule and counts for each source are recorded in `input_qc.json`. Missing both
measures makes a site unknown. This general caller is functional on known-input
fixtures; real-data heterozygote sensitivity and assembly/read-route differences
still require biological acceptance before histories are treated interchangeably.

## Scheduling, cache and evidence

Mapping, calling, masks/consensus, observed fit, split input, each bootstrap,
curve export and R/ggplot presentation are separate processes. Each replicate is
scheduled by Nextflow. Upstream splitfa defines blocks; a small adapter uses
recorded independent seeds and the upstream nearest-total-length stopping rule
to sample blocks with replacement. PSMC fits the resampled input without native
`-b`; verified `-I0` disables unseeded initialization perturbation. No upstream
source patch or runtime download is required. The selected block indices and
lengths are retained. More replicates reuse earlier stable replicate IDs.

Changing mutation rate/generation time or the selected plotting iteration reuses
genotypes, masked input and fits. Changing thresholds invalidates the relevant
mask/calling stage; changing fitting settings invalidates fits. Ordinary Nextflow
cache rules apply. Native package pins are in `environments/psmc.yaml`; plotting
pins are separate in `environments/psmc_plot.yaml`.

Published outputs are under `post_assembly/psmc/`:

- `README.md`, `summary.tsv`, `psmc_plan.json`: current-run results/skips and input selection;
- `<taxid>/<sample>/genotypes/`: VCF/index, Dipcall confident BED or read-depth table;
- `<taxid>/<sample>/input/`: masked diploid FASTA, callable BED, PSMCFA,
  per-contig QC, input status and provenance;
- `<taxid>/<sample>/fits/`: observed/bootstrap PSMC outputs and resampling records;
- `<taxid>/<sample>/report/`: Markdown, PNG/PDF, curve table, interval diagnostics
  and explicit fit/iteration/scaling manifest.

Large BAMs remain in Nextflow work. The main report links the cohort and individual
reports and shows input mode, outcome, callability and scaling. Disabled analyses
do not add a main-report section. The integration review archive collects compact
PSMC fits/QC/provenance/figures, excluding large consensus/VCF/BAM inputs.

Bootstrap spread reflects block resampling conditional on the supplied diploid
input. It does not include uncertainty in mutation rate, generation time,
heterozygote detection or the PSMC model. No confidence ribbon is invented from
ten pilot curves.

## First Crest pilot

Use the accepted core baseline and normal main-pipeline launcher/cache. This runs
all eligible individuals through automatic mode, so this project's five paired
assemblies use Dipcall. Scaling parameters remain null. It is a biological pilot,
not a completed PSMC acceptance claim.

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-freeze-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/core-psmc-01 \
  --run_report true --report_chromosome_painting true \
  --run_psmc true --psmc_bootstraps 10 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

Assembly, scaffolding, accepted cuts, graph construction and existing read evidence
should stay cached. New jobs are PSMC planning, per-individual Dipcall and input
conversion, observed fits and up to eight simultaneous bootstrap tasks, exports,
figures and dependent reports/provenance. Mappings/fits may be substantial; inspect
the first trace before making runtime promises. Return
`comparisons/core-psmc-01.review.tar.gz`.

After the paired-route pilot, compare a manifest-selected individual's complete
HiFi route against its assembly route in a distinct output directory. Do not
concurrently resume two runs using the same shared launch/session directory.

## Methods sources

- [PSMC input, fitting, scaling and bootstrap documentation](https://github.com/lh3/psmc)
- [Pinned PSMC resampling and output definitions](https://github.com/lh3/psmc/blob/0.6.5/aux.c)
- [Dipcall confident-region and genotype semantics](https://github.com/lh3/dipcall)
- [bcftools calling documentation](https://samtools.github.io/bcftools/bcftools.html)
