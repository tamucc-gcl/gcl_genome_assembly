# Optional PSMC demographic history

Implementation status: **accepted**, following native checks, the completed
HiFi pilot and full five-individual main-pipeline cohort acceptance. Assembly-route
defaults are `asm20` and a 10 kb minimum alignment; both remain configurable.

The accepted default policy is Dipcall for two finalized haplotype assemblies
from the same diploid individual, and complete-HiFi mapping/calling against an
own reference when only one assembly representation is available. Both routes
remain supported; `psmc_mode` can explicitly select either applicable method.
Here a single haplotype assembly still represents a **diploid biological
individual**. Standard PSMC is not enabled for biologically haploid samples.

The difference between the pilot's assembly-derived and read-derived histories
is an **open biological/methodological question**, deferred rather than an
implementation acceptance blocker. Neither route is asserted to be a truth set
or interchangeable with the other. Matched-callable-scope fits may be revisited
later; no additional PSMC investigation is required to close this implementation
pass. Preserve method, mask and parameter provenance when reporting results.
Accepted assembly, chimera and core graph procedures are unchanged.

## Enable and select inputs

All settings below live in the normal `nextflow.config` params block and may be
set on the command line or in a Nextflow parameter file. Key current defaults:

| Parameter | Default | Effect |
| --- | --- | --- |
| `run_psmc` | `false` | Enable PSMC independently of graph/post-assembly switches |
| `psmc_mode` | `auto` | Same-individual paired assemblies, or supported complete-HiFi fallback |
| `psmc_assembly_preset` | `asm20` | Minimap2 preset for both PAF/masks and SAM/variant calls |
| `psmc_assembly_min_alignment_bp` | `10000` | Minimum assembly-alignment length in both filters |
| `psmc_min_contig_bp` | `1000000` | Reference scaffold inclusion; distinct from alignment length |
| `psmc_bootstraps` | `10` | Pilot replicate count; explicitly request 100 for a full analysis |
| `psmc_bootstrap_max_forks` | `8` | Concurrent bootstrap Slurm tasks; increase within available resources |
| `psmc_read_call_workers` | `8` | Concurrent whole-scaffold callers for the HiFi route |
| `psmc_mutation_rate`, `psmc_generation_time` | `null`, `null` | Leave both unset for unscaled results |

The remaining mask, depth, quality, interval-pattern and fitting controls remain
in the same params block. For example, `--psmc_assembly_preset asm10
--psmc_assembly_min_alignment_bp 20000` overrides the new defaults without
editing pipeline code. A normal fresh run does not need comparison directories,
cache manifests or a project-specific PSMC selection sheet.

Set `--run_psmc true`. PSMC is independent of `run_pangenome`, `qc_mode`,
`run_post_assembly` and contact-map/visualization controls. Disabled PSMC registers
no analysis tasks. Each individual uses one of its **own finalized assemblies**
as its reference. With a finalized hap1/hap2 pair, hap1 is chosen deterministically;
with a sole assembly, that assembly is used. This does not assert that hap1 has
better quality. PSMC does not use the species graph reference, graph VCF or GREF.
The `PSMC_REFERENCE` jobs verify the selected FASTA against its FAI and record
the assembly ID, individual, selection reason, SHA256 and sequence dictionary.
The planner rejects a reference belonging to another individual. Changing the
species graph's reference does not change PSMC reference selection.

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
on that individual's exact selected final reference. Recheck any BED prepared
against a different reference: identical chromosome names do not imply identical
coordinates. Wrong names or out-of-range intervals fail.

## Scaling and fitting controls

HiFi pileup/calling now runs whole reference scaffolds concurrently within each
`PSMC_READ_CALL` task. `psmc_read_call_workers` defaults to 8; actual concurrency
is capped by scaffold count and half the allocated CPUs (one pileup/call pair
per worker). The process requests up to twice the requested workers in CPUs,
bounded by `max_cpus`, and budgets 24 GB per concurrent worker after an observed
18.1 GB serial peak on Crest. Normal process configuration can override this
resource request for other data. Entire
scaffolds preserve calling context and avoid internal chunk-boundary artifacts.
All mapping, quality, depth and diploid-calling options are unchanged. Outputs
are concatenated in reference-FAI order and indexed. The published
`genotypes/calling_manifest.json` and `calling.log` record scaffold commands,
timings, concurrency and stderr. Individual shard VCFs stay internal to the task.
Native tests require record-for-record equality with the previous serial command,
including an empty scaffold and nonalphabetical reference order. Missing or
failed shards fail the task rather than producing a partial successful VCF.

Changing this implementation invalidates read calling and downstream results on
the next resume, while unchanged HiFi mappings remain reusable. It does not
alter an already-running task. The parallel approach is the default for future
runs; setting one worker provides a serial comparison with the same logic.

`--psmc_assembly_preset asm5|asm10|asm20` selects the minimap2 preset for the
paired-assembly route. The default is `asm20`, following the Mat reliability and length pilot.
The selected preset applies consistently to both PAF/mask and SAM/VCF targets,
without editing installed Dipcall files or changing its confidence thresholds.
The exact make override is retained in `genotypes/dipcall_settings.json` and the
command log; assembly-route metadata also records the preset. Changing it
invalidates the assembly caller and its downstream inputs/fits, not HiFi mapping.
See the [bounded preset comparisons](psmc-mask-investigation.md) for validation.

`--psmc_assembly_min_alignment_bp` defaults to **10000**. It sets both the
paftools mask/variant alignment threshold and the Dipcall SAM-filter alignment
threshold. Values must be integers at least 10000. The 10 kb coverage/ambiguity
threshold and MAPQ 5 remain unchanged. The helper validates and edits exactly
two mask and two SAM-filter recipes in the run's generated makefile; installed
package files are untouched. `genotypes/dipcall.makefile.diff` and
`dipcall_settings.json` retain the precise changes and threshold values.
Never enlarge a mask alone while retaining calls made with a different threshold.
The assembly caller requests 192 GB with asm20 and 96 GB otherwise, calibrated
from observed Crest peaks; normal Nextflow process configuration can override
resource requests for other genome sizes. The defaults are general-purpose
parameters, not sample-specific rules; use alternative values for sensitivity
checks or other datasets. The five-individual cohort acceptance passed with these defaults.

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
| `psmc_pattern` | `1+1+1+1+25*2+4+6` | Fit the four youngest atomic intervals separately (31 free size parameters) |
| `psmc_max_time` | 15 | Initial maximum scaled coalescent time |
| `psmc_theta_rho_ratio` | 5 | Initial theta/rho ratio |
| `psmc_plot_iteration` | 20 | Explicit selected complete iteration |
| `psmc_show_last_interval` | false | Omit terminal free-parameter block, as upstream plotting does |
| `psmc_bootstraps` | 10 | Pilot replicates; request 100 explicitly for a full set |
| `psmc_bootstrap_block_bp` | 5,000,000 | Physical block size, converted to PSMC bins for splitfa |
| `psmc_bootstrap_max_forks` | 8 | Concurrent bootstrap tasks |
| `psmc_seed` | 1729 | Derives stable per-individual/per-replicate resampling seeds |

Guard thresholds are operational floors, not guarantees of demographic power.
The upstream converter also requires **at least 20% callable bases per whole
input sequence**, independently of `psmc_min_good_bases`. This fixed rule is
audited in `conversion_audit.tsv`; lowering a global guard does not override it.
Insufficient callable data or an input without both homozygous/heterozygous bins
produces a **not fitted** report with retained evidence. A failed requested tool
or malformed/truncated fit fails the pipeline rather than becoming an empty plot.
Per-parameter-group expected recombination diagnostics are exported; groups below
the upstream advisory of about ten are highlighted for review, not automatically
removed. Iteration/pattern choices need biological review after the pilot.

### Individual observed-fit diagnostics

Each individual report includes `psmc_qc.tsv` with input callability, masking,
heterozygosity and bin occupancy. Reference coverage by usable bins uses the
entire selected reference assembly as the denominator and accounts for partial
terminal bins. Bins marked K or T can contain masked bases. In assembly mode the
base-level audit counts any callable consensus bases outside the Dipcall confident
BED before conversion. The input table is also available for individuals that
do not pass the fit guards.

Fitted individuals also get `psmc_intervals.tsv`, `psmc_iterations.tsv` and three
PNGs embedded in the report showing the observed history colored by interval
support, recombinations across time and log-likelihood by iteration. Interval
metrics use the selected plotting iteration and include the open-ended terminal
group. Expected recombinations are summed across each free interval group.
Counts below 10 are insufficient, 10 to below 20 warrant caution, and counts of
20 or more define the reported temporal endpoints. These endpoints can span
intervening groups with weaker support. The observed fit is retained in full.
Convergence metrics compare the final iteration to five rounds earlier, or to
the earliest available round for shorter fits. PSMC's recorded likelihood
precedes the parameter update for each round. No convergence cutoff is applied.

Screening flags are for review. Input flags use 50% selected-scope callability,
50% masking among retained bins and a per-base heterozygosity range of 0.01–2%.
The recent-peak screen checks the first five groups for a local maximum at
least twice the size of both immediate neighbors. For demographic interpretation,
compare plausible alternate patterns on the same input when a peak matters.
Population structure can also produce peak-and-decline patterns.
See [Hilgers et al.](https://doi.org/10.1016/j.cub.2024.09.028).

## Genotypes and callable sequence

Paired assemblies use pinned Dipcall 0.3 with `-a` so scaffold naming does not
silently exclude sequence. Its confident BED defines where homozygous reference
can be assumed. Called SNPs supply homozygous alternate or heterozygous IUPAC
bases. Filtered/missing calls and indel/complex reference spans plus five-base
flanks remain unknown. Gaps are always unknown.

Dipcall's confident BED is stricter than alignment coverage: each haplotype
must have one qualifying >=50 kb alignment, with no other >=10 kb alignment
covering that position. Repeats, fragmented mappings and competing placements
can therefore remove aligned sequence. `mask_audit.tsv` records, per scaffold,
selected non-gap scope, confidence-mask losses, simple-genotype losses,
complex-variant losses and duplicate-position losses. These are sequential,
non-overlapping loss counts; no mapping/quality thresholds are relaxed to improve
retention. Read-mode audits report callable bases after genotype checks rather
than inventing a Dipcall confidence mask.

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
- `combined_history.png/pdf`: observed individual curves together, with logarithmic
  time and size axes; species/scaling groups have separate panels, and own-reference
  IDs do not split comparable individuals into separate panels;
- `individual_histories.png/pdf`: companion panels with bootstraps and individual
  size scales; unfitted/skipped individuals are explicitly listed without curves;
- `combined_curves.tsv`: unchanged exported curve values with individual identities;
- `<taxid>/<sample>/genotypes/`: VCF/index, Dipcall confident BED or read-depth table;
- `<taxid>/<sample>/reference/`: verified individual reference JSON and Markdown;
- `<taxid>/<sample>/input/`: masked diploid FASTA, callable BED, PSMCFA,
  per-contig QC, input status and provenance;
- `<taxid>/<sample>/input/mask_audit.tsv`, `conversion_audit.tsv`: stage-by-stage
  mask losses, retained scaffold names and converter rejection reasons;
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

## First Crest pilot (historical)

`core-psmc-01` completed with a shared species reference. That reference policy
departed from the intended own-individual workflow and has been corrected.
Retain that pilot as comparison evidence; use a new output directory below.

After the paired-route pilot, compare a manifest-selected individual's complete
HiFi route against its assembly route in a distinct output directory. Do not
concurrently resume two runs using the same shared launch/session directory.

## Pilot callability audit and combined figures

The first pilot completed operationally but retained only one scaffold in three
individuals and none in another. The initial 9–21% callability needs attribution;
it is not accepted as an intrinsic limitation of these assemblies. See the
[assessment and implemented diagnostics](../planning_docs/core-psmc-01-assessment.md).
After syncing, resume the production pipeline with own-individual references,
mask-loss counters and cohort figures, retaining current quality thresholds:

```bash
cd /work/birdlab/GCL/spratelloides_delicatulus_genome
sbatch gcl_genome_assembly/scripts/comparisons/run_core_integration.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-psmc-01/results \
  data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv \
  comparisons/core-psmc-02 \
  --run_report true --report_chromosome_painting true \
  --run_psmc true --psmc_bootstraps 10 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json"
```

Assembly/scaffolding, accepted manual cuts and core graph work should stay cached.
PSMC reference verification runs per individual. Dipcall calls must rerun for
changed reference inputs, followed by masking, fits and plots; ordinary Nextflow
cache rules still apply to unchanged inputs. Scaling parameters remain null.
Return `comparisons/core-psmc-02.review.tar.gz`. Compare pre/post mask losses and
retention before interpreting histories. Own references remove an unnecessary
between-individual comparison but do not guarantee adequate callability.

## Methods sources

- [PSMC input, fitting, scaling and bootstrap documentation](https://github.com/lh3/psmc)
- [Pinned PSMC resampling and output definitions](https://github.com/lh3/psmc/blob/0.6.5/aux.c)
- [Dipcall confident-region and genotype semantics](https://github.com/lh3/dipcall)
- [bcftools calling documentation](https://samtools.github.io/bcftools/bcftools.html)
