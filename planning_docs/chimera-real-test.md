# Isolated CTlk real-data workflow test

From `/work/birdlab/GCL/spratelloides_delicatulus_genome`, after syncing the repository:

```bash
mkdir -p logs
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_real_test.sbatch ctlk_analysis
```

The two array tasks run concurrently, each requesting 16 CPUs and 96 GB for up to 24 hours. The activated analysis environment must contain Python 3.11, numpy, matplotlib, cooler, minimap2, samtools, and tidk. Nextflow 23.10.1 is loaded as a module. Optional arguments after the environment override the foundation assessment directory and retained BAM directory.

The default inputs are `comparisons/chimera-foundation-01` and `genome_assembly/assembly/chimeras/sequence_context`. Preflight checks require the exact retained assessment FASTAs, mapping provenance, BAM indexes, source scaffolding pairs/readset manifests, AGPs, native primary graphs, and harmonization outputs. A mismatch fails rather than remapping or changing the assessment silently.

## Parallel cases

- **auto:** Runs the actual CHIMERA workflow for both CTlk haplotypes, discovers candidates and nearby review-only gap hypotheses, gathers peer/HiFi/each-library Hi-C/native-graph evidence, publishes decisions, and applies only automatically eligible cuts. Reuses HiFi BAMs. No assembly rerun is needed.
- **manual:** Applies only the previously reviewed H01 gap cut to hap1, then reruns chromosome assignment on the corrected species cohort with the same reference. The named H01 fixture belongs exclusively to the regression test; generic detection and decision code do not know H01.

All outputs go under `comparisons/chimera-real-JOBID/{auto,manual}`. Production assemblies are inputs and are not replaced. This test ends at corrected/reassigned pre-finishing assemblies; it does not run downstream polishing or finishing.

The manual case must preserve 953,570,870 hap1 bases and produce 363 records, with scaffold 1 split into 63,051,325 and 75,395,654 bp pieces. The 100 gap Ns remain on the left. Expected hap1 N50 is 69,062,049 bp. The verifier reconstructs every original parent from output pieces, checks complete coordinates and sequences, and rejects pending chromosome assignments. Reassignment may still classify a remaining composite as composite: that is different from leaving its identity pending.

The auto case has no prescribed biological cut outcome. Review local support, graph uncertainty, independent peer relationships, matched-control availability, and both Hi-C libraries. Unknown graph evidence and insufficient controls block automatic cutting. Zero broad-interval bridges do not establish an unsupported join. Different scaffold names representing fragments of the same chromosome do not count as different chromosomes. Contact-loss thresholds are explicit conservative heuristics requiring validation, not probabilities of misassembly.

## Return packets

Return `auto.review.tar.gz` and `manual.review.tar.gz` from the job directory. They include exit status, execution trace, decision tables/audits, raw control measurements, candidate plots, IGV sessions, native graph neighborhoods, chromosome assignment maps, and sequence verification. Large FASTAs/BAMs and work directories are excluded from packets but remain on Crest. Evidence coordinates refer to the original assessment FASTA; IGV sessions use that FASTA and its matching BAM.

Local verification covers Python measurement/decision/action behavior, mocked tool orchestration, and sequence reconstruction. Actual Nextflow scheduling, tools, plot rendering, and biological measurements require this Crest run.
