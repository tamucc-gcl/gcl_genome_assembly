# Full-cohort evidence and report review

Run the same production CHIMERA workflow on all assemblies in the existing assessment cohort. Each assembly supplies AGP, reference alignment, graph, Hi-C pairs and its coordinate source. This run creates evidence and Markdown review packets and verifies that no assemblies were changed. It neither reassembles genomes nor supplies a historical cut fixture.

From the Crest project directory, after syncing the updated repository:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_cohort_evidence.sbatch ctlk_analysis
```

Optional second and third arguments specify the assessment root and retained HiFi BAM root. Defaults match the existing real-test inputs. Outputs are in `comparisons/chimera-cohort-JOBID/`, with a sibling `.review.tar.gz` packet. Open `results/assembly/chimeras/review/README.md`.

The launcher explicitly permits unavailable retained HiFi BAMs so the whole cohort can be assessed. `inputs/evidence-availability.tsv`, the job log and candidate measurement tables expose that limitation. Missing BAMs do not silently become zero read support. Existing BAMs must have valid coordinate-stage/checksum provenance. This harness has no raw HiFi input channel; without a retained exact BAM, local HiFi evidence is unavailable for that assembly. Supply those BAMs for complete sequence-evidence coverage.

The completion check requires every cohort assembly in the review registry, including assemblies without candidates, a nonempty assessment checksum, an individual Markdown report and a passing unchanged-assembly verification. Review selections must all remain NO.

Workshop the reports using actual cohort results: chromosome pair assignments, supported continuity, alignment observability, candidate ranges, exact gap cuts, control adequacy and user-added cuts. Ranges are measured transition intervals, not automatically chosen midpoint cuts. Poorly measured comparisons are distinct from demonstrated absence of expected homologous sequence; evidence for the latter requires checking available sequence and alignment sensitivity.
