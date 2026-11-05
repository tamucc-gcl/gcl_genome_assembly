# main-report-02 assessment — 2026-11-05

## Crest outcome

Nextflow completed successfully: 939 tasks, 928 cached and 11 completed, with no
failed tasks or retries. Cactus construction remained cached. The outer launcher
returned 1 because the post-publication validator rejected the valid `contact_maps/`
directory link in assembly_details.md. Its recorded result checked 361 links and
reported only that error. The subsequent core integration audit therefore did not
execute; its missing output is not a completed PASS.

The validator now accepts existing directories for navigation hyperlinks while
requiring actual nonempty files for images and registered artifacts. Missing
directories still fail. Eight validator tests pass.

## Retained results reviewed

- All 133 expected report-output checks are PRESENT.
- Ten final name maps, ten painting tables, the complete QC metrics CSV,
  assembly eligibility JSON and three pangenome audit JSONs are byte-identical to
  main-report-01. This does not independently hash full FASTAs or graphs, which are
  excluded from the review packet.
- Four manual cuts are applied and verified; two chr7/chr12 events are explicitly
  retained by review. Main-report decision tables preserve that distinction.
- Five Needs attention rows concern pangenome graph representation only, with the
  actual <80% screen and measured ratios. Pla scope limitations do not appear there.
- Five mitogenome audits each contain 37 loci: 13 CDS, 22 tRNAs and two rRNAs.
- The combined chromosome structure figure was visually inspected. Chromosome
  ordering, paired haplotypes, selected-final-reference painting, displayed assembly
  fractions and separately scaled unplaced totals are present. The rejected
  fragmented-assembly overview was not introduced.
- Quality/download tables and themed detailed-evidence navigation follow the agreed
  layout; contact resolutions are grouped by assembly and ribbons follow dotplots.

An additional navigation check exposed a missing chimera link in the assembly run
summary. Its staged-path name comparison was false despite an available README.
Matching the staged path suffix fixes it; the real Nextflow cohort wiring test now
passes, including the summary link and plots-disabled branch.

## Finish the native postchecks without recomputing assemblies

After syncing the validator fix, run from the Crest project directory:

```bash
python3 gcl_genome_assembly/py_scripts/validate_report.py \
  --root comparisons/main-report-02/results

python3 gcl_genome_assembly/scripts/comparisons/audit_core_integration.py \
  --baseline comparisons/main-report-01/results \
  --results comparisons/main-report-02/results \
  --trace comparisons/main-report-02/trace.tsv \
  --out comparisons/main-report-02
```

These commands check the existing published results. The small run-summary link
fix and new detailed figures will appear on the next normal cached report run.

## Detailed figures started

The first addition is an R/ggplot2 junction-placement figure in each event's detailed
chimera report. It consumes retained measurements; it performs no mapping, candidate
discovery, ranking or cut classification. The adapter verifies credible counts against
the existing assay's predicate and retained summary rather than silently replacing
those measurements.

The figure pairs regional primary-spanning curves with probe-level counts for all
primary spans, high-MAPQ spans, credible spans, unique-arm bridges, same-haplotype
competition and other-haplotype matches. Categories overlap and are not stacked.
Missing assays are NA; measured zeros stay zero. Local read-length feasibility is
shown, and other-haplotype matches are explicitly context rather than a veto.
Probe markers and submitted cuts are distinct. Counts, PNG/PDF, R session and package
versions are retained. The plot enablement behavior is preserved.

All six actual CTlk events rendered using the pinned environment. The report CLI
also rendered both actual CTlk packets with reviewed decisions and verification
records. Figures were inspected, and PDF Unicode font warnings were corrected with
the Cairo device. Generic tests cover missing assays, measured zeros, same- versus
other-haplotype competition and disagreement with retained credible counts.

This starts the detailed pass; it does not finish migration. Existing overview and
contact figures remain until their R replacements preserve their accepted tracks
and evidence. Next: the combined chromosome/HiFi/Hi-C/telomere overview and separate
contact maps, followed by graph-sharing figures and remaining detailed QC layout.
Do not request an assembly or graph rebuild solely for presentation changes.
