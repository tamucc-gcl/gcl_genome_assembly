# gcl_genome_assembly

Nextflow DSL2 pipeline for multi-individual eukaryotic genome assembly, source-bound
manual misassembly review, assembly QC and optional core pangenome construction.

Author: Jason Selwyn · TAMUCC CORE / birdlab (GCL).

The core interface is frozen as **1.0 (2026-11-05)** against the accepted Crest
reporting integration. Execution acceptance does not imply biological completeness.

## Start here

- [Supported inputs, controls, outputs and manual cuts](docs/core-interface.md)
- [Configured default inventory](docs/parameter-defaults.tsv)
- [Verification and native acceptance](docs/verification.md)
- [Freeze changes and final cached check](docs/core-freeze-2026-11-05.md)
- [Maintained source layout and retired components](docs/repository-layout.md)
- [Pangenome tool corrections and biological rationale](docs/pangenome-tool-changes.md)
- [Canonical completion checklist](planning_docs/remaining-pipeline-work-2026-10-15.md)
- [Agreed optional extensions](planning_docs/post-core-optional-analyses-plan.md)
- [R/ggplot2 plotting convention](planning_docs/r-plotting-conventions.md)

## Quick start

Use Nextflow 23.10.1, conda and the configured container runtime. The `slurm` profile
is tailored to Crest; adapt resource, scratch, cache and database paths for another site.
Per-task software requirements are declared in configuration and environment files.
Mutable image defaults are recorded in provenance; override with audited immutable
images or SIF files when reproducibility requires it.

```bash
nextflow run main.nf -profile slurm -entry VALIDATE_INPUTS \
  --sample_sheet /project/data/samples.csv --hic_readsets /project/data/hic.csv

nextflow run main.nf -profile slurm -resume \
  --sample_sheet /project/data/samples.csv --hic_readsets /project/data/hic.csv \
  -params-file /project/data/assembly-params.json --outdir /project/results
```

Omit the optional Hi-C table when unused. Taxid/species and biological settings belong
in the sample table or project parameter file. Relative read paths resolve against
that input table. Keep project data, cut sheets and run outputs outside this repository.

| Input | Assembly route | Notes |
| --- | --- | --- |
| HiFi + Hi-C | hifiasm, YaHS, refinement/finishing | Library identity retained for evidence |
| HiFi only | hifiasm, refinement/finishing | No Hi-C scaffolding |
| Paired shotgun reads | SPAdes, supported conditioning | Collapsed/ploidy handling and organelle status are explicit |

Hi-C-only rows are explicit skips. Linked-read scaffolding and shotgun + Hi-C
scaffolding are not supported. Phased fragmented assemblies may contribute to a
species graph without serving as chromosome voters or reference candidates.

## Evidence, review and reports

The first pass produces an editable **unselected** review TSV and Markdown evidence
reports. To make a cut, copy the TSV outside generated output, review the evidence,
provide exact original-scaffold coordinates and reviewer/reason fields, then resume
with `--chimera_break /absolute/path/reviewed.tsv`. Undetected manual cuts are allowed.
The cutter validates the source fingerprint and verifies exact sequence reconstruction;
finishing and graph inputs consume the corrected sequences. Automatic cutting is deferred.

The main report is `assembly_report.md`. It links final FASTAs/name maps, organelles,
individual-based quality tables, chromosome structure, decisions and detailed reports.
The combined regional evidence figure and the placement-confidence figure both remain;
per-library junction contact maps are separate. The interactive QC rankings HTML remains
available. Disabled analyses do not create empty report sections.

Core pangenome outputs include construction/identity audits, sequence accounting,
sharing/growth, regional attribution, descriptive comparisons and variant export audits.
Graph representation flags are not automatic assembly-edit recommendations. Biological
interpretation of private sequence and structural variants remains an extension.

## Resume and validation

Use normal Nextflow cache handling. Preserve the launch session, work directory and
published-symlink targets. Do not expect recovery after deleting cached scientific inputs.
Select a new database snapshot directory to refresh stored databases.

```bash
python3 scripts/check_repository.py
python3 -m unittest discover -s tests -p 'test_*.py'
python3 py_scripts/validate_report.py --root /project/results
```

The [verification guide](docs/verification.md) distinguishes controlled fixtures, actual
rendering checks and native scientific-output acceptance. The latest accepted reporting
run is recorded in [detailed-report-01](planning_docs/detailed-report-01-assessment.md).
The codebase no longer carries executable obsolete detector branches or report mockups;
retired sources remain available in Git history.
