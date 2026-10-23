# Pinned regional assignment source

These are unchanged source files from ComparativeGenomicsToolkit/cactus-gfa-tools,
commit `f9b999e54c49a05ad6e9d44e1832c4a8421f28f2`. The upstream MIT license is included.
`manifest.json` records the archive checksum and checksums of each source file.

The production sequence policy builds a corrected tool. The optional regional assignment
diagnostic also builds the original tool for comparison. The sole source correction in `smooth_query_intervals` compares
`interval.value.second` (reference chromosome ID) with `ref`, rather than
`interval.value.first` (aligned-base count). No thresholds or interval lengths change.
The exact applied diff, build provenance and real CLI regression results are published.

The build uses `environments/rgfa_split_patch.yaml`. It does not replace the executable
inside the Cactus container. Production uses task-local wrappers, chromosome-isolated
deletion filtering and the separately configured regional chunk size (25 kb by default).
The diagnostic flags and replay modules remain independently removable. The production
build, environment and this source directory remain required until an accepted upstream
release supplies the correction and passes the same production checks.
No species, sample, scaffold or breakpoint is encoded here.
