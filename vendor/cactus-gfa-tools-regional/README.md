# Pinned regional assignment source

These are unchanged source files from ComparativeGenomicsToolkit/cactus-gfa-tools,
commit `f9b999e54c49a05ad6e9d44e1832c4a8421f28f2`. The upstream MIT license is included.
`manifest.json` records the archive checksum and checksums of each source file.

The optional regional assignment diagnostic builds the original tool and a corrected
copy. The sole source correction in `smooth_query_intervals` compares
`interval.value.second` (reference chromosome ID) with `ref`, rather than
`interval.value.first` (aligned-base count). No thresholds or interval lengths change.
The exact applied diff, build provenance and real CLI regression results are published.

The build uses `environments/rgfa_split_patch.yaml`. It does not replace the executable
inside the Cactus container. The diagnostic flag, module, adapter, environment and this
directory can be removed together when an accepted upstream release supplies the fix.
No species, sample, scaffold or breakpoint is encoded here.
