# Pipeline-wide R and plotting conventions

User decision, 2026-10-01: use tidyverse-style R, ggplot2 graphics and patchwork for combining plots throughout the pipeline.

- Use tidy data, named columns, native pipes, and dplyr/tidyr/readr/purrr/stringr as appropriate. Base R remains appropriate for argument handling, filesystem operations and numerical/statistical APIs without a tidyverse equivalent.
- New or revised plots use ggplot2 and patchwork rather than base graphics or gridExtra.
- Use named consistent palettes across figures. Derive sharing categories from reported thresholds and precedence, never hard-coded cohort sizes. Keep growth-quorum criteria distinct from sharing tiers.
- Plot validated outputs, keep presentation independent of expensive computations, retain PNG/PDF and relative Markdown links, and record package versions.
- Migrate retained legacy scripts when their analyses are revamped. Dormant scripts are not yet all converted; do not reactivate retired analyses just to restyle them.

The CLIP sharing figures now follow this convention. Histogram color assignments are checked against the independently generated partition table before plotting.
