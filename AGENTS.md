# Plotting

Use R and ggplot2 for new or revised pipeline figures wherever practical. Prefer
established R plotting packages (for example pafr and gggenomes for synteny)
instead of replacing them with Python renderers. Pin and test plotting dependencies
without installing packages dynamically during a pipeline task.

Keep evidence generation separate from presentation. A plotting change must not
silently change alignment mode, filtering or biological interpretation. Verify
real rendered figures and retained evidence, not just successful process exits.

Existing Python figures are candidates for migration during their scheduled
reporting/plotting revision. Preserve their accepted evidence calculations while
migrating the renderer; do not claim the whole pipeline is R-only until audited.
