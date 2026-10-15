"""Prepare a small production-finishing test with explicit external-tool doubles."""
import runpy
from pathlib import Path
runpy.run_path(str(Path(__file__).with_name("setup_misassembly_integration.py")))
root = Path("/tmp/chimera-nxf-runtime")
for sample in ("mixed", "fragmented", "unclassified"):
    (root / (sample + ".fa")).write_text(">big\nAAAACCCC\n>small\nGGGG\n>uncertain\nTTTT\n")
    (root / (sample + ".tsv")).write_text(
        "old_name\tchromosome_member\nbig\t" + ("yes" if sample == "mixed" else "no") +
        "\nsmall\tno\nuncertain\tunresolved\n")
(root / "NO_HARMONIZE").touch()
