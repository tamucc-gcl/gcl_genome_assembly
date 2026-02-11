#!/usr/bin/env python3
"""
Inject BUSCO summary metrics into an assembly-stats JSON.

- Input: assembly-stats JSON produced by asm2stats.pl
- BUSCO input: BUSCO short_summary*.txt OR short_summary*.json
- Output: updated assembly-stats JSON with a top-level "busco" object:
    { "C": <pct>, "D": <pct>, "F": <pct>, "M": <pct>, "n": <total_buscos> }

Usage:
  python busco_inject.py \
    --assembly-stats assembly.assembly-stats.json \
    --busco short_summary.something.txt \
    --out assembly.assembly-stats.withbusco.json

  python busco_inject.py \
    --assembly-stats assembly.assembly-stats.json \
    --busco short_summary.something.json \
    --out assembly.assembly-stats.withbusco.json
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Tuple, Any


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, obj: Any) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=False)
        f.write("\n")


def pct(numer: float, denom: float) -> float:
    return round((numer / denom) * 100.0, 3)


def parse_busco_short_summary_txt(path: Path) -> Dict[str, float]:
    """
    Parse BUSCO short_summary*.txt from BUSCO v3-v5.

    We try, in order:
      1) Counts line like:
         C:1234[S:1200,D:34],F:12,M:45,n:1291
      2) Percent line like:
         C:98.1%[S:95.0%,D:3.1%],F:0.9%,M:1.0%,n:1291
         (then we compute n from the n field; if only percents are present,
          n is still present; counts might not be)
      3) Older lines like:
         "Complete BUSCOs (C)           : 1234"
         "Complete and single-copy (S)  : 1200"
         "Complete and duplicated (D)   : 34"
         "Fragmented BUSCOs (F)         : 12"
         "Missing BUSCOs (M)            : 45"
         "Total BUSCO groups searched   : 1291"
         (then compute percentages)
    Returns dict with keys C, D, F, M, n (C/D/F/M are percentages).
    """
    text = path.read_text(encoding="utf-8", errors="replace")

    # --- Pattern 1: counts in the compact line ---
    # Example: C:1234[S:1200,D:34],F:12,M:45,n:1291
    m_counts = re.search(
        r"\bC:(\d+)\s*\[\s*S:(\d+)\s*,\s*D:(\d+)\s*\]\s*,\s*F:(\d+)\s*,\s*M:(\d+)\s*,\s*n:(\d+)\b",
        text,
    )
    if m_counts:
        Cc, Sc, Dc, Fc, Mc, n = map(int, m_counts.groups())
        # sanity: C counts should be S + D
        if Cc != (Sc + Dc):
            # don't fail; BUSCO formats can vary, but usually equals
            Cc = Sc + Dc
        return {
            "C": pct(Cc, n),
            "D": pct(Dc, n),
            "F": pct(Fc, n),
            "M": pct(Mc, n),
            "n": int(n),
        }

    # --- Pattern 2: percents in the compact line ---
    # Example: C:98.1%[S:95.0%,D:3.1%],F:0.9%,M:1.0%,n:1291
    m_pcts = re.search(
        r"\bC:([0-9]*\.?[0-9]+)%\s*\[\s*S:([0-9]*\.?[0-9]+)%\s*,\s*D:([0-9]*\.?[0-9]+)%\s*\]\s*,\s*F:([0-9]*\.?[0-9]+)%\s*,\s*M:([0-9]*\.?[0-9]+)%\s*,\s*n:(\d+)\b",
        text,
    )
    if m_pcts:
        C, S, D, F, M, n = m_pcts.groups()
        return {
            "C": round(float(C), 3),
            "D": round(float(D), 3),
            "F": round(float(F), 3),
            "M": round(float(M), 3),
            "n": int(n),
        }

    # --- Pattern 3: older tabular lines (counts) ---
    def find_int(label_regex: str) -> int | None:
        mm = re.search(label_regex, text, flags=re.IGNORECASE | re.MULTILINE)
        if not mm:
            return None
        return int(mm.group(1))

    Cc = find_int(r"Complete\s+BUSCOs\s*\(C\)\s*:\s*(\d+)")
    Dc = find_int(r"Complete\s+and\s+duplicated\s*\(D\)\s*:\s*(\d+)")
    Fc = find_int(r"Fragmented\s+BUSCOs\s*\(F\)\s*:\s*(\d+)")
    Mc = find_int(r"Missing\s+BUSCOs\s*\(M\)\s*:\s*(\d+)")
    n = find_int(r"Total\s+BUSCO\s+groups\s+searched\s*:\s*(\d+)")

    if all(v is not None for v in [Cc, Dc, Fc, Mc, n]):
        return {
