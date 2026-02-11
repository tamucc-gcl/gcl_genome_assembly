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
            "C": pct(float(Cc), float(n)),
            "D": pct(float(Dc), float(n)),
            "F": pct(float(Fc), float(n)),
            "M": pct(float(Mc), float(n)),
            "n": int(n),
        }

    raise ValueError(
        f"Could not parse BUSCO short_summary TXT: {path}\n"
        "Tried compact summary line and older tabular summary formats."
    )


def parse_busco_short_summary_json(path: Path) -> Dict[str, float]:
    """
    Parse BUSCO short_summary*.json (BUSCO v5).

    BUSCO's JSON keys have varied slightly across versions.
    We try to locate counts for:
      complete, duplicated, fragmented, missing, total
    Then convert to percentages.
    """
    obj = read_json(path)

    # Common structure examples (varies):
    # obj["results"]["complete_busco"] etc
    # obj["results"]["complete"] etc
    # obj["results"]["summary"]["complete"] etc
    # or nested under obj["results"]["lineage_dataset"]["..."]
    #
    # We'll search recursively for the first dict that contains the needed fields.
    wanted_sets = [
        ("complete", "duplicated", "fragmented", "missing", "total"),
        ("complete_busco", "duplicated_busco", "fragmented_busco", "missing_busco", "total_busco_groups_searched"),
        ("Complete", "Duplicated", "Fragmented", "Missing", "Total"),
    ]

    def iter_dicts(x: Any):
        if isinstance(x, dict):
            yield x
            for v in x.values():
                yield from iter_dicts(v)
        elif isinstance(x, list):
            for v in x:
                yield from iter_dicts(v)

    found = None
    found_keys = None

    for d in iter_dicts(obj):
        for ks in wanted_sets:
            if all(k in d for k in ks):
                found = d
                found_keys = ks
                break
        if found:
            break

    if not found:
        # Special case: BUSCO sometimes reports complete/single/dup separately
        # and has "n_markers"/"total" elsewhere. We'll try a second heuristic.
        # Look for a dict with complete/single-copy/duplicated/fragmented/missing,
        # plus total (often "n_markers" or similar).
        raise ValueError(
            f"Could not find BUSCO summary fields in JSON: {path}\n"
            "If you paste a small snippet of the JSON top-level keys, I can adapt the parser."
        )

    c_key, d_key, f_key, m_key, n_key = found_keys

    complete = float(found[c_key])
    duplicated = float(found[d_key])
    fragmented = float(found[f_key])
    missing = float(found[m_key])
    total = float(found[n_key])

    # Some BUSCO JSONs store percentages already; detect that:
    # If total is <= 100 and values look like percentages, treat as percent inputs.
    # Otherwise treat as counts.
    if total <= 100.0 and (complete <= 100.0 and duplicated <= 100.0 and fragmented <= 100.0 and missing <= 100.0):
        # Percent inputs; need an integer n somewhere else—try to find it.
        # If not found, set n to 0 and still inject percents.
        n = None
        # common places
        for d in iter_dicts(obj):
            for cand in ("n_markers", "total_buscos", "total_busco_groups_searched", "total", "n"):
                if cand in d and isinstance(d[cand], (int, float)) and float(d[cand]) > 100:
                    n = int(d[cand])
                    break
            if n is not None:
                break
        if n is None:
            n = 0
        return {
            "C": round(complete, 3),
            "D": round(duplicated, 3),
            "F": round(fragmented, 3),
            "M": round(missing, 3),
            "n": int(n),
        }

    # Count inputs
    n = int(round(total))
    return {
        "C": pct(complete, total),
        "D": pct(duplicated, total),
        "F": pct(fragmented, total),
        "M": pct(missing, total),
        "n": n,
    }


def parse_busco(path: Path) -> Dict[str, float]:
    if path.suffix.lower() == ".json":
        return parse_busco_short_summary_json(path)
    # .txt or anything else -> treat as text
    return parse_busco_short_summary_txt(path)


def inject_busco(assembly_stats_json: Dict[str, Any], busco: Dict[str, float]) -> Dict[str, Any]:
    # assembly-stats expects a top-level 'busco' key per docs
    assembly_stats_json["busco"] = {
        "C": float(busco["C"]),
        "D": float(busco["D"]),
        "F": float(busco["F"]),
        "M": float(busco["M"]),
        "n": int(busco["n"]),
    }
    return assembly_stats_json


def main() -> None:
    ap = argparse.ArgumentParser(description="Inject BUSCO summary into assembly-stats JSON.")
    ap.add_argument("--assembly-stats", required=True, type=Path, help="assembly-stats JSON from asm2stats.pl")
    ap.add_argument("--busco", required=True, type=Path, help="BUSCO short_summary*.txt or short_summary*.json")
    ap.add_argument("--out", required=True, type=Path, help="Output JSON path")
    args = ap.parse_args()

    asm = read_json(args.assembly_stats)
    busco = parse_busco(args.busco)
    asm2 = inject_busco(asm, busco)
    write_json(args.out, asm2)

    print(
        "Injected BUSCO into assembly-stats JSON:\n"
        f"  C={busco['C']}%  D={busco['D']}%  F={busco['F']}%  M={busco['M']}%  n={busco['n']}\n"
        f"Wrote: {args.out}"
    )


if __name__ == "__main__":
    main()
