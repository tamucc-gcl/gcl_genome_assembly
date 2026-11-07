#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
source_run="${1:?Supply the completed preset-mask comparison directory}"
out="${2:?Supply a new full-comparison output directory}"
python3 "$repo/scripts/comparisons/prepare_psmc_full_presets.py" --source "$source_run" --out "$out"
out=$(realpath "$out")
mkdir -p "$out/logs"
for preset in asm10 asm20; do
    job=$(sbatch --parsable --output="$out/logs/$preset-%j.out" \
        "$repo/scripts/comparisons/run_psmc_full_preset.sbatch" "$out" "$preset")
    printf '%s workflow: %s\nReview archive: %s/%s.review.tar.gz\n' "$preset" "$job" "$out" "$preset"
done
