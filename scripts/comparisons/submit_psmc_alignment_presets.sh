#!/usr/bin/env bash
# Removable preset-only diagnostic; no production workflow settings are changed.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
source_run="${1:?Supply the completed mask-investigation directory}"
out="${2:?Supply a new preset-comparison output directory}"
sample="${3:?Supply the individual to compare}"
python3 "$repo/scripts/comparisons/psmc_mask_investigation.py" prepare-presets \
    --source "$source_run" --out "$out" --sample "$sample" "${@:4}"
manifest=$(realpath "$out/manifest.json")
out=$(dirname "$manifest")
mkdir -p "$out/logs"
count=$(python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["tasks"]))' "$manifest")
array=$(sbatch --parsable --array="0-$((count-1))" --output="$out/logs/preset-%A_%a.out" \
    "$repo/scripts/comparisons/run_psmc_mask_investigation.sbatch" "$manifest" task)
array=${array%%;*}
collector=$(sbatch --parsable --dependency="afterany:$array" --cpus-per-task=1 --mem=8G --time=2:00:00 \
    --output="$out/logs/report-%j.out" \
    "$repo/scripts/comparisons/run_psmc_mask_investigation.sbatch" "$manifest" report-presets)
printf 'Preset array: %s\nReport job: %s\nReview archive: %s.review.tar.gz\n' "$array" "$collector" "$out"
