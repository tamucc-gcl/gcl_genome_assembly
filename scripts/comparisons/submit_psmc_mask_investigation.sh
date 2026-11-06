#!/usr/bin/env bash
# Removable diagnostic: ordinary Bash/Slurm arrays, independent of Nextflow sessions.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
run="${1:?Supply the completed run directory (containing results and trace.tsv)}"
out="${2:?Supply a new investigation output directory}"
sample="${3:?Supply the individual for the separate complete-HiFi comparison}"
python3 "$repo/scripts/comparisons/psmc_mask_investigation.py" prepare \
    --run "$run" --out "$out" --hifi-sample "$sample" "${@:4}"
manifest=$(realpath "$out/manifest.json")
out=$(dirname "$manifest")
mkdir -p "$out/logs"
count=$(python3 -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["tasks"]))' "$manifest")
array=$(sbatch --parsable --array="0-$((count-1))" --output="$out/logs/haplotype-%A_%a.out" \
    "$repo/scripts/comparisons/run_psmc_mask_investigation.sbatch" "$manifest" task)
array=${array%%;*}
printf 'Submitted haplotype array: %s\n' "$array"
collector=$(sbatch --parsable --dependency="afterany:$array" --cpus-per-task=1 --mem=8G --time=2:00:00 \
    --output="$out/logs/report-%j.out" \
    "$repo/scripts/comparisons/run_psmc_mask_investigation.sbatch" "$manifest" report)
printf 'Report job: %s\nHiFi manifest: %s/hifi_manifest.tsv\nReview archive: %s.review.tar.gz\n' "$collector" "$out" "$out"
