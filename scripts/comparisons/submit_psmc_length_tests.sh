#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
source_run="${1:?Supply completed full-preset comparison directory}"
out="${2:?Supply a new length-test output directory}"
hifi=$(realpath "${3:?Supply the HiFi comparison run with completed mapping}")
python3 "$repo/scripts/comparisons/prepare_psmc_length_tests.py" --source "$source_run" --out "$out" --hifi-run "$hifi"
out=$(realpath "$out")
mkdir -p "$out/logs"
dependencies=()
for label in asm20_20kb asm20_10kb; do
    job=$(sbatch --parsable --output="$out/logs/$label-%j.out" \
        "$repo/scripts/comparisons/run_psmc_full_preset.sbatch" "$out" "$label")
    job=${job%%;*};dependencies+=("$job")
    audit=$(sbatch --parsable --dependency="afterany:$job" --output="$out/logs/reliability-$label-%j.out" \
        "$repo/scripts/comparisons/run_psmc_reliability.sbatch" "$out" "$out/reliability_$label" \
        --left-run asm20_50kb --right-run "$label" --hifi-run "$hifi")
    audit=${audit%%;*};dependencies+=("$audit")
    printf '%s: workflow %s; reliability audit %s\n' "$label" "$job" "$audit"
done
dependency=$(IFS=:; printf '%s' "${dependencies[*]}")
collector=$(sbatch --parsable --dependency="afterany:$dependency" --output="$out/logs/summary-%j.out" \
    "$repo/scripts/comparisons/run_psmc_length_summary.sbatch" "$out")
printf 'Summary job: %s\nCombined review archive: %s.review.tar.gz\n' "$collector" "$out"
