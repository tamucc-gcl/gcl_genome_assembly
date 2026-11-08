#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_reliability
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/psmc-reliability-%j.out
set -euo pipefail
module load miniconda3
comparison=$(realpath "${1:?Supply the completed full-preset comparison directory}")
out=$(realpath -m "${2:?Supply a new reliability output directory}")
if [[ -e "$out" ]]; then
    printf 'Use a new reliability output directory: %s\n' "$out" >&2
    exit 2
fi
repo=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["repository"])' "$comparison/comparison.json")
source_manifest=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source"]+"/manifest.json")' "$comparison/comparison.json")
prefix=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["controls"][0]["env"])' "$source_manifest")
plot_prefix=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["plot_env"])' "$source_manifest")
finish() {
    rc=$?
    trap - EXIT
    if [[ -d "$out" ]]; then
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        tar -czf "$out.review.tar.gz" -C "$out" . || { [[ "$rc" -ne 0 ]] || rc=1; }
    fi
    exit "$rc"
}
trap finish EXIT
export PATH="$prefix/bin:$PATH"
python3 "$repo/scripts/comparisons/psmc_reliability.py" --comparison "$comparison" --out "$out" "${@:3}"
"$plot_prefix/bin/Rscript" "$repo/r_scripts/plot_psmc_reliability.R" "$out"
