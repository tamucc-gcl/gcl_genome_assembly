#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_input_support
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
set -euo pipefail
module load miniconda3
root=$(realpath "${1:?Supply full pipeline test directory}")
sample="${2:?Supply individual ID}"
repo="${3:?Supply repository path}"
out="$root/psmc_checks/support_$sample"
finish() { rc=$?; trap - EXIT; mkdir -p "$out"; printf '%s\n' "$rc" > "$out/exit_status.txt"; exit "$rc"; }
trap finish EXIT
prefix=$(python3 "$repo/scripts/comparisons/psmc_cohort.py" environment --root "$root" --sample "$sample")
export PATH="$prefix/bin:$PATH"
python3 "$repo/scripts/comparisons/psmc_input_support.py" --results "$root/results" --sample "$sample" --out "$out"
