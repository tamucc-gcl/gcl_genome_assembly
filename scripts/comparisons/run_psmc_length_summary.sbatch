#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_length_summary
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail
module load miniconda3
out=$(realpath "${1:?Supply prepared length-test directory}")
repo=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["repository"])' "$out/comparison.json")
python3 "$repo/scripts/comparisons/psmc_length_summary.py" --comparison "$out"
