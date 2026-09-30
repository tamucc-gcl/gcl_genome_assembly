#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_region
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/chimera-region-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
stamp=$(date +%Y%m%d-%H%M%S)
out="comparisons/chimera-region-${stamp}-${SLURM_JOB_ID:-manual}"
python3 gcl_genome_assembly/scripts/comparisons/explore_chimera_region.py "$@" --out "$out"
tar -czf "${out}.tar.gz" "$out"
printf '\nDownload %s.tar.gz and this Slurm log.\n' "$out"
