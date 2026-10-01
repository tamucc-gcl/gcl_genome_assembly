#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=gref_inspection
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/gref-inspection-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
# Activate the environment from the successful variant-audit task before submitting,
# or pass the absolute path to its bcftools executable with --bcftools.
out="comparisons/gref-inspection-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
python3 gcl_genome_assembly/scripts/comparisons/inspect_gref_catalog.py "$@" --out "$out"
tar -czf "${out}.tar.gz" "$out"
printf 'Return %s.tar.gz and this Slurm log.\n' "$out"
