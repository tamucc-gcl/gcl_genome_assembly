#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=contam_preflight
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=comparisons/contamination-staged-preflight-01/slurm-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the project root}"
module load miniconda3

python3 -u "$PWD/gcl_genome_assembly/scripts/comparisons/contamination_staged_preflight.py" \
    --project "$PWD" --out "$PWD/comparisons/contamination-staged-preflight-01" \
    --expected-assemblies 10 --workers "${SLURM_CPUS_PER_TASK:-4}"
