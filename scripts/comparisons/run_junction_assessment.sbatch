#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=junction_assessment
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/junction-assessment-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Supply tool environment followed by diagnostic arguments}
shift
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
out="comparisons/junction-assessment-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
finish() {
    rc=$?
    trap - EXIT TERM HUP INT
    if [[ -d "$out" ]]; then
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        if ! tar --exclude='*.fa' --exclude='*.fastq' --exclude='*.sam' \
            -czf "$out.tar.gz" "$out"; then
            [[ "$rc" -ne 0 ]] || rc=1
        fi
        printf 'Return %s.tar.gz and logs/junction-assessment-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
    fi
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
python3 gcl_genome_assembly/scripts/comparisons/assess_junction_batch.py \
    --out "$out" --threads "${SLURM_CPUS_PER_TASK:-8}" "$@"
