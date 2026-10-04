#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_history
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/chimera-history-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Supply Conda environment then comparison arguments}
shift
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
out="comparisons/chimera-history-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
finish() {
    rc=$?
    trap - EXIT TERM INT HUP
    if [[ -d "$out" ]]; then
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        if ! tar --exclude='*.fa' -czf "$out.tar.gz" "$out"; then
            [[ "$rc" -ne 0 ]] || rc=1
        fi
        printf 'Return %s.tar.gz and logs/chimera-history-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
    fi
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
python3 gcl_genome_assembly/scripts/comparisons/compare_chimera_history.py \
    --out "$out" --threads "${SLURM_CPUS_PER_TASK:-4}" "$@"
