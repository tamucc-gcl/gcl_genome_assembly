#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=ctlk_boundaries
#SBATCH --array=0-5
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/ctlk-boundaries-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Supply analysis environment with Python, samtools and minimap2}
previous=${2:-comparisons/ctlk-adjudication-1511274}
assessment=${3:-comparisons/chimera-foundation-01}
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
tasks=(H01 J01 J02 J05 J06 J07)
task=${tasks[$SLURM_ARRAY_TASK_ID]}
out="comparisons/ctlk-boundaries-${SLURM_ARRAY_JOB_ID}/${task}"
finish() {
    rc=$?
    trap - EXIT TERM INT HUP
    if [[ -d "$out" ]]; then
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        tar --exclude='*.fa' --exclude='*.sam' --exclude='*.bam' \
            -czf "$out.review.tar.gz" "$out" || rc=1
    fi
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
python -m unittest discover -s gcl_genome_assembly/tests -p test_ctlk_boundary_batch.py
python gcl_genome_assembly/scripts/comparisons/run_ctlk_boundary_batch.py \
    --task "$task" --previous "$previous" --assessment "$assessment" \
    --out "$out" --threads "${SLURM_CPUS_PER_TASK:-16}"
