#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=ctlk_adjudicate
#SBATCH --array=0-7
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/ctlk-adjudicate-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Supply analysis environment with Python, samtools and minimap2}
packet=${2:-comparisons/junction-assessment-20261005-084623-1506302}
assessment=${3:-comparisons/chimera-foundation-01}
fusion=${4:-comparisons/ctlk-fusion-batch-1511216}
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
lanes=(anchors graph hifi hic)
hap=$(( SLURM_ARRAY_TASK_ID % 2 + 1 ))
lane=${lanes[$(( SLURM_ARRAY_TASK_ID / 2 ))]}
out="comparisons/ctlk-adjudication-${SLURM_ARRAY_JOB_ID}/${lane}_hap${hap}"
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
python -m unittest discover -s gcl_genome_assembly/tests -p test_ctlk_adjudication.py
python gcl_genome_assembly/scripts/comparisons/run_ctlk_adjudication.py \
    --lane "$lane" --hap "$hap" --packet "$packet" --assessment "$assessment" \
    --fusion "$fusion" --out "$out" --threads "${SLURM_CPUS_PER_TASK:-16}"
