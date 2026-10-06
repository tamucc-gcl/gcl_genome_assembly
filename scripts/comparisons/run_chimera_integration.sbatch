#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_integration
#SBATCH --array=0-2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/chimera-integration-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
tool_env=${1:?Supply analysis conda environment}
repo="${PWD}/gcl_genome_assembly"
module load miniconda3
module load nextflow/23.10.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
command -v nextflow
export PYTHONDONTWRITEBYTECODE=1
scenarios=(unresolved eligible manual)
scenario=${scenarios[$SLURM_ARRAY_TASK_ID]}
out="${PWD}/comparisons/chimera-integration-${SLURM_ARRAY_JOB_ID}/${scenario}"
mkdir -p "$out"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    tar --dereference --exclude=work -czf "$out.review.tar.gz" -C "$(dirname "$out")" "$scenario"
    exit "$rc"
}
trap finish EXIT
python "$repo/tests/prepare_chimera_integration.py" --out "$out/fixture" --scenario "$scenario"
mode=auto
if [[ "$scenario" == manual ]]; then mode=file; fi
cd "$out"
nextflow run "$repo/tests/chimera_integration.nf" \
    --fixture "$out/fixture" --scenario "$scenario" --chimera_break "$mode" \
    --chimera_min_piece_bp 4 --publish_dir_mode copy \
    --outdir "$out/results" -work-dir "$out/work" -with-trace "$out/trace.tsv" \
    -with-report "$out/execution.html"
test -s "$out/results/verification/PASS.json"
