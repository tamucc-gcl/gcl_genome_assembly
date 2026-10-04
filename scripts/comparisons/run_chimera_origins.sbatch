#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_origins
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/chimera-origins-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Usage: sbatch run_chimera_origins.sbatch TOOL_ENV RESULTS [extra diagnostic arguments]}
results=${2:?Supply comparison results directory}
shift 2
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
for tool in python3 samtools minimap2; do command -v "$tool" >/dev/null; done
out="comparisons/chimera-origins-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
finish() {
    rc=$?
    trap - EXIT TERM INT HUP
    if [[ -d "$out" ]]; then
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        # Never dereference stage FASTA symlinks or include large sequence inputs.
        if ! tar --exclude='*.fa' -czf "$out.tar.gz" "$out"; then
            echo 'Archive creation failed; outputs remain in the comparison directory.' >&2
            [[ "$rc" -ne 0 ]] || rc=1
        fi
        printf 'Return %s.tar.gz and logs/chimera-origins-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
    fi
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
python3 gcl_genome_assembly/scripts/comparisons/trace_chimera_origins.py \
    --results "$results" --work-root work --out "$out" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" "$@"
cp gcl_genome_assembly/scripts/comparisons/trace_chimera_origins.py "$out/"
cp gcl_genome_assembly/py_scripts/chimera_origin.py "$out/"
