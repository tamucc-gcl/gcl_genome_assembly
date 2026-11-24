#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=fcs_assembly_test
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=comparisons/fcs-assembly-acceptance-01/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
module load nextflow/23.10.1 miniconda3
pipeline="$PWD/gcl_genome_assembly"
run="$PWD/comparisons/${FCS_ACCEPTANCE_RUN:-fcs-assembly-acceptance-01}"
mkdir -p "$run"
helper="$pipeline/scripts/comparisons/prepare_fcs_assembly_acceptance.py"
python3 "$helper" configure --project "$PWD" --out "$run" "$@"
cat > "$run/acceptance.config" <<'CONFIG'
process {
    withLabel: fcs { maxForks = 1 }
}
CONFIG
image=()
if [[ -n "${FCS_ADAPTOR_IMAGE:-}" ]]; then
    image=(--fcs_adaptor_container "$FCS_ADAPTOR_IMAGE")
fi
# Keep the controller cache beside these results, independent of other runs.
cd "$run"
nextflow -log "$run/nextflow-${SLURM_JOB_ID}.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/fcs_assembly_acceptance.nf" -profile slurm \
    -c "$run/acceptance.config" -params-file "$run/params.json" \
    -w "$run/work" --outdir "$run/results" "${image[@]}" -resume
python3 "$helper" verify --out "$run"
tar -czf "$run.review.tar.gz" --exclude='*.fasta' --exclude='*.fa' \
    --exclude='*.bam' --exclude='*.bai' \
    -C "$run" results inputs.json params.json acceptance.json acceptance.config \
    "nextflow-${SLURM_JOB_ID}.log" "slurm-${SLURM_JOB_ID}.out"
printf 'Review archive: %s\n' "$run.review.tar.gz"
