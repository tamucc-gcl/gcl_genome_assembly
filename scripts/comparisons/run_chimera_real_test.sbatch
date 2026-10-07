#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_real
#SBATCH --array=0-1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=logs/chimera-real-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
tool_env=${1:?Supply analysis conda environment}
repo="${PWD}/gcl_genome_assembly"
assessment="${2:-${PWD}/comparisons/chimera-foundation-01}"
bam_root="${3:-${PWD}/genome_assembly/assembly/chimeras/sequence_context}"
scenarios=(auto manual)
scenario=${scenarios[$SLURM_ARRAY_TASK_ID]}
out="${PWD}/comparisons/chimera-real-${SLURM_ARRAY_JOB_ID}/${scenario}"
mkdir -p "$out"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    tar --dereference --exclude=work --exclude='*.fa' --exclude='*.fasta' \
        --exclude='*.bam' --exclude='*.bai' --exclude='*.csi' --exclude='*.sam' \
        -czf "$out.review.tar.gz" -C "$(dirname "$out")" "$scenario"
    exit "$rc"
}
trap finish EXIT
module load miniconda3
module load nextflow/23.10.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
export PYTHONDONTWRITEBYTECODE=1
required_tools=(python samtools minimap2 nextflow)
if [[ "$scenario" == auto ]]; then required_tools+=(tidk); fi
missing=()
for tool in "${required_tools[@]}"; do
    if ! command -v "$tool"; then missing+=("$tool"); fi
done
if (( ${#missing[@]} )); then
    printf 'ERROR: missing tools in environment %s: %s\n' "$tool_env" "${missing[*]}" >&2
    printf 'Missing tools: %s\n' "${missing[*]}" > "$out/preflight_error.txt"
    exit 1
fi
python -c 'import numpy, matplotlib, cooler'

python "$repo/scripts/comparisons/prepare_chimera_real_test.py" \
    --assessment "$assessment" --bam-root "$bam_root" --species 373251 --sample Sde-CTlk_104 \
    --reviewed-template "$repo/tests/fixtures/ctlk_reviewed_action.json" --out "$out/inputs"
mode=auto
context=true
evidence=true
if [[ "$scenario" == manual ]]; then
    mode="$out/inputs/reviewed-actions.tsv"
    context=false
    evidence=false
fi
cd "$out"
nextflow run "$repo/chimera_real_test.nf" \
    -c "$repo/nextflow.config" -c "$repo/scripts/comparisons/chimera_real_test.config" \
    --manifest "$out/inputs/manifest.json" --reviewed_actions_json "$out/inputs/reviewed-actions.json" \
    --scenario "$scenario" --chimera_break "$mode" --chimera_sequence_context "$context" \
    --chimera_evidence "$evidence" --chimera_hifi_bam_manifest "$out/inputs/bams.json" \
    --outdir "$out/results" -work-dir "$out/work" -with-trace "$out/trace.tsv" \
    -with-report "$out/execution.html"
test -s "$out/results/verification/Sde-CTlk_104_hap1.verification.json"
test -s "$out/results/verification/Sde-CTlk_104_hap2.verification.json"

if [[ "$scenario" == auto ]]; then
    for hap in hap1 hap2; do
        test -s "$out/results/assembly/chimeras/sequence_context/Sde-CTlk_104_${hap}.sequence_context/decision_measurements.json"
        test -s "$out/results/assembly/chimeras/sequence_context/Sde-CTlk_104_${hap}.sequence_context/chromosome_blocks.tsv"
        test -s "$out/results/assembly/chimeras/sequence_context/Sde-CTlk_104_${hap}.sequence_context/candidate_1.farther_contacts.png"
        test -s "$out/results/assembly/chimeras/decisions/Sde-CTlk_104_${hap}.decisions/actions.tsv"
    done
fi
