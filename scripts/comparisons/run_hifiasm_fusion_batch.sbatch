#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=ctlk_fusion_batch
#SBATCH --array=0-4%3
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/ctlk-fusion-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
tool_env=${1:?Supply existing environment with hifiasm 0.25.0-r726 and gfatools}
input_root=${2:-genome_assembly}
assessment=${3:-comparisons/chimera-foundation-01}
analysis_env=${4:-$tool_env}
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$analysis_env"
for tool in python3 samtools minimap2; do command -v "$tool" >/dev/null; done
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_hifiasm_fusion_batch.py'
out="comparisons/ctlk-fusion-batch-${SLURM_ARRAY_JOB_ID}"
mkdir -p "$out"
labels=(both A_only B_only hifi_only both_no_postjoin)
label=${labels[$SLURM_ARRAY_TASK_ID]}
finish() {
    rc=$?
    trap - EXIT TERM INT HUP
    if [[ -d "$out/$label" ]]; then
        printf '%s\n' "$rc" > "$out/$label/exit_status.txt"
        tar --exclude='*.bin' --exclude='*.gfa' --exclude='*.fa' \
            -czf "$out/$label.review.tar.gz" "$out/$label" || rc=1
    fi
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
python3 gcl_genome_assembly/scripts/comparisons/run_hifiasm_fusion_batch.py \
    --input-root "$input_root" --out "$out" --treatment "$SLURM_ARRAY_TASK_ID" \
    --assembler-env "$tool_env" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --reference "$assessment/assembly/scaffold/yahs_round2/Sde-CMat_203_hap1_round2_scaffolds.fa" \
    --reference-name-map "$assessment/assembly/harmonization/Sde-CMat_203_hap1.harmonized_name_map.tsv" \
    --baseline-hap1 "$assessment/assembly/contig/hifiasm/Sde-CTlk_104.hap1.p_ctg.fasta" \
    --baseline-hap2 "$assessment/assembly/contig/hifiasm/Sde-CTlk_104.hap2.p_ctg.fasta" \
    --target-contigs gcl_genome_assembly/planning_docs/investigations/ctlk-fusion-contigs.tsv
