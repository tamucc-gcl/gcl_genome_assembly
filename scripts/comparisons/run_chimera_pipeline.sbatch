#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_pipeline
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/chimera-pipeline-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$root/gcl_genome_assembly"
sample_sheet="${1:-data/assembly_samplesheet.csv}"
hic_readsets="${2:-data/hic_readsets.csv}"
out="${3:-$root/comparisons/chimera-pipeline-${SLURM_JOB_ID}}"
[[ "$out" == /* ]] || out="$root/$out"
mkdir -p "$out" "$root/logs"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    if [[ -f "$root/logs/nextflow_assembly_${SLURM_JOB_ID}.log" ]]; then
        cp "$root/logs/nextflow_assembly_${SLURM_JOB_ID}.log" "$out/nextflow.log"
    fi
    # Package reports and evidence, not the full assembly output tree.
    files=(exit_status.txt)
    for name in nextflow.log trace.tsv results/assembly_report.md results/assembly/chimeras; do
        [[ -e "$out/$name" ]] && files+=("$name")
    done
    tar --dereference --exclude='*.fa' --exclude='*.fasta' --exclude='*.fa.gz' --exclude='*.fasta.gz' \
        --exclude='*.bam' --exclude='*.bai' --exclude='*.csi' --exclude='*.sam' \
        -czf "$out.review.tar.gz" -C "$out" "${files[@]}"
    exit "$rc"
}
trap finish EXIT

# This is the normal main.nf launch, with its full sample/read channels, process
# environments and SLURM resource requests. Reuse the production launch/cache
# and work paths; publish this test's outputs in a separate results directory.
bash "$repo/run_assembly_checkpoint.sbatch" "$sample_sheet" "$hic_readsets" assembly \
    --outdir "$out/results" --publish_dir_mode copy \
    --chimera_detect true --chimera_break false \
    --chimera_sequence_context true --chimera_evidence true --chimera_hifi_context true \
    -with-trace "$out/trace.tsv" -with-report "$out/execution.html" \
    "${@:4}"

test -s "$out/results/assembly/chimeras/review/README.md"
test -s "$out/results/assembly/chimeras/review/assembly_registry.tsv"
test -s "$out/results/assembly/chimeras/review/chimera_review.tsv"
