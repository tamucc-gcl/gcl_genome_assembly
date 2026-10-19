#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=join_investigation
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/join-investigation-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$root/gcl_genome_assembly"
tool_env="${1:-ctlk_analysis}"
source_packet="${2:-$root/comparisons/chimera-observable-contact-v7/results/assembly/chimeras}"
final_dir="${3:-$root/comparisons/core-graph-assignment-01/results/assembly/final}"
out="${4:-$root/comparisons/four-join-investigation-01}"
for key in source_packet final_dir out; do
    [[ "${!key}" == /* ]] || printf -v "$key" '%s/%s' "$root" "${!key}"
done
mkdir -p "$out" "$root/logs"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    files=(exit_status.txt)
    for name in review trace.tsv nextflow.log inputs/manifest.json inputs/input_resolution.json preflight.log; do
        [[ -e "$out/$name" ]] && files+=("$name")
    done
    tar --dereference --exclude='*.sam' --exclude='*.bam' --exclude='*.fa' \
        -czf "$out.review.tar.gz" -C "$out" "${files[@]}"
    exit "$rc"
}
trap finish EXIT
module load miniconda3
module load nextflow/23.10.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
export PYTHONDONTWRITEBYTECODE=1
for tool in python3 samtools minimap2 nextflow; do command -v "$tool"; done
python3 -c 'import numpy, matplotlib'
python3 -m unittest discover -s "$repo/tests" -p test_join_investigation.py

# Validate retained sources on each launch. Unchanged manifest files retain their
# timestamps so ordinary Nextflow resume can reuse completed diagnostic tasks.
inputs="$out/inputs"
python3 "$repo/scripts/comparisons/join_investigation.py" prepare \
    --source "$source_packet" --final-dir "$final_dir" \
    --focus "$repo/scripts/comparisons/ctlk_four_join_focus.tsv" --out "$inputs" \
    2>&1 | tee "$out/preflight.log"
cd "$out"
nextflow -log "$out/nextflow.log" run "$repo/scripts/comparisons/join_investigation.nf" \
    -c "$repo/scripts/comparisons/join_investigation.config" \
    -resume -work-dir "$out/work" \
    --repo "$repo" --manifest "$inputs/manifest.json" --inputs "$inputs" \
    --outdir "$out" -with-trace "$out/trace.tsv" -with-report "$out/execution.html"
test -s "$out/review/README.md"
printf '\nReview: %s/review/README.md\nReturn: %s.review.tar.gz\n' "$out" "$out"
