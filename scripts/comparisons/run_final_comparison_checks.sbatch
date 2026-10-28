#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=final_comparisons
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/final-comparisons-%j.out
set -euo pipefail
module load nextflow/23.10.1
module load miniconda3
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$(realpath gcl_genome_assembly)"
manifest="$(realpath "${1:?Supply accepted assembly_eligibility.json}")"
out="${2:?Supply comparison output directory}"
[[ "$out" == /* ]] || out="$root/$out"
mkdir -p "$out/launch" "$root/logs"
entry="$repo/.native-final-comparison-checks.nf"
[[ ! -e "$entry" ]] || { echo 'Another final comparison check is using this entry' >&2; exit 2; }
finish() {
    rc=$?
    trap - EXIT
    rm -f "$entry"
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    files=(exit_status.txt source_assembly_manifest.json)
    for name in nextflow.log trace.tsv results; do
        [[ -e "$out/$name" ]] && files+=("$name")
    done
    tar --dereference --exclude='*.paf.gz' -czf "$out.review.tar.gz" -C "$out" "${files[@]}" || {
        [[ "$rc" -ne 0 ]] || rc=1
    }
    exit "$rc"
}
cp "$manifest" "$out/source_assembly_manifest.json"
trap finish EXIT
cp "$repo/tests/native_final_comparisons.nf" "$entry"
cd "$out/launch"
nextflow -log "$out/nextflow.log" run "$entry" -profile slurm -resume \
    -w "$out/work" --outdir "$out/results" --publish_dir_mode symlink \
    --assembly_manifest "$manifest" --pairwise_alignment_mode within_sample \
    -with-trace "$out/trace.tsv" "${@:3}"
test -s "$out/results/comparison_acceptance.md"
