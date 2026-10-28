#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=identity_reference
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/identity-reference-%j.out
set -euo pipefail
module load nextflow/23.10.1
module load miniconda3
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$(realpath gcl_genome_assembly)"
out="${1:?Supply comparison directory}"
[[ "$out" == /* ]] || out="$root/$out"
taxid="${2:?Supply taxid}"
species="${3:?Supply scientific species name}"
mkdir -p "$out/launch" "$root/logs"
# This removable diagnostic imports the actual production modules. Stage at the
# repository root so projectDir-based helper paths retain their normal meaning.
entry="$repo/.native-identity-reference-checks.nf"
[[ ! -e "$entry" ]] || { echo 'Another native identity check entry already exists' >&2; exit 2; }
finish() {
    rc=$?
    trap - EXIT
    rm -f "$entry"
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    files=(exit_status.txt)
    for name in nextflow.log trace.tsv results; do
        [[ -e "$out/$name" ]] && files+=("$name")
    done
    tar --dereference -czf "$out.review.tar.gz" -C "$out" "${files[@]}" || {
        [[ "$rc" -ne 0 ]] || rc=1
    }
    exit "$rc"
}
trap finish EXIT
cp "$repo/tests/native_identity_reference_checks.nf" "$entry"
cd "$out/launch"
nextflow -log "$out/nextflow.log" run "$entry" -profile slurm -resume \
    -w "$out/work" --outdir "$out/results" --publish_dir_mode symlink \
    --taxid "$taxid" --species "$species" -with-trace "$out/trace.tsv"
test -s "$out/results/identity_reference_acceptance.md"
