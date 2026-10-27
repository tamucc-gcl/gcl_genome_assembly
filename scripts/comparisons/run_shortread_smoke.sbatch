#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=shortread_smoke
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/shortread-smoke-%j.out
set -euo pipefail
module load nextflow/23.10.1
module load miniconda3
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$(realpath gcl_genome_assembly)"
testdir="$(realpath "${1:?Supply prepared smoke-test directory}")"
mode="${2:-validate}"
[[ -s "$testdir/subset.json" && -s "$testdir/samples.csv" ]] || { echo 'Subset preparation incomplete' >&2; exit 2; }
case "$mode" in
  validate) launch=( -entry VALIDATE_INPUTS ); rundir="$testdir/validation" ;;
  assembly) launch=( -resume ); rundir="$testdir/launch" ;;
  *) echo 'Use validate or assembly' >&2; exit 2 ;;
esac
mkdir -p "$rundir" "$testdir/logs"
runlog="$testdir/logs/nextflow-${mode}-${SLURM_JOB_ID}.log"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$testdir/exit_status-${mode}.txt"
    files=(subset.json samples.csv "exit_status-${mode}.txt")
    for name in "$runlog" "$testdir/trace-${mode}.tsv" "$testdir/effective-test.config"; do
        [[ -f "$name" ]] && files+=("${name#"$testdir/"}")
    done
    # Review native tool outcomes without collecting reads, assemblies or graphs.
    if [[ -d "$testdir/results" ]]; then
        while IFS= read -r -d '' name; do files+=("${name#"$testdir/"}"); done < <(
            find -L "$testdir/results" -type f -size -20M \( \
                -name '*.tsv' -o -name '*.json' -o -name '*.md' -o \
                -name '*.log' -o -name '*.log.txt' -o -name '*.gb' -o \
                -name '*.gbk' -o -name '*.png' -o -name '*.fai' -o \
                -name '*.fcs_gx_report.txt' -o -name '*.taxonomy.rpt' \) -print0)
    fi
    tar --dereference -czf "$testdir.${mode}.review.tar.gz" -C "$testdir" "${files[@]}" || {
        echo 'Review archive creation failed' >&2
        [[ "$rc" -ne 0 ]] || rc=1
    }
    exit "$rc"
}
trap finish EXIT
cp "$repo/scripts/comparisons/shortread_smoke.config" "$testdir/effective-test.config"
cd "$rundir"
nextflow -log "$runlog" \
  run "$repo/main.nf" -profile slurm "${launch[@]}" \
  -c "$repo/scripts/comparisons/shortread_smoke.config" \
  -w "$testdir/work" --sample_sheet "$testdir/samples.csv" \
  --outdir "$testdir/results" --publish_dir_mode symlink \
  -with-trace "$testdir/trace-${mode}.tsv" \
  "${@:3}"
if [[ "$mode" == assembly ]]; then
    if [[ ! -s "$testdir/results/assembly_run_summary.md" ]]; then
        echo '[SHORTREAD_PILOT] Nextflow completed but the required assembly_run_summary.md is missing' >&2
        exit 1
    fi
fi
