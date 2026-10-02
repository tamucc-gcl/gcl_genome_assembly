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
cd "$rundir"
nextflow -log "$testdir/logs/nextflow-${mode}-${SLURM_JOB_ID}.log" \
  run "$repo/main.nf" -profile slurm "${launch[@]}" \
  -c "$repo/scripts/comparisons/shortread_smoke.config" \
  -w "$testdir/work" --sample_sheet "$testdir/samples.csv" \
  --outdir "$testdir/results"
