#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=sv_diagnostics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/sv-diagnostics-%j.out
set -euo pipefail
export TZ=UTC
cd "${SLURM_SUBMIT_DIR:-$PWD}"
envdir="$(realpath "${1:?Supply R environment and SyRI result directory}")"
source_dir="$(realpath "${2:?Supply SyRI result directory}")"
[[ -x "$envdir/bin/Rscript" ]] || { echo 'Rscript missing' >&2; exit 2; }
for file in reference.fa query.fa syri.out out.filtered.coords; do
  [[ -s "$source_dir/$file" ]] || { echo "Missing $file" >&2; exit 2; }
done
script="$(realpath gcl_genome_assembly/r_scripts/sv_pilot_diagnostics.R)"
out="comparisons/sv-diagnostics-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
printf 'role\tbases\n' > "$out/sequence_lengths.tsv"
for role in reference query; do
  awk -v role="$role" '!/^>/ {gsub(/[[:space:]]/, ""); n+=length($0)} END {printf "%s\t%.0f\n",role,n}' \
    "$source_dir/$role.fa" >> "$out/sequence_lengths.tsv"
done
printf '%s\n' "$source_dir" > "$out/source.txt"
cp "$script" "$out/plot_script.R"
cp "$0" "$out/submitted_runner.sh"
"$envdir/bin/Rscript" "$script" "$source_dir" "$out" > "$out/plot.log" 2>&1 || {
  cat "$out/plot.log"; tar -czf "$out.tar.gz" "$out"; exit 1;
}
tar -czf "$out.tar.gz" "$out"
printf 'Return %s.tar.gz and logs/sv-diagnostics-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
