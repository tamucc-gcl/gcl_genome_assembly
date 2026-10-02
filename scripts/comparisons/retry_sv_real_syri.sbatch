#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=sv_syri_retry
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/sv-syri-retry-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
envdir="$(realpath "${1:?Supply environment and completed alignment directory}")"
previous="$(realpath "${2:?Supply completed alignment directory}")"
export PATH="$envdir/bin:$PATH"
for exe in syri show-snps; do
  [[ -x "$envdir/bin/$exe" ]] || { echo "Missing $exe" >&2; exit 2; }
done
for file in reference.fa query.fa out.filtered.delta out.filtered.coords; do
  [[ -s "$previous/$file" ]] || { echo "Missing or empty $previous/$file" >&2; exit 2; }
done
out="comparisons/sv-syri-retry-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
cp "$0" "$out/submitted_runner.sh"
cp "$envdir"/conda-meta/*.json "$out/"
printf '%s\n' "$previous" > "$out/alignment_source.txt"
for file in reference.fa query.fa out.filtered.delta out.filtered.coords; do
  cp "$previous/$file" "$out/$file"
done
for file in inputs.tsv anchor_mode.txt selected_sequence_sha256.txt alignment.resources.txt filtering.resources.txt; do
  if [[ -f "$previous/$file" ]]; then cp "$previous/$file" "$out/$file"; fi
done
completed=0
finish() {
  rc=$?
  trap - EXIT TERM INT HUP
  if [[ "$rc" -eq 0 && "$completed" -ne 1 ]]; then rc=1; fi
  printf '%s\n' "$rc" > "$out/exit_status.txt"
  tar --exclude='*.fa' --exclude='*.delta' -czf "${out}.tar.gz" -C "$(dirname "$out")" "$(basename "$out")"
  printf 'Return %s.tar.gz and logs/sv-syri-retry-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
  exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
trap 'exit 129' HUP
cd "$out"
sha256sum reference.fa query.fa out.filtered.delta out.filtered.coords > retry_input_sha256.txt
if [[ -f selected_sequence_sha256.txt ]]; then sha256sum -c selected_sequence_sha256.txt; fi
/usr/bin/time -v -o syri.resources.txt syri -c out.filtered.coords -d out.filtered.delta -r reference.fa -q query.fa > commands.log 2>&1
[[ -s syri.out && -s syri.summary && -s syri.vcf ]]
completed=1
