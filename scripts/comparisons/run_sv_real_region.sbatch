#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=sv_real_chromosome
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/sv-real-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
envdir="${1:?Supply environment, reference FASTA, query FASTA, reference chromosome, query chromosome}"
ref="${2:?Reference FASTA}"
query="${3:?Query FASTA}"
refchr="${4:?Reference chromosome ID}"
querychr="${5:?Query chromosome ID}"
anchors="${6:-maxmatch}"
case "$anchors" in maxmatch|mumreference) ;; *) echo 'Anchor mode must be maxmatch or mumreference' >&2; exit 2;; esac
anchor_args=()
if [[ "$anchors" = maxmatch ]]; then anchor_args=(--maxmatch); fi
[[ "$envdir" = /* ]] || { echo 'Environment must be absolute' >&2; exit 2; }
export PATH="$envdir/bin:$PATH"
for exe in samtools nucmer delta-filter show-coords show-snps syri; do
  [[ -x "$envdir/bin/$exe" ]] || { echo "Missing $exe" >&2; exit 2; }
done
[[ -x /usr/bin/time ]] || { echo 'GNU time is required' >&2; exit 2; }
ref="$(realpath "$ref")"
query="$(realpath "$query")"
# Require existing indexes; never create files in production result directories.
for source in "$ref" "$query"; do
  [[ -f "$source.fai" ]] || { echo "Missing existing index: $source.fai" >&2; exit 2; }
done
check_chromosome() {
  awk -F '\t' -v id="$2" '$1==id {n++; len=$2} END {if(n!=1 || len<=0 || len>100000000) exit 1}' "$1.fai"
}
check_chromosome "$ref" "$refchr" || { echo 'Reference ID missing/duplicated or exceeds 100 Mb pilot cap' >&2; exit 2; }
check_chromosome "$query" "$querychr" || { echo 'Query ID missing/duplicated or exceeds 100 Mb pilot cap' >&2; exit 2; }
out="comparisons/sv-real-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
printf '%s\n' "$anchors" > "$out/anchor_mode.txt"
cp "$0" "$out/submitted_runner.sh"
cp "$envdir"/conda-meta/*.json "$out/"
printf 'role\tsource\toriginal_id\tcomparison_id\nreference\t%s\t%s\tchr1\nquery\t%s\t%s\tchr1\n' \
  "$ref" "$refchr" "$query" "$querychr" > "$out/inputs.tsv"
# Trap retains diagnostic results even on a normal command failure.
completed=0
stage=extraction
finish() {
  rc=$?
  trap - EXIT TERM INT HUP
  # A terminated foreground process can leave a misleading shell status.
  # Success requires reaching the explicit end of the analysis as well.
  if [[ "$rc" -eq 0 && "$completed" -ne 1 ]]; then rc=1; fi
  printf '%s\n' "$rc" > "$out/exit_status.txt"
  printf 'stage\tcompleted\texit_status\n%s\t%s\t%s\n' "$stage" "$completed" "$rc" > "$out/run_status.tsv"
  tar --exclude='*.fa' --exclude='*.delta' -czf "${out}.tar.gz" -C "$(dirname "$out")" "$(basename "$out")"
  printf 'Return %s.tar.gz and logs/sv-real-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
  exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
trap 'exit 129' HUP
# Normalize names only in local copies of the selected, user-confirmed homologs.
samtools faidx "$ref" "$refchr" | sed '1s/.*/>chr1/' > "$out/reference.fa"
samtools faidx "$query" "$querychr" | sed '1s/.*/>chr1/' > "$out/query.fa"
cd "$out"
sha256sum reference.fa query.fa > selected_sequence_sha256.txt
exec > >(tee -a commands.log) 2>&1
set -x
stage=alignment
printf '%s\n' "$stage" > stage.txt
/usr/bin/time -v -o alignment.resources.txt nucmer --threads "${SLURM_CPUS_PER_TASK:-4}" "${anchor_args[@]}" -c 100 -b 500 -l 50 reference.fa query.fa
stage=filtering
printf '%s\n' "$stage" > stage.txt
show-coords -THrd out.delta > out.raw.coords
/usr/bin/time -v -o filtering.resources.txt delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
show-coords -THrd out.filtered.delta > out.filtered.coords
stage=syri
printf '%s\n' "$stage" > stage.txt
/usr/bin/time -v -o syri.resources.txt syri -c out.filtered.coords -d out.filtered.delta -r reference.fa -q query.fa
[[ -s syri.out && -s syri.summary && -s syri.vcf ]]
completed=1
stage=complete
printf '%s\n' "$stage" > stage.txt
