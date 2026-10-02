#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=sv_benchmark
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/sv-benchmark-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
# One explicitly selected environment per submission; no installation in jobs.
tool="${1:?Usage: sbatch run_sv_benchmark.sbatch syri|syri-mummer|svim-asm /absolute/conda/environment}"
envdir="${2:?Supply an existing tool environment}"
suite="${3:-basic}"
anchors="${4:-maxmatch}"
case "$anchors" in maxmatch|mumreference) ;; *) echo 'Anchor mode must be maxmatch or mumreference' >&2; exit 2;; esac
# MUMmer 4 uses reference-unique anchors by default. It does not accept
# the older --mumreference spelling; retain that name only as our mode label.
anchor_args=()
if [[ "$anchors" = maxmatch ]]; then anchor_args=(--maxmatch); fi
if [[ "$tool" != syri-mummer && "$anchors" != maxmatch ]]; then
  echo 'Anchor selection is only supported for syri-mummer' >&2; exit 2
fi
case "$suite" in basic|repeats) ;; *) echo 'Suite must be basic or repeats' >&2; exit 2;; esac
case "$tool" in syri|syri-mummer|svim-asm) ;; *) echo 'Unknown tool' >&2; exit 2;; esac
[[ "$envdir" = /* ]] || { echo 'Environment must be an absolute path' >&2; exit 2; }
export PATH="$envdir/bin:$PATH"
caller="$tool"
extra=()
if [[ "$tool" = syri-mummer ]]; then
  caller=syri
  extra=(nucmer delta-filter show-coords show-snps)
fi
for exe in python minimap2 samtools "$caller" "${extra[@]}"; do
  [[ -x "$envdir/bin/$exe" ]] || { echo "Missing $envdir/bin/$exe" >&2; exit 2; }
done
[[ -x /usr/bin/time ]] || { echo 'GNU time is required' >&2; exit 2; }
out="comparisons/sv-benchmark-${tool}-${suite}-${anchors}-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
printf 'tool\tsuite\tanchor_mode\n%s\t%s\t%s\n' "$tool" "$suite" "$anchors" > "$out/settings.tsv"
python gcl_genome_assembly/scripts/comparisons/make_sv_benchmark.py --suite "$suite" --out "$out/fixtures"
cp "$0" "$out/submitted_runner.sh"
cp gcl_genome_assembly/scripts/comparisons/make_sv_benchmark.py "$out/fixture_generator.py"
cp "$envdir/conda-meta/history" "$out/environment_history.txt"
cp "$envdir"/conda-meta/*.json "$out/" # Exact installed builds; small metadata only.
minimap2 --version > "$out/minimap2_version.txt"
samtools --version > "$out/samtools_version.txt"
if [[ "$tool" = syri-mummer ]]; then nucmer --version > "$out/nucmer_version.txt" 2>&1; fi
printf 'case\ttool\tstatus\tseconds\n' > "$out/status.tsv"
failed=0
for name in control inversion tandem_duplication inverted_duplication translocation fragmented_control; do
  if [[ "$caller" = syri && "$name" = fragmented_control ]]; then
    printf '%s\t%s\tNOT_ELIGIBLE_chromosome_correspondence\t0\n' "$name" "$tool" >> "$out/status.tsv"
    continue
  fi
  mkdir "$out/$name"
  started=$SECONDS
  # Run each case in a subshell so failed tools retain diagnostics and other cases run.
  set +e
  (
    set -euo pipefail
    cd "$out/$name"
    set -x
    ref="$out/fixtures/reference.fa"
    query="$out/fixtures/$name.fa"
    if [[ "$caller" = syri ]]; then
      minimap2 -t "${SLURM_CPUS_PER_TASK:-4}" -ax asm5 --eqx "$ref" "$query" > alignment.sam
      # All alignment records, including secondary/supplementary; no bases/qualities.
      printf 'query\tflag\treference\tposition_1based\tmapq\tcigar\n' > minimap2_alignments.tsv
      awk 'BEGIN {OFS="\t"} !/^@/ {print $1,$2,$3,$4,$5,$6}' alignment.sam >> minimap2_alignments.tsv
      if [[ "$tool" = syri-mummer ]]; then
        # SyRI's documented MUMmer workflow; preserve raw and filtered blocks.
        /usr/bin/time -v -o alignment.resources.txt nucmer "${anchor_args[@]}" -c 100 -b 500 -l 50 "$ref" "$query"
        show-coords -THrd out.delta > out.raw.coords
        delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
        show-coords -THrd out.filtered.delta > out.filtered.coords
        /usr/bin/time -v -o syri.resources.txt syri -c out.filtered.coords -d out.filtered.delta -r "$ref" -q "$query"
      else
        syri -c alignment.sam -r "$ref" -q "$query" -k -F S
      fi
    else
      minimap2 -t "${SLURM_CPUS_PER_TASK:-4}" -ax asm5 --cs -r2k "$ref" "$query" |
        samtools sort -@ 1 -m 1G -o alignment.bam -
      samtools index alignment.bam
      svim-asm haploid calls alignment.bam "$ref"
    fi
  ) > "$out/$name/commands.log" 2>&1
  rc=$?
  set -e
  printf '%s\t%s\t%s\t%s\n' "$name" "$tool" "$rc" "$((SECONDS - started))" >> "$out/status.tsv"
  if [[ "$rc" -ne 0 ]]; then failed=1; fi
done
# Keep alignments/fixtures locally, return caller outputs and provenance only.
tar --exclude='*.sam' --exclude='*.bam' --exclude='*.bai' --exclude='*.fa' \
  -czf "${out}.tar.gz" -C "$(dirname "$out")" "$(basename "$out")"
printf 'Return %s.tar.gz and logs/sv-benchmark-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
exit "$failed"
