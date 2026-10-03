#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_support
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/chimera-support-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
results=${1:-genome_assembly}
samtools=${2:-comparisons/environments/syri-1.7.1/bin/samtools}
[[ -x "$samtools" ]] || { echo "Supply an existing samtools executable as argument 2: $samtools" >&2; exit 2; }
samtools=$(realpath "$samtools")
out="comparisons/chimera-support-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
python3 gcl_genome_assembly/scripts/comparisons/scan_chimera_support.py \
  "$results/assembly/chimeras/sequence_context" "$out" --samtools "$samtools"
cp gcl_genome_assembly/scripts/comparisons/scan_chimera_support.py "$out/"
"$samtools" --version > "$out/samtools.version.txt"
tar -czf "$out.tar.gz" "$out"
printf 'Return %s.tar.gz and logs/chimera-support-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
