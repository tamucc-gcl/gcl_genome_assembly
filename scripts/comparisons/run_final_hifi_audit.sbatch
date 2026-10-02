#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=final_hifi_audit
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/final-hifi-audit-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
envdir="$(realpath "${1:?Supply existing environment containing samtools}")"
results="$(realpath "${2:-genome_assembly}")"
out="comparisons/final-hifi-audit-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
cp "$0" "$out/"
cp gcl_genome_assembly/scripts/comparisons/audit_final_hifi.py "$out/"
cp "$results/assembly_eligibility.json" "$out/"
rc=0
python3 gcl_genome_assembly/scripts/comparisons/audit_final_hifi.py "$results" "$out" \
  --samtools "$envdir/bin/samtools" > "$out/audit.log" 2>&1 || rc=$?
cat "$out/audit.log"
tar -czf "$out.tar.gz" "$out"
printf 'Return %s.tar.gz and this job log\n' "$out"
exit "$rc"
