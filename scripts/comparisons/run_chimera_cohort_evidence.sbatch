#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_cohort
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=logs/chimera-cohort-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
tool_env=${1:?Supply analysis conda environment}
repo="${PWD}/gcl_genome_assembly"
assessment="${2:-${PWD}/comparisons/chimera-foundation-01}"
bam_root="${3:-${PWD}/genome_assembly/assembly/chimeras/sequence_context}"
out="${PWD}/comparisons/chimera-cohort-${SLURM_JOB_ID}"
mkdir -p "$out"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    tar --dereference --exclude=work --exclude='*.fa' --exclude='*.fasta' \
        --exclude='*.bam' --exclude='*.bai' --exclude='*.csi' --exclude='*.sam' \
        -czf "$out.review.tar.gz" -C "$(dirname "$out")" "$(basename "$out")"
    exit "$rc"
}
trap finish EXIT
module load miniconda3
module load nextflow/23.10.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
export PYTHONDONTWRITEBYTECODE=1
for tool in python samtools minimap2 tidk nextflow; do command -v "$tool"; done
python -c 'import numpy, matplotlib, cooler'

# Existing assemblies are assessed; no historical cut fixture is supplied.
# Missing retained BAMs are visible in availability.tsv and candidate reports.
python "$repo/scripts/comparisons/prepare_chimera_real_test.py" \
    --assessment "$assessment" --bam-root "$bam_root" --species 373251 --sample all \
    --allow-missing-bams --out "$out/inputs"
cat "$out/inputs/evidence-availability.tsv"
cd "$out"
nextflow run "$repo/chimera_real_test.nf" \
    -c "$repo/nextflow.config" -c "$repo/scripts/comparisons/chimera_real_test.config" \
    --manifest "$out/inputs/manifest.json" --reviewed_actions_json "$out/inputs/reviewed-actions.json" \
    --scenario evidence --chimera_break false --chimera_sequence_context true \
    --chimera_evidence true --chimera_hifi_bam_manifest "$out/inputs/bams.json" \
    --outdir "$out/results" -work-dir "$out/work" -with-trace "$out/trace.tsv" \
    -with-report "$out/execution.html"

python - "$out" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1])
manifest=json.loads((root/'inputs/manifest.json').read_text())
expected={r['meta']['id'] for r in manifest['cohort']}
review=root/'results/assembly/chimeras/review'
with (review/'assembly_registry.tsv').open() as handle:
    registry=list(csv.DictReader(handle,delimiter='\t'))
assert {r['assembly'] for r in registry}==expected,'Cohort assessment incomplete'
assert all(r['assessment_sha256'] for r in registry),'Missing assessed FASTA identity'
for assembly in expected:
    assert (review/(assembly+'.review')/'report.md').is_file(),assembly
    check=json.loads((root/'results/verification'/(assembly+'.verification.json')).read_text())
    assert check['status']=='PASS',assembly
with (review/'chimera_review.tsv').open() as handle:
    assert all(r['selected']=='NO' for r in csv.DictReader(handle,delimiter='\t')),'Unexpected selected cut'
assert (review/'README.md').is_file()
print('Full-cohort evidence and unchanged assemblies verified:',len(expected))
PY
