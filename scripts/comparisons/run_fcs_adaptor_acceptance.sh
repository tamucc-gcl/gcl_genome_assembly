#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=fcs_adaptor_test
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=comparisons/fcs-adaptor-acceptance-01/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
module load nextflow/23.10.1 miniconda3
pipeline="$PWD/gcl_genome_assembly"
run="$PWD/comparisons/fcs-adaptor-acceptance-01"
nextflow -log "$run/nextflow-${SLURM_JOB_ID}.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/fcs_adaptor_acceptance.nf" -profile slurm \
    -w "$run/work" --outdir "$run/results" -resume
tar -czf "$run.review.tar.gz" --exclude='*.fasta' --exclude='*.fa' --exclude='*.fa.gz' \
    -C "$run" results "nextflow-${SLURM_JOB_ID}.log" "slurm-${SLURM_JOB_ID}.out"
printf 'Review archive: %s\n' "$run.review.tar.gz"
