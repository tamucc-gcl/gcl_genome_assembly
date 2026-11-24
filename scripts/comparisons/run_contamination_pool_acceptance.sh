#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=contam_pool_test
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --exclude=crest-c018,crest-c021,crest-c028,crest-c029
#SBATCH --output=comparisons/contamination-pool-acceptance-01/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
module load nextflow/23.10.1 miniconda3
pipeline="$PWD/gcl_genome_assembly"
run="$PWD/comparisons/contamination-pool-acceptance-01"
node="${POOL_TEST_NODE:-crest-c018}"
nextflow -log "$run/nextflow-${SLURM_JOB_ID}.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_pool_acceptance.nf" -profile slurm \
    -w "$run/work" --outdir "$run/results" --pool_test_root "$run/fixture" --pool_test_node "$node"
tar -czf "$run.review.tar.gz" -C "$run" results fixture/acceptance.json \
    fixture/first fixture/resume "nextflow-${SLURM_JOB_ID}.log" "slurm-${SLURM_JOB_ID}.out"
printf 'Review archive: %s\n' "$run.review.tar.gz"
