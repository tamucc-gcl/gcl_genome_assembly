#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=assembly_handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
set -euo pipefail
project="${SLURM_SUBMIT_DIR:-$PWD}"
bash "$project/gcl_genome_assembly/run_assembly_checkpoint.sbatch" \
    "${1:-data/assembly_samplesheet.csv}" "${2:-NONE}" assembly "${@:3}"
