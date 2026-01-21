#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=assembly_handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

#git clone git@github.com:tamucc-gcl/gcl_genome_assembly.git

module load nextflow/23.10.1
module load miniconda3

sample_sheet=data/samplesheet.csv
outdir=test_nf


nextflow run gcl_genome_assembly/main.nf \
    -profile slurm \
    -resume \
    --gxdb.profile test-only \
    --decon.source_taxid 373251 \
    --sample_sheet ${sample_sheet} \
    --outdir ${outdir}
