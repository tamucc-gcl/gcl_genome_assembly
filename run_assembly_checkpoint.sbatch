#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=assembly_checkpoint
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/assembly-%j.out

#git clone git@github.com:tamucc-gcl/gcl_genome_assembly.git

set -euo pipefail

module load nextflow/23.10.1
module load miniconda3

# Capture the project root (where sbatch was invoked)
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# Resolve sample sheet to absolute path
sample_sheet="${1:-data/assembly_samplesheet.csv}"
[[ "$sample_sheet" != /* ]] && sample_sheet="${PROJECT_ROOT}/${sample_sheet}"

# The first argument remains the main sample sheet; the second is the optional Hi-C table.
hic_readsets="${2:-data/hic_readsets.csv}"
[[ "$hic_readsets" != /* ]] && hic_readsets="${PROJECT_ROOT}/${hic_readsets}"
mode="${3:-validate}"
case "$mode" in
    validate) launch_args=(-profile slurm -entry VALIDATE_INPUTS) ;;
    assembly) launch_args=(-profile slurm -resume) ;;
    *) echo "Mode must be validate or assembly" >&2; exit 2 ;;
esac
# Validation uses a separate launch directory so it cannot become the assembly resume target.
launch_dir="${PROJECT_ROOT}/.nf/assembly"
[[ "$mode" == validate ]] && launch_dir="${PROJECT_ROOT}/.nf/assembly-input-check"

outdir="${PROJECT_ROOT}/genome_assembly"

echo "Project root  : ${PROJECT_ROOT}"
echo "Sample sheet  : ${sample_sheet}"
echo "Output dir    : ${outdir}"


mkdir -p "$launch_dir" "${PROJECT_ROOT}/logs"
cd "$launch_dir"
nextflow -log "${PROJECT_ROOT}/logs/nextflow_${mode}_${SLURM_JOB_ID:-manual}.log" \
    run "${PROJECT_ROOT}/gcl_genome_assembly/main.nf" \
        "${launch_args[@]}" \
        -w "${PROJECT_ROOT}/work" \
        --publish_dir_mode 'link' \
        --taxid 373251 \
        --ploidy 2 \
        --qc_mode 'none' \
        --run_purge_dups true \
        --run_hic_balance false \
        --run_blobtools_evidence false \
        --yahs_round1_enzyme 'GATC,CTNAG,TTAA,GANTC' \
        --yahs_round2_enzyme 'GATC,CTNAG,TTAA,GANTC' \
        --inspector_contig_min_contig_bp 25000 \
        --inspector_contig_assemblyerror_min_bp 200000 \
        --inspector_contig_min_assembly_error_bp 5000 \
        --inspector_contig_max_assembly_error_bp 50000000 \
        --inspector_scaffold_min_contig_bp 25000 \
        --inspector_scaffold_assemblyerror_min_bp 2000000 \
        --inspector_scaffold_min_assembly_error_bp 25000 \
        --inspector_scaffold_max_assembly_error_bp 200000000 \
        --pairwise_alignment_mode 'all' \
        --chimera_break false \
        --run_pangenome false \
        --run_post_assembly false \
        --pangenome_2d_viz false \
        --harmonize_graph_rule 'plurality' \
        --hic_readsets "${hic_readsets}" \
        --sample_sheet "${sample_sheet}" \
        --outdir "${outdir}"

cd "${PROJECT_ROOT}"
#bash scripts/utils/fasta_to_webshare.sh
