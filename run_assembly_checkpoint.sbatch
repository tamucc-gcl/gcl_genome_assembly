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
hic_readsets="${2:-}"
[[ "$hic_readsets" == NONE ]] && hic_readsets=''
hic_args=()
if [[ -n "$hic_readsets" ]]; then
    [[ "$hic_readsets" == /* ]] || hic_readsets="${PROJECT_ROOT}/${hic_readsets}"
    hic_args=(--hic_readsets "$hic_readsets")
fi
mode="${3:-validate}"
graph_enabled=false
graph_validate_only=false
case "$mode" in
    validate) launch_args=(-profile slurm -entry VALIDATE_INPUTS) ;;
    assembly) launch_args=(-profile slurm -resume) ;;
    graph_validate) launch_args=(-profile slurm -resume); graph_enabled=true; graph_validate_only=true ;;
    graph_build) launch_args=(-profile slurm -resume); graph_enabled=true ;;
    *) echo "Mode must be validate, assembly, graph_validate or graph_build" >&2; exit 2 ;;
esac
# Validation uses a separate launch directory so it cannot become the assembly resume target.
launch_dir="${PROJECT_ROOT}/.nf/assembly"
[[ "$mode" == validate ]] && launch_dir="${PROJECT_ROOT}/.nf/assembly-input-check"

outdir="${PROJECT_ROOT}/genome_assembly"
report_outdir="$outdir"
extra_arguments=("${@:4}")
for ((i=0; i<${#extra_arguments[@]}; i++)); do
    case "${extra_arguments[i]}" in
        --outdir) report_outdir="${extra_arguments[i+1]:?Missing --outdir value}" ;;
        --outdir=*) report_outdir="${extra_arguments[i]#--outdir=}" ;;
    esac
done
[[ "$report_outdir" == /* ]] || report_outdir="${PROJECT_ROOT}/$report_outdir"
outdir="$report_outdir"
# Nextflow launches below .nf/assembly; keep caller-supplied output paths relative
# to the project rather than accidentally publishing under that launch directory.
for ((i=0; i<${#extra_arguments[@]}; i++)); do
    case "${extra_arguments[i]}" in
        --outdir) extra_arguments[i+1]="$report_outdir" ;;
        --outdir=*) extra_arguments[i]="--outdir=$report_outdir" ;;
    esac
done

echo "Project root  : ${PROJECT_ROOT}"
echo "Sample sheet  : ${sample_sheet}"
echo "Output dir    : ${outdir}"


mkdir -p "$launch_dir" "${PROJECT_ROOT}/logs"
cd "$launch_dir"
nextflow -log "${PROJECT_ROOT}/logs/nextflow_${mode}_${SLURM_JOB_ID:-manual}.log" \
    run "${PROJECT_ROOT}/gcl_genome_assembly/main.nf" \
        "${launch_args[@]}" \
        -w "${PROJECT_ROOT}/work" \
        --publish_dir_mode 'symlink' \
        --run_pangenome "$graph_enabled" \
        --pangenome_validate_only "$graph_validate_only" \
        "${hic_args[@]}" \
        --sample_sheet "${sample_sheet}" \
        --outdir "${outdir}" \
        "${extra_arguments[@]}"

cd "${PROJECT_ROOT}"
if [[ "$mode" != validate ]]; then
    python3 "${PROJECT_ROOT}/gcl_genome_assembly/py_scripts/validate_report.py" --root "$report_outdir"
fi
#bash scripts/utils/fasta_to_webshare.sh
