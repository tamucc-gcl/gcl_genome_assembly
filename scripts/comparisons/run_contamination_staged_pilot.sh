#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=contam_staged
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=comparisons/contamination-staged-pilot-01/slurm-%j.out

# Submit from the Crest project root after creating the output directory.
# Only whole contigs totalling <=30 Mb are searched. No new DB download/build.
set -euo pipefail
printf "%s\n" "Heuristic sampling is retired. Use the standard FCS workflow; the full BlobTools branch remains optional. See gcl_genome_assembly/docs/contamination.md." >&2
exit 2

cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the project root}"
module load nextflow/23.10.1 miniconda3

project="$PWD"
pipeline="$project/gcl_genome_assembly"
run="$project/comparisons/contamination-staged-pilot-01"
protein="$project/comparisons/contamination-full-protein-01/db/diamond/uniprot_reference_proteomes/2026_04"
taxonomy="$project/comparisons/contamination-full-protein-01/db/ncbi_taxonomy/2026_04"
nt="$project/comparisons/contamination-ctlk-whole-02/db/blast/nt/2026-11-08/database/nt"

for required in "$pipeline/main.nf" "$protein/reference_proteomes.dmnd" \
    "$protein/protein_database.json" "$taxonomy/nodes.dmp" "$taxonomy/names.dmp" \
    "${nt}.nal" "$project/genome_assembly/pipeline/pipeline_trace.txt"; do
    if [[ ! -s "$required" ]]; then
        printf 'Missing required retained input: %s\n' "$required" >&2
        exit 1
    fi
done
mkdir -p "$run"
if [[ ! -s "$run/inputs/selection.json" ]]; then
    python3 "$pipeline/scripts/comparisons/prepare_contamination_real_pilot.py" contigs \
        --project "$project" --assembly Sde-CTlk_104_hap1 --max-bp 30000000 \
        --use-trace-coverage --out "$run/inputs"
fi
test -s "$run/inputs/contigs.fasta"
test -s "$run/inputs/coverage.tsv"

# Keep this pilot's Nextflow cache independent of the production run.
cd "$run"
resume=()
if [[ -d .nextflow/cache ]]; then
    resume=(-resume)
fi
nextflow -log "$run/nextflow.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_evidence.nf" \
    -profile slurm -w "$run/work" -qs 24 "${resume[@]}" \
    -with-trace "$run/trace.tsv" \
    --search_mode staged --stage post_fcs \
    --assembly Sde-CTlk_104_hap1 --taxid 373251 \
    --fasta "$run/inputs/contigs.fasta" --coverage_tsv "$run/inputs/coverage.tsv" \
    --protein_db_dir "$protein" --taxonomy_dir "$taxonomy" --nucleotide_db "$nt" \
    --contamination_megablast_execution chunks \
    --contamination_blastx_max_forks 8 --contamination_blastn_max_forks 8 \
    --targeted_blastn false --outdir "$run/results"

tar -czf "$run.review.tar.gz" --exclude='*.fasta' --exclude='*.bam' --exclude='*.bai' \
    -C "$run" results/review results/diagnostics inputs/selection.json nextflow.log trace.tsv
printf '\nReview archive: %s\n' "$run.review.tar.gz"
printf 'This is a subset pilot, not a production removal sheet. No BAM read-support panels were generated.\n'
