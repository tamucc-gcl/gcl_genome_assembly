#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=contam_haplotype
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=crest-c018
#SBATCH --output=comparisons/contamination-staged-haplotype-01/slurm-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the project root}"
module load nextflow/23.10.1 miniconda3

project="$PWD"
pipeline="$project/gcl_genome_assembly"
run="$project/comparisons/contamination-staged-haplotype-01"
preflight="$project/comparisons/contamination-staged-preflight-01"
helper="$pipeline/scripts/comparisons/prepare_contamination_staged_haplotype.py"
mkdir -p "$run"
python3 "$helper" configure --project "$project" --preflight "$preflight" --out "$run"

# c018 has previously verified 1.8 TB local scratch. Keep the controller off it;
# otherwise a controller placed there could block its own exclusive worker.
cat > "$run/worker.config" <<'CONFIG'
process {
    withName: '.*CONTAMINATION_MEGABLAST_NODE' {
        clusterOptions = '--exclusive --mem=0 --nodelist=crest-c018'
    }
}
CONFIG

cd "$run"
resume=()
if [[ -d .nextflow/cache ]]; then
    resume=(-resume)
fi
nextflow -log "$run/nextflow.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/scripts/comparisons/contamination_evidence.nf" -profile slurm \
    -c "$run/worker.config" -params-file "$run/params.json" \
    -w "$run/work" -qs 32 "${resume[@]}" -with-trace "$run/trace.tsv" \
    --outdir "$run/results"

python3 "$helper" verify --expected "$run/expected.json" --evidence "$run/results/review/evidence.json"
tar -czf "$run.review.tar.gz" --exclude='*.fasta' --exclude='*.bam' --exclude='*.bai' \
    -C "$run" results/review results/diagnostics \
    results/pipeline/contamination_megablast_completion.json \
    params.json expected.json worker.config nextflow.log trace.tsv
printf '\nReview archive: %s\n' "$run.review.tar.gz"
