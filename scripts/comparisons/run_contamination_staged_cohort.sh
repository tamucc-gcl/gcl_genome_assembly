#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=contam_cohort
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=crest-c018,crest-c021,crest-c028,crest-c029
#SBATCH --output=comparisons/contamination-staged-cohort-01/slurm-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the project root}"
module load nextflow/23.10.1 miniconda3
project="$PWD"
pipeline="$project/gcl_genome_assembly"
run="$project/comparisons/contamination-staged-cohort-01"
session=747f1a27-488b-4b45-b991-201d684489ed
helper="$pipeline/scripts/comparisons/prepare_contamination_staged_cohort.py"

# Refuse a fresh assembly run if the previous native cache has disappeared.
launch=''
for candidate in "$project/.nf/assembly" "$project" "$pipeline"; do
    if [[ -d "$candidate/.nextflow/cache/$session" ]]; then
        launch="$candidate"
        break
    fi
done
if [[ -z "$launch" ]]; then
    printf 'Previous assembly cache %s was not found. No pipeline work was launched.\n' "$session" >&2
    exit 1
fi
mkdir -p "$run"
python3 "$helper" configure --project "$project" \
    --preflight "$project/comparisons/contamination-staged-preflight-01" --out "$run"

cat > "$run/workers.config" <<'CONFIG'
process {
    withName: '.*CONTAMINATION_MEGABLAST_COLLECT' {
        clusterOptions = {
            def controller = System.getenv('SLURMD_NODENAME') ?: System.getenv('HOSTNAME')
            if (!controller) throw new IllegalStateException('Controller hostname unavailable')
            "--nodelist=${controller}"
        }
    }
    withName: '.*CONTAMINATION_MEGABLAST_NODE' {
        clusterOptions = {
            def hosts = ['crest-c018', 'crest-c021', 'crest-c028', 'crest-c029']
            "--exclusive --mem=0 --nodelist=${hosts[(worker as int) - 1]}"
        }
    }
}
CONFIG

cd "$launch"
nextflow -log "$run/nextflow-${SLURM_JOB_ID}.log" -c "$pipeline/nextflow.config" run \
    "$pipeline/main.nf" -profile slurm -resume "$session" -w "$project/work" -qs 100 \
    -c "$run/workers.config" -params-file "$run/params.json" \
    --outdir "$run/results"

python3 "$helper" verify --expected "$run/expected.json" --results "$run/results" --out "$run/acceptance.json"
tar -czf "$run.review.tar.gz" --exclude='*.bam' --exclude='*.bai' \
    --exclude='*.fasta' --exclude='*.fasta.gz' --exclude='*.fa' --exclude='blobdir' \
    -C "$run" results/assembly/contamination results/pipeline \
    params.json expected.json workers.config acceptance.json \
    "nextflow-${SLURM_JOB_ID}.log" "slurm-${SLURM_JOB_ID}.out"
printf '\nReview archive: %s\n' "$run.review.tar.gz"
