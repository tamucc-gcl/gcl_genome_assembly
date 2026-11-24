#!/bin/bash
# Run interactively from the project root AFTER syncing the cache-guard fix.
# Arguments: controller job ID, pending worker job ID, completed PLAN's plan.json.
set -euo pipefail
if [[ $# -ne 3 || ! "$1" =~ ^[0-9]+$ || ! "$2" =~ ^[0-9]+$ ]]; then
    echo 'Usage: bash recover_contamination_staged_cohort.sh CONTROLLER WORKER PLAN_JSON' >&2
    exit 1
fi
controller="$1"
pending="$2"
plan="$3"
project="$PWD"
module load miniconda3
for job in "$controller" "$pending"; do
    owner=$(squeue -h -j "$job" -o '%u')
    if [[ "$owner" != "$USER" ]]; then
        echo "Job $job is absent or does not belong to $USER; stopping." >&2
        exit 1
    fi
done
controller_name=$(squeue -h -j "$controller" -o '%j')
worker_name=$(squeue -h -j "$pending" -o '%j')
if [[ "$controller_name" != contam_cohort || "$worker_name" != *CONTAMINATION_MEGABLAST_NODE* ]]; then
    echo 'Job names do not match the cohort controller and node worker; stopping.' >&2
    exit 1
fi
state=$(squeue -h -j "$pending" -o '%T')
if [[ "$state" != PENDING ]]; then
    echo 'Worker is no longer pending; recheck its status before recovery.' >&2
    exit 1
fi
python3 "$project/gcl_genome_assembly/scripts/comparisons/check_contamination_cohort_recovery.py" \
    --project "$project" --plan "$plan" --controller "$controller"

scancel --signal=TERM "$controller"
scancel "$pending" || true
for attempt in {1..30}; do
    if [[ -z "$(squeue -h -j "$controller" -o '%i' 2>/dev/null || true)" ]]; then
        sbatch --parsable "$project/gcl_genome_assembly/scripts/comparisons/run_contamination_staged_cohort.sh"
        exit 0
    fi
    sleep 2
done
echo 'Controller still stopping. Wait for it to exit, then submit run_contamination_staged_cohort.sh.' >&2
exit 1
