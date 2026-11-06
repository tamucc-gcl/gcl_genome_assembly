#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_mask
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
set -euo pipefail
if type module >/dev/null 2>&1; then module load miniconda3; fi
manifest=$(realpath "${1:?Supply manifest.json}")
repo=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["repository"])' "$manifest")
mode="${2:-task}"
if [[ "$mode" == report ]]; then
    out=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["out"])' "$manifest")
    finish() {
        rc=$?
        trap - EXIT
        printf '%s\n' "$rc" > "$out/exit_status.txt"
        # Keep complete compact evidence even if a task or the renderer failed.
        # Raw PAFs/FASTA links stay outside the review archive.
        tar --exclude='*/raw' --exclude='*.paf.gz' --exclude='*.fa' --exclude='*.fasta' \
            -czf "$out.review.tar.gz.tmp" -C "$out" . && mv "$out.review.tar.gz.tmp" "$out.review.tar.gz" || {
                [[ "$rc" -ne 0 ]] || rc=1
            }
        exit "$rc"
    }
    trap finish EXIT
    python3 "$repo/scripts/comparisons/psmc_mask_investigation.py" report --manifest "$manifest"
else
    index="${SLURM_ARRAY_TASK_ID:-${3:?Supply a task index when running interactively}}"
    prefix=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["tasks"][int(sys.argv[2])]["env"])' "$manifest" "$index")
    export PATH="$prefix/bin:$PATH"
    python3 "$repo/scripts/comparisons/psmc_mask_investigation.py" task \
        --manifest "$manifest" --index "$index" --cpus "${SLURM_CPUS_PER_TASK:-8}"
fi
