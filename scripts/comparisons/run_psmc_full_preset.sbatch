#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_preset
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
module load nextflow/23.10.1
module load miniconda3
comparison=$(realpath "${1:?Supply the prepared full-comparison directory}")
preset="${2:?Supply the prepared run label}"
[[ "$preset" =~ ^[A-Za-z0-9_-]+$ ]]
repo=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["repository"])' "$comparison/comparison.json")
out="$comparison/$preset"
test -s "$out/params.json"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    files=(params.json execution.config exit_status.txt)
    for name in results trace.tsv nextflow.log; do [[ ! -e "$out/$name" ]] || files+=("$name"); done
    tar --dereference --exclude='*.fa' --exclude='*.fasta' --exclude='*.fa.gz' --exclude='*.fasta.gz' \
        --exclude='*.bam*' --exclude='*.bai' --exclude='*.vcf*' --exclude='*.psmcfa' \
        -czf "$out.review.tar.gz" -C "$out" "${files[@]}" || { [[ "$rc" -ne 0 ]] || rc=1; }
    exit "$rc"
}
trap finish EXIT
cat > "$out/execution.config" <<'EOF'
trace.fields = 'task_id,hash,native_id,name,status,exit,submit,duration,realtime,%cpu,peak_rss,workdir'
EOF
# Independent launch/cache directories permit both comparisons and the existing
# HiFi workflow to run concurrently without sharing a live Nextflow session.
mkdir -p "$out/launch"
cd "$out/launch"
nextflow -log "$out/nextflow.log" run "$repo/psmc_from_final.nf" -profile slurm \
    -params-file "$out/params.json" -c "$out/execution.config" -w "$out/work" \
    -with-trace "$out/trace.tsv" -with-report "$out/execution.html"
