#!/usr/bin/env bash
# Read-only collection of small evidence files. Run from the cluster project root.
set -euo pipefail
job=${1:?Usage: bash collect_chimera_review.sh JOB_ID [results] [samplesheet] [hic_readsets]}
results=${2:-genome_assembly}
samples=${3:-data/assembly_samplesheet.csv}
readsets=${4:-data/hic_readsets.csv}
[[ "$job" =~ ^[0-9]+$ ]] || { echo 'JOB_ID must be numeric' >&2; exit 1; }
for entry in "$results" "$samples" "$readsets"; do
    [[ "$entry" != /* && "$entry" != *'..'* && -e "$entry" ]] || {
        echo "Expected an existing project-relative path: $entry" >&2; exit 1;
    }
done
prefix="comparisons/chimera-review-$(date +%Y%m%d-%H%M%S)-${job}"
mkdir -p "$prefix"
list="$prefix/archive-files.nul"
find -L "$results" -type f -printf '%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS\n' > "$prefix/inventory.tsv"
find -L "$results" -type f \( \
    -path '*/harmonization/*' -o -path '*/chimeras/*' -o \
    -name '*harmoniz*' -o -name '*chimera*' -o -name '*name_map*' -o \
    -name '*.agp' -o -name '*.fai' -o -name 'assembly_eligibility.*' -o \
    -name 'assembly_run_summary.md' -o -name '*trace*.txt' -o -name '*trace*.tsv' \
    \) \( -name '*.tsv' -o -name '*.csv' -o -name '*.json' -o -name '*.md' -o \
    -name '*.txt' -o -name '*.log' -o -name '*.png' -o -name '*.agp' -o -name '*.fai' \) \
    -print0 > "$list"
for entry in "$samples" "$readsets" scripts/run_assembly.sbatch \
    "logs/assembly-${job}.out" "logs/nextflow_assembly_${job}.log" \
    gcl_genome_assembly/nextflow.config; do
    if [[ -f "$entry" ]]; then
        printf '%s\0' "$entry" >> "$list"
    else
        printf 'Not found: %s\n' "$entry" >> "$prefix/missing-files.txt"
    fi
done
printf 'job_id\t%s\nresults\t%s\nsamplesheet\t%s\nhic_readsets\t%s\n' \
    "$job" "$results" "$samples" "$readsets" > "$prefix/requested-run.tsv"
git -C gcl_genome_assembly rev-parse HEAD > "$prefix/current-checkout.txt"
git -C gcl_genome_assembly status --short >> "$prefix/current-checkout.txt"
git -C gcl_genome_assembly diff > "$prefix/current-checkout.diff"
# Current checkout/published files are not themselves proof of run provenance.
find "$prefix" -maxdepth 1 -type f ! -name 'archive-files.nul' -print0 >> "$list"
sort -zu "$list" -o "$list"
tar --dereference --null -czf "${prefix}.tar.gz" -T "$list"
ls -lh "${prefix}.tar.gz"
echo 'Return this archive. Keep work/cache, FASTAs, BAMs and alignment files on the cluster.'
