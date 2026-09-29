#!/usr/bin/env bash
# Collect existing outputs only. Run from the project root.
set -euo pipefail
old_results="${1:-tst/genome_assembly_store}"
new_results="${2:-genome_assembly}"
sample_sheet="${3:-data/assembly_samplesheet.csv}"
hic_readsets="${4:-data/hic_readsets.csv}"
for input in "$sample_sheet" "$hic_readsets"; do
    [[ "$input" != /* && "$input" != *".."* && -f "$input" ]] || {
        echo "Supply existing project-relative input tables (arguments 3 and 4): $input" >&2
        exit 1
    }
done
for root in "$old_results" "$new_results"; do
    [[ "$root" != /* && "$root" != *".."* && -d "$root" ]] || {
        echo "Supply existing project-relative result directories: $root" >&2
        exit 1
    }
done
stamp=$(date +%Y%m%d-%H%M%S)
mkdir -p comparisons
prefix="comparisons/chimera-comparison-${stamp}"
mkdir "$prefix"
reports=$(mktemp)
alignments=$(mktemp)
trap 'rm -f "$reports" "$alignments"' EXIT
printf 'role\tresults_directory\nold\t%s\nnew\t%s\n' "$old_results" "$new_results" > "$prefix/sources.tsv"
printf 'role\tprovided_input_file\nsample_sheet\t%s\nhic_readsets\t%s\n' "$sample_sheet" "$hic_readsets" > "$prefix/input_files.tsv"
git -C gcl_genome_assembly rev-parse HEAD > "$prefix/current_checkout_commit.txt"
git -C gcl_genome_assembly status --short > "$prefix/current_checkout_status.txt"
git -C gcl_genome_assembly diff > "$prefix/current_checkout.diff"
# Current checkout does not prove which software generated stored results.
for root in "$old_results" "$new_results"; do
    find -L "$root" -type f -printf '%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS\n' >> "$prefix/inventory.tsv"
    find -L "$root" -type f \( \
        -path '*/harmonization/*' -o -path '*/chimeras/*' -o \
        -name '*chimera*' -o -name '*harmoniz*' -o \
        -name '*name_map*' -o -name '*.agp' -o -name '*.fai' -o \
        -name 'assembly_run_summary.md' -o -name 'assembly_eligibility.*' -o \
        -name '*trace*.txt' -o -name '*trace*.tsv' \
    \) \( -name '*.tsv' -o -name '*.json' -o -name '*.md' -o \
           -name '*.txt' -o -name '*.tidk.log' -o -name '*.png' -o -name '*.agp' -o -name '*.fai' \
    \) -print0 >> "$reports"
    find -L "$root" -type f -name '*.ref.paf.gz' -print0 >> "$alignments"
done
# Collect all assembly logs; the latest job may have reused earlier results.
if [[ -d logs ]]; then
    find logs -maxdepth 1 -type f \( -name '*assembly*.log*' -o \
        -name 'assembly-*.out' \) -print0 >> "$reports"
fi
for f in "$sample_sheet" "$hic_readsets" \
         gcl_genome_assembly/nextflow.config gcl_genome_assembly/run_assembly_checkpoint.sbatch; do
    [[ ! -f "$f" ]] || printf '%s\0' "$f" >> "$reports"
done
printf '%s\0' "$prefix" >> "$reports"
sort -zu "$reports" -o "$reports"
tar --dereference --null -czf "${prefix}-reports.tar.gz" -T "$reports"
if [[ -s "$alignments" ]]; then
    sort -zu "$alignments" -o "$alignments"
    tar --dereference --null -czf "${prefix}-alignments.tar.gz" -T "$alignments"
fi
ls -lh "${prefix}"-*.tar.gz
echo "Attach the reports archive first; retain the alignments archive for follow-up."

