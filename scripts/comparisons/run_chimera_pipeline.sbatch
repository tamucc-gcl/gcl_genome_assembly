#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_pipeline
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/chimera-pipeline-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$root/gcl_genome_assembly"
sample_sheet="${1:-data/assembly_samplesheet.csv}"
hic_readsets="${2:-data/hic_readsets.csv}"
out="${3:-$root/comparisons/chimera-pipeline-${SLURM_JOB_ID}}"
[[ "$out" == /* ]] || out="$root/$out"
mkdir -p "$out" "$root/logs"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    if [[ -f "$root/logs/nextflow_assembly_${SLURM_JOB_ID}.log" ]]; then
        cp "$root/logs/nextflow_assembly_${SLURM_JOB_ID}.log" "$out/nextflow.log"
    fi
    # Package reports and evidence, not the full assembly output tree.
    files=(exit_status.txt)
    for name in nextflow.log trace.tsv results/assembly_run_summary.md results/assembly_report.md results/assembly_eligibility.json results/assembly_eligibility.tsv results/assembly_eligibility.md results/references results/assembly/chimeras results/assembly/scaffold/teloclip results/qc/assembly results/qc/hifi results/qc/hic results/qc/shortread results/snail_plots results/reports; do
        [[ -e "$out/$name" ]] && files+=("$name")
    done
    # Final names are needed to review the chromosome inventory after manual cuts.
    # Include these small tables without packaging the final FASTAs or mappings.
    for name in "$out"/results/assembly/final/*.name_map.tsv; do
        [[ -f "$name" ]] && files+=("${name#"$out/"}")
    done
    # Review contact-map images/statistics without packaging large cooler matrices.
    for name in "$out"/results/contact_maps/*_final_*bp_contact_map.png "$out"/results/contact_maps/*_final_contact_stats.txt; do
        [[ -f "$name" ]] && files+=("${name#"$out/"}")
    done
    # Compact graph reports, provenance, summary tables and figures. Graphs, indexes,
    # alignments, VCFs and large per-node tables remain in the selected results directory.
    if [[ -d "$out/results/pangenome" ]]; then
        while IFS= read -r -d '' name; do files+=("${name#"$out/"}"); done < <(
            find -L "$out/results/pangenome" -type f \( -size -20M -o \
                -name 'coverage_*.tsv.gz' -o -name 'surviving_intervals.tsv.gz' \) \( \
                -name '*.md' -o -name '*.json' -o -name '*.tsv' -o -name 'coverage_*.tsv.gz' -o -name 'surviving_intervals.tsv.gz' -o -name '*.csv' -o \
                -name '*.txt' -o -name '*filter.log' -o -name '*.xml' -o -name '*.png' -o -name '*.svg' -o -name '*.pdf' -o \
                -name '*.nwk' -o -name '*.stats.tgz' \) -print0)
    fi
    tar --dereference --exclude='*.fa' --exclude='*.fasta' --exclude='*.fa.gz' --exclude='*.fasta.gz' \
        --exclude='*.bam' --exclude='*.bai' --exclude='*.csi' --exclude='*.sam' \
        -czf "$out.review.tar.gz" -C "$out" "${files[@]}"
    exit "$rc"
}
trap finish EXIT

# This is the normal main.nf launch, with its full sample/read channels, process
# environments and SLURM resource requests. Reuse the production launch/cache
# and work paths; publish this test's outputs in a separate results directory.
# Symbolically link cached files AND directories rather than copying them.
# Nextflow's `link` mode means hard links and still uses its transfer executor.
# The review archive above dereferences links so it remains self-contained.
bash "$repo/run_assembly_checkpoint.sbatch" "$sample_sheet" "$hic_readsets" assembly \
    --outdir "$out/results" --publish_dir_mode symlink \
    --qc_mode all_stages \
    --chimera_detect true --chimera_break false \
    --chimera_plots true --chimera_hifi_context true \
    -with-trace "$out/trace.tsv" -with-report "$out/execution.html" \
    "${@:4}"

test -s "$out/results/assembly/chimeras/README.md"
test -s "$out/results/assembly/chimeras/assembly_registry.tsv"
test -s "$out/results/assembly/chimeras/review.tsv"
test -s "$out/results/assembly_run_summary.md"
