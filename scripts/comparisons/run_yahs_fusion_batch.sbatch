#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=ctlk_yahs_batch
#SBATCH --array=0-3%3
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/ctlk-yahs-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
tool_env=${1:?Supply environment with yahs and samtools}
assessment=${2:-comparisons/chimera-foundation-01}
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
hap=$(( SLURM_ARRAY_TASK_ID % 2 + 1 ))
variant=$(( SLURM_ARRAY_TASK_ID / 2 ))
labels=(ec_on ec_on_cap10Mb)
label=${labels[$variant]}
sample="Sde-CTlk_104_hap${hap}"
out="comparisons/ctlk-yahs-batch-${SLURM_ARRAY_JOB_ID}/${label}_hap${hap}"
mkdir -p "$(dirname "$out")"
mkdir "$out"
root=$(pwd)
fasta="$root/$assessment/assembly/scaffold/misassembly_correction/${sample}_corrected.fasta"
bam="$root/$assessment/bam/hic/scaffold/filtered/${sample}.filtered.sorted.bam"
[[ -s "$fasta" && -s "$bam" ]] || { echo "Missing exact retained input: $fasta or $bam"; exit 1; }
finish() {
    rc=$?
    trap - EXIT TERM INT HUP
    printf '%s\n' "$rc" > "$root/$out/exit_status.txt"
    tar --exclude='*.fa' --exclude='*.fasta' --exclude='*.bam' --exclude='*.bin' \
        -czf "$root/$out.review.tar.gz" -C "$root" "$out" || rc=1
    exit "$rc"
}
trap finish EXIT
trap 'exit 143' TERM HUP
trap 'exit 130' INT
cd "$out"
ln -s "$fasta" input.fa
samtools faidx input.fa
resolutions='50000,100000,200000,500000,1000000,2000000,5000000,10000000'
if [[ "$variant" == 0 ]]; then resolutions+=',20000000,50000000,100000000'; fi
# Both variants ENABLE scaffold EC by omitting --no-scaffold-ec.
# The second is an experimental coarse-resolution cap, not a production default.
args=(yahs input.fa "$bam" -o "$sample" -l 100000 -q 20 -r "$resolutions" -e 'GATC,CTNAG,TTAA,GANTC')
printf '%q ' "${args[@]}" > command.txt
printf '\n' >> command.txt
yahs --version > versions.txt 2>&1
"${args[@]}" 2>&1 | tee yahs.log
