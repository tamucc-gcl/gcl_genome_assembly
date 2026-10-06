#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=ctlk_yahs_r1
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --output=logs/ctlk-yahs-r1-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
tool_env=${1:?Supply environment with yahs and samtools}
assessment=${2:-comparisons/chimera-foundation-01}
# Keep enzyme correction identical across treatments. Supply 'none' if the
# retained baseline round-1 command used no enzyme correction.
enzyme=${3:-GATC,CTNAG,TTAA,GANTC}
module load miniconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$tool_env"
hap=$(( SLURM_ARRAY_TASK_ID % 2 + 1 ))
variant=$(( SLURM_ARRAY_TASK_ID / 2 ))
labels=(baseline cap10Mb cap10Mb_min100kb)
label=${labels[$variant]}
sample="Sde-CTlk_104_hap${hap}"
root=$(pwd)
out="comparisons/ctlk-yahs-round1-${SLURM_ARRAY_JOB_ID}/${label}_hap${hap}"
fasta="$root/$assessment/assembly/contig/decontam/${sample}.decontaminated.fasta"
bam="$root/$assessment/bam/hic/contig/filtered/${sample}.filtered.sorted.bam"
[[ -s "$fasta" && -s "$bam" ]] || { echo "Missing exact round-1 input: $fasta or $bam"; exit 1; }
mkdir -p "$(dirname "$out")"
mkdir "$out"
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
# Verify BAM reference names AND lengths against this unscaffolded FASTA.
samtools view -H "$bam" | awk -F '\t' '$1=="@SQ" {sn=""; ln=""; for(i=2;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4);if($i~/^LN:/)ln=substr($i,4)} print sn "\t" ln}' | LC_ALL=C sort > bam.references.tsv
cut -f1,2 input.fa.fai | LC_ALL=C sort > fasta.references.tsv
diff -u fasta.references.tsv bam.references.tsv > reference_check.txt || { echo 'BAM/FASTA reference mismatch'; exit 1; }
resolutions='20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000'
if [[ "$variant" == 0 ]]; then resolutions+=',20000000,50000000,100000000,200000000'; fi
min_len=20000
if [[ "$variant" == 2 ]]; then min_len=100000; fi
args=(yahs input.fa "$bam" -o "$sample" -l "$min_len" -q 10 -r "$resolutions")
if [[ "$enzyme" != none ]]; then args+=(-e "$enzyme"); fi
# Current round-1 baseline already enables both EC stages.
printf '%q ' "${args[@]}" > command.txt
printf '\n' >> command.txt
printf 'fasta\t%s\nbam\t%s\n' "$fasta" "$bam" > input_paths.tsv
stat -Lc '%n\t%s\t%y' "$fasta" "$bam" > input_metadata.tsv
yahs --version > versions.txt 2>&1
"${args[@]}" 2>&1 | tee yahs.log
