#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_sequence
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/chimera-sequence-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
run= tool_env= region= motif=CCCTAA
labels=(); fastas=()
while (( $# )); do
    case "$1" in
        --mapping-run) run="$2"; shift 2 ;;
        --tool-env) tool_env="$2"; shift 2 ;;
        --region) region="$2"; shift 2 ;;
        --motif) motif="$2"; shift 2 ;;
        --comparison) labels+=("$2"); fastas+=("$3"); shift 3 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -s "$run/hifi.bam" && -s "$run/assessment.fa" && "$region" =~ ^[^:]+:[0-9]+-[0-9]+$ ]] || {
    echo 'Required: --mapping-run DIR --region scaffold:start-end [--tool-env PREFIX] [--comparison LABEL FASTA] [--motif CCCTAA]' >&2; exit 2;
}
if [[ -n "$tool_env" ]]; then
    module load miniconda3
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$tool_env"
fi
for tool in samtools minimap2 tidk sha256sum; do command -v "$tool" >/dev/null; done
for i in "${!labels[@]}"; do
    [[ "${labels[$i]}" =~ ^[A-Za-z0-9_-]+$ && -s "${fastas[$i]}" ]] || exit 2
    for ((j=0;j<i;j++)); do [[ "${labels[$i]}" != "${labels[$j]}" ]] || { echo 'Duplicate label' >&2; exit 2; }; done
done
expected=$(awk -F '\t' '$1=="sha256"{print $2}' "$run/reports/provenance.tsv")
actual=$(sha256sum "$run/assessment.fa" | awk '{print $1}')
[[ -n "$expected" && "$actual" == "$expected" ]] || { echo 'Assessment checksum mismatch' >&2; exit 1; }
samtools quickcheck "$run/hifi.bam"
out="comparisons/chimera-sequence-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out"
printf 'mapping_run\t%s\nregion_1based_inclusive\t%s\nmotif\t%s\nassessment_sha256\t%s\n' "$(realpath "$run")" "$region" "$motif" "$actual" > "$out/provenance.tsv"
samtools --version > "$out/samtools.version.txt"
minimap2 --version > "$out/minimap2.version.txt"
tidk --version > "$out/tidk.version.txt"
cp "$run/reports/assessment.fa.fai" "$out/assessment.fa.fai"
chrom=${region%%:*}; interval=${region#*:}; start=${interval%-*}; end=${interval#*-}
awk -v c="$chrom" -v s="$start" -v e="$end" 'BEGIN{ok=0} $1==c && s>=1 && e>=s && e<=$2{ok=1} END{exit !ok}' "$out/assessment.fa.fai"
# Whole scaffold profile preserves coordinates and includes its terminal controls.
samtools faidx "$run/assessment.fa" "$chrom" > "$out/scaffold.fa"
samtools faidx "$run/assessment.fa" "$region" > "$out/region.fa"
tidk search --string "$motif" --window 1000 --output assessment --dir "$out/tidk" "$out/scaffold.fa" > "$out/tidk-assessment.log" 2>&1
# Region-only motif profile at finer resolution. Coordinates are relative to region start.
tidk search --string "$motif" --window 100 --output region --dir "$out/tidk" "$out/region.fa" > "$out/tidk-region.log" 2>&1
samtools view "$run/hifi.bam" "$region" | cut -f1 | LC_ALL=C sort -u > "$out/read_names.txt"
# Fetch complete primary sequences for every read with any alignment in the region.
# Reading the complete BAM preserves reads whose primary placement is elsewhere.
samtools view -u -N "$out/read_names.txt" -F 2304 "$run/hifi.bam" |
    samtools fasta - > "$out/reads.fa"
tidk search --string "$motif" --window 100 --output reads --dir "$out/tidk" "$out/reads.fa" > "$out/tidk-reads.log" 2>&1
samtools view -N "$out/read_names.txt" "$run/hifi.bam" |
    awk 'BEGIN{FS=OFS="\t";print "qname","flag","rname","pos_1based","mapq","cigar","rnext","pnext","tlen","tags"} {for(i=1;i<=9;i++)printf "%s\t",$i;for(i=12;i<=NF;i++)printf "%s%s",$i,(i<NF?";":"");print ""}' | gzip > "$out/alignments.tsv.gz"
for i in "${!labels[@]}"; do
    label=${labels[$i]}; fasta=${fastas[$i]}
    printf 'comparison_%s\t%s\n' "$label" "$(realpath "$fasta")" >> "$out/provenance.tsv"
    sha256sum "$fasta" > "$out/${label}.sha256.txt"
    ln -s "$(realpath "$fasta")" "$out/${label}.fa"
    samtools faidx "$out/${label}.fa"
    total=$(awk '{n+=$2} END{printf "%.0f",n}' "$out/${label}.fa.fai")
    (( total < 8000000000 )) || { echo 'Comparison exceeds single-index budget' >&2; exit 2; }
    # Full comparison assemblies retain alternative sequence placements.
    minimap2 -x asm5 -I 8G -c --cs=long --secondary=yes -N 20 -t "${SLURM_CPUS_PER_TASK:-8}" "$fasta" "$out/region.fa" |
        gzip > "$out/${label}.region.paf.gz"
    tidk search --string "$motif" --window 1000 --output "$label" --dir "$out/tidk" "$fasta" > "$out/tidk-${label}.log" 2>&1
done
gzip "$out/reads.fa" "$out/region.fa"
# Retain the large scaffold locally; archive sequence of the requested region and reads.
printf 'SUCCESS\n' > "$out/status.txt"
tar --exclude='*.fa' -czf "${out}.tar.gz" "$out"
printf '\nDownload %s.tar.gz and logs/chimera-sequence-%s.out\n' "$out" "${SLURM_JOB_ID:-manual}"
