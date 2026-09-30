#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=chimera_hifi
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/chimera-hifi-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONDONTWRITEBYTECODE=1
fasta= reads= expected= region= tool_env=
while (( $# )); do
    case "$1" in
        --fasta) fasta="$2"; shift 2 ;;
        --reads) reads="$2"; shift 2 ;;
        --sha256) expected="$2"; shift 2 ;;
        --region) region="$2"; shift 2 ;;
        --tool-env) tool_env="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -s "$fasta" && -s "$reads" && "$expected" =~ ^[[:xdigit:]]{64}$ && "$region" =~ ^[^:]+:[0-9]+-[0-9]+$ ]] || {
    echo 'Required: --fasta complete-assessment.fa --reads reads.fastq.gz --sha256 HASH --region scaffold:start-end [--tool-env PREFIX]' >&2; exit 2;
}
if [[ -n "$tool_env" ]]; then
    module load miniconda3
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$tool_env"
fi
for tool in minimap2 samtools sha256sum awk gzip; do command -v "$tool" >/dev/null; done
actual=$(sha256sum "$fasta" | awk '{print $1}')
[[ "$actual" == "$expected" ]] || { echo 'Assessment checksum mismatch' >&2; exit 1; }
threads=${SLURM_CPUS_PER_TASK:-16}
(( threads >= 4 )) || { echo 'Request at least four CPUs' >&2; exit 2; }
out="comparisons/chimera-hifi-$(date +%Y%m%d-%H%M%S)-${SLURM_JOB_ID:-manual}"
mkdir -p "$out/reports" "$out/tmp"
exec > >(tee "$out/reports/diagnostic.log") 2>&1
printf 'fasta\t%s\nsha256\t%s\nreads\t%s\nregion_1based_inclusive\t%s\n' "$(realpath "$fasta")" "$actual" "$(realpath "$reads")" "$region" > "$out/reports/provenance.tsv"
minimap2 --version > "$out/reports/minimap2.version.txt"
samtools --version > "$out/reports/samtools.version.txt"
# Local symlink avoids creating an index beside a published/hard-linked FASTA.
ln -s "$(realpath "$fasta")" "$out/assessment.fa"
samtools faidx "$out/assessment.fa"
cp "$out/assessment.fa.fai" "$out/reports/assessment.fa.fai"
chrom=${region%%:*}; interval=${region#*:}; start=${interval%-*}; end=${interval#*-}
awk -v c="$chrom" -v s="$start" -v e="$end" 'BEGIN {ok=0} $1==c && s>=1 && e>=s && e<=$2 {ok=1} END {exit !ok}' "$out/assessment.fa.fai"
# Reject multi-part indices: mapping qualities are not reliable across parts.
total=$(awk '{n+=$2} END {printf "%.0f",n}' "$out/assessment.fa.fai")
(( total < 8000000000 )) || { echo 'Assembly exceeds this diagnostic single-index budget; adjust resources/index explicitly.' >&2; exit 2; }
printf '%s\n' 'minimap2 -ax map-hifi -I 8G --secondary=yes -N 20; samtools sort' > "$out/reports/mapping_recipe.txt"
minimap2 -ax map-hifi -I 8G --secondary=yes -N 20 -t "$((threads-3))" "$out/assessment.fa" "$reads" 2> "$out/reports/minimap2.log" |
    samtools sort -@ 2 -m 2G -T "$out/tmp/sort" -o "$out/hifi.bam" -
samtools quickcheck -v "$out/hifi.bam"
samtools index -@ 2 "$out/hifi.bam"
samtools flagstat "$out/hifi.bam" > "$out/reports/flagstat.txt"
samtools idxstats "$out/hifi.bam" > "$out/reports/idxstats.tsv"
samtools view "$out/hifi.bam" "$region" | cut -f1 | LC_ALL=C sort -u > "$out/reports/regional_read_names.txt"
# All reported placements for regional reads, including elsewhere in the assembly.
# This is a TSV, not SAM: SEQ and QUAL are deliberately omitted to keep transfer small.
samtools view -N "$out/reports/regional_read_names.txt" "$out/hifi.bam" |
    awk 'BEGIN {FS=OFS="\t"; print "qname","flag","rname","pos_1based","mapq","cigar","rnext","pnext","tlen","tags"} {for(i=1;i<=9;i++) printf "%s\t",$i; for(i=12;i<=NF;i++) printf "%s%s",$i,(i<NF?";":""); print ""}' |
    gzip > "$out/reports/regional_read_alignments.tsv.gz"
# Primary alignments only: supplementary records remain available above for inspection.
# depth's -q is base quality; -Q is mapping quality. Report two MAPQ filters.
for mapq in 0 20; do
    samtools depth -aa -r "$region" -q 0 -Q "$mapq" -G 0xF04 "$out/hifi.bam" |
        awk 'BEGIN {OFS="\t"; print "chrom","start_0based","end_exclusive","bases","mean_depth","zero_depth_bases"} function emit(){if(n) print c,first,last,n,sum/n,z} {b=int(($2-1)/1000); if(n && (b!=prev || $1!=c)){emit();n=sum=z=0} if(!n){c=$1;first=$2-1} last=$2;prev=b;n++;sum+=$3;if($3==0)z++} END{emit()}' > "$out/reports/depth_mapq${mapq}_1kb.tsv"
done
printf 'SUCCESS\n' > "$out/reports/status.txt"
tar -czf "${out}-reports.tar.gz" "$out/reports"
printf '\nDownload %s-reports.tar.gz and logs/chimera-hifi-%s.out\nFull BAM retained at %s/hifi.bam\n' "$out" "${SLURM_JOB_ID:-manual}" "$out"
