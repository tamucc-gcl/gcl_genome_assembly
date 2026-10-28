#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=core_integration
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/core-integration-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
samples="${1:?Supply the sample sheet}"
hic="${2:?Supply the Hi-C readsets table}"
export CORE_INTEGRATION_BASELINE="$(realpath "${3:?Supply accepted results directory}")"
cuts="$(realpath "${4:?Supply the accepted manual cut sheet}")"
out="${5:?Supply the integration output directory}"
[[ "$out" == /* ]] || out="$root/$out"
mkdir -p "$out" logs
cat > "$out/integration.config" <<EOF
trace {
  fields = 'task_id,hash,native_id,name,status,exit,submit,duration,realtime,%cpu,peak_rss,peak_vmem,rchar,wchar,workdir,container'
}
EOF
bash "$root/gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch" \
    "$samples" "$hic" "$out" \
    --chimera_break "$cuts" \
    --run_pangenome true --pangenome_validate_only false --pangenome_regional_chunk_bp 25000 \
    --run_post_assembly true --run_pairwise_alignments true --pairwise_alignment_mode within_sample \
    --pangenome_2d_viz false -c "$out/integration.config" "${@:6}"
