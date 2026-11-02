#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=scaf_replay
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/scaffolding-replay-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
repo="$root/gcl_genome_assembly"
assembly="${1:?assembly ID}"
old_fasta="$(realpath "${2:?old decontaminated FASTA}")"
new_fasta="$(realpath "${3:?new decontaminated FASTA}")"
old_bin="$(realpath "${4:?old YaHS BIN}")"
new_bin="$(realpath "${5:?new YaHS BIN}")"
old_agp="$(realpath "${6:?accepted round1 AGP}")"
out="${7:?output directory}"
mkdir -p "$out"
out="$(realpath "$out")"
finish() {
    rc=$?
    trap - EXIT
    printf '%s\n' "$rc" > "$out/exit_status.txt"
    files=(exit_status.txt)
    for f in input_identity.json replay_check.json nextflow.log trace.tsv results; do
        [[ -e "$out/$f" ]] && files+=("$f")
    done
    tar --dereference --exclude='*.fa' --exclude='*.fasta' --exclude='*.bin' \
        -czf "$out.review.tar.gz" -C "$out" "${files[@]}"
    exit "$rc"
}
trap finish EXIT
python3 - "$old_fasta" "$new_fasta" "$old_bin" "$new_bin" "$out" <<'PY'
import hashlib,json,sys
from pathlib import Path
def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
paths=sys.argv[1:5]
hashes=[digest(p) for p in paths]
result=dict(paths=paths,sha256=hashes,fasta_match=hashes[0]==hashes[1],links_match=hashes[2]==hashes[3])
(Path(sys.argv[5])/'input_identity.json').write_text(json.dumps(result,indent=2)+'\n')
if not result['fasta_match'] or not result['links_match']:
    raise SystemExit('Old/new inputs differ; inspect input_identity.json before replay')
PY
if ! command -v nextflow >/dev/null; then module load nextflow; fi
cd "$out"
nextflow -log "$out/nextflow.log" run "$repo/scripts/comparisons/scaffolding_replay.nf" \
    -c "$repo/nextflow.config" -profile slurm -w "$out/work" \
    --outdir "$out/results" --publish_dir_mode symlink \
    --replay_assembly "$assembly" --replay_fasta "$new_fasta" \
    --replay_links "$new_bin" --replay_expected_agp "$old_agp" \
    -with-trace "$out/trace.tsv" "${@:8}"
