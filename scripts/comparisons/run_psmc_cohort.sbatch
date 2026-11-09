#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --job-name=psmc_cohort
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/psmc-cohort-%j.out
set -euo pipefail
root="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$root"
repo="$root/gcl_genome_assembly"
samples="${1:?Supply sample sheet}"
hic="${2:?Supply Hi-C readsets table}"
export CORE_INTEGRATION_BASELINE=$(realpath "${3:?Supply accepted main-run results}")
cuts=$(realpath "${4:?Supply reviewed manual cut sheet}")
out=$(realpath -m "${5:?Supply a new cohort output directory}")
length_run=$(realpath "${6:?Supply completed length-test directory}")
hifi_run=$(realpath "${7:?Supply completed HiFi comparison directory}")
if [[ -e "$out" ]]; then printf 'Use a new cohort output directory: %s\n' "$out" >&2; exit 2; fi
mkdir -p "$out" logs
finish() {
    rc=$?; trap - EXIT
    python3 "$repo/scripts/comparisons/psmc_cohort.py" collect --root "$out" || { [[ "$rc" -ne 0 ]] || rc=1; }
    printf '%s\n' "$rc" > "$out/cohort_exit_status.txt"
    exit "$rc"
}
trap finish EXIT
python3 "$repo/scripts/comparisons/psmc_cohort.py" prepare --baseline "$CORE_INTEGRATION_BASELINE" --out "$out"
cat > "$out/integration.config" <<'EOF'
trace.fields = 'task_id,hash,native_id,name,status,exit,submit,duration,realtime,%cpu,peak_rss,peak_vmem,rchar,wchar,workdir,container'
EOF
bash "$repo/scripts/comparisons/run_chimera_pipeline.sbatch" "$samples" "$hic" "$out" \
    -params-file "$out/baseline-params.json" --chimera_break "$cuts" \
    --run_pangenome true --pangenome_validate_only false --pangenome_regional_chunk_bp 25000 \
    --run_post_assembly true --run_pairwise_alignments true --pairwise_alignment_mode within_sample \
    --pangenome_2d_viz false -c "$out/integration.config" "${@:8}"
python3 "$repo/scripts/comparisons/psmc_cohort.py" audit --root "$out"
python3 "$repo/scripts/comparisons/psmc_cohort.py" context --root "$out" --length-run "$length_run" --hifi-run "$hifi_run"
pids=()
while IFS= read -r sample; do
    sbatch --wait --parsable --output="$out/psmc_checks/support-$sample-%j.out" \
        "$repo/scripts/comparisons/run_psmc_input_support.sbatch" "$out" "$sample" "$repo" &
    pids+=("$!")
done < <(python3 -c 'import json,sys; print("\n".join(j["sample"] for j in json.load(open(sys.argv[1]))["jobs"]))' "$out/results/post_assembly/psmc/psmc_plan.json")
for left in control_20kb hifi; do
    label=20kb_vs_cohort
    [[ "$left" != hifi ]] || label=hifi_vs_cohort
    sbatch --wait --parsable --output="$out/psmc_checks/$label-%j.out" \
        "$repo/scripts/comparisons/run_psmc_reliability.sbatch" "$out/psmc_checks" "$out/psmc_checks/$label" \
        --left-run "$left" --right-run cohort --hifi-run "$hifi_run" &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
[[ "$failed" == 0 ]]
