#!/usr/bin/env bash
set -euo pipefail
runtime=/tmp/chimera-nxf-runtime
mkdir -p "$runtime" "$runtime/assets"
touch "$runtime/assets/NO_PAF" "$runtime/assets/NO_PAIRS"
mkdir -p "$runtime/context"
printf '{"sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"}\n' > "$runtime/context/provenance.json"
printf 'metric\tvalue\nalignment_transition_assessment\tassessed\n' > "$runtime/audit.tsv"
downloads="$(dirname "$0")/../comparisons/nxf-runtime"
mkdir -p "$downloads"
[[ -s "$downloads/java.tar.gz" ]] || curl -fsSL https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jre/hotspot/normal/eclipse -o "$downloads/java.tar.gz"
cp "$downloads/java.tar.gz" "$runtime/java.tar.gz"
tar -xzf "$runtime/java.tar.gz" -C "$runtime"
[[ -s "$downloads/nextflow.jar" ]] || curl -fsSL https://github.com/nextflow-io/nextflow/releases/download/v23.10.1/nextflow-23.10.1-all -o "$downloads/nextflow.jar"
cp "$downloads/nextflow.jar" "$runtime/nextflow.jar"
if [[ "${1:-}" != chimera_cohort_wiring_check.nf ]]; then
printf 'id\trole\trole_reasons\nsample_hap1\tvoter\tqualified\n' > "$runtime/quality.tsv"
printf 'old_name\tnew_name\tchromosome_member\ns\tchr1_1\tyes\n' > "$runtime/names.tsv"
touch "$runtime/NO_PAIRS" "$runtime/a.fa" "$runtime/reads.fq" "$runtime/sample.hap1.p_ctg.gfa" "$runtime/lib.tsv" "$runtime/a.agp" "$runtime/a.pairs" "$runtime/calls.tsv"
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
python3 - "$(dirname "$0")/${1:-chimera_channel_reproduction.nf}" "$runtime/reproduction.nf" "$repo" <<'PY'
import sys
from pathlib import Path
source = Path(sys.argv[1]).read_text()
if '__POST_CUT_QC__' in source:
    main = (Path(sys.argv[3]) / 'main.nf').read_text()
    start = main.index('    if (run_all_qc && chimera_on) {')
    end = main.index('    if (qc_on) ch_staged_assemblies', start)
    source = source.replace('__POST_CUT_QC__', main[start:end])
if '__FINAL_HIC_QC__' in source:
    main = (Path(sys.argv[3]) / 'main.nf').read_text()
    start = main.index('    // Final Hi-C evidence belongs to QC')
    end = main.index('    // String-id view', start)
    source = source.replace('__FINAL_HIC_QC__', main[start:end])
Path(sys.argv[2]).write_text(source.replace('__REPO__',sys.argv[3]))
PY
cd "$runtime"
extra=()
if [[ "${1:-}" == chimera_cohort_wiring_check.nf || "${1:-}" == teloclip_scope_wiring_check.nf ]]; then
    [[ "${3:-}" == resume ]] || python3 "$repo/tests/setup_misassembly_integration.py" "${2:-}"
    [[ "${1:-}" != teloclip_scope_wiring_check.nf ]] || python3 "$repo/tests/setup_teloclip_scope_integration.py"
    cp reproduction.nf "$repo/.misassembly-wiring.nf"
    trap 'rm -f "$repo/.misassembly-wiring.nf"' EXIT
    extra=(-c "$runtime/check.config")
    [[ "${3:-}" == resume ]] && extra+=(-resume)
    script="$repo/.misassembly-wiring.nf"
    # Keep the Nextflow cache on Linux: its database locks require native filesystem semantics.
    launch="$runtime/nxf-launch"
    mkdir -p "$repo/comparisons/misassembly-validation"
    mkdir -p "$launch"
    cd "$launch"
elif [[ "${1:-}" == deletion_replay_checks.nf ]]; then
    python3 "$repo/tests/setup_deletion_replay_checks.py" "$runtime"
    cp reproduction.nf "$repo/.deletion-replay-checks.nf"
    trap 'rm -f "$repo/.deletion-replay-checks.nf"' EXIT
    script="$repo/.deletion-replay-checks.nf"
    extra=(-c "$runtime/deletion-check.config")
    mkdir -p "$runtime/deletion-launch"
    cd "$runtime/deletion-launch"
elif [[ "${1:-}" == assignment_replay_checks.nf ]]; then
    python3 "$repo/tests/setup_assignment_replay_checks.py" "$runtime"
    cp reproduction.nf "$repo/.assignment-replay-checks.nf"
    trap 'rm -f "$repo/.assignment-replay-checks.nf"' EXIT
    script="$repo/.assignment-replay-checks.nf"
    extra=(-c "$runtime/assignment-check.config")
    mkdir -p "$runtime/assignment-launch"
    cd "$runtime/assignment-launch"
elif [[ "${1:-}" == core_contract_checks.nf ]]; then
    python3 "$repo/tests/setup_core_contract_checks.py" "$runtime" "${2:-success}"
    cp reproduction.nf "$repo/.core-contract-checks.nf"
    trap 'rm -f "$repo/.core-contract-checks.nf"' EXIT
    script="$repo/.core-contract-checks.nf"
    extra=(-c "$runtime/core-check.config")
    mkdir -p "$runtime/core-launch"
    cd "$runtime/core-launch"
elif [[ "${1:-}" == qc_post_cut_wiring_check.nf || "${1:-}" == final_hic_qc_wiring_check.nf ]]; then
    touch "$runtime/b.fa"
    cp reproduction.nf "$repo/.qc-wiring.nf"
    trap 'rm -f "$repo/.qc-wiring.nf"' EXIT
    script="$repo/.qc-wiring.nf"
    cat > "$runtime/qc-check.config" <<'EOF'
conda.enabled = false
docker.enabled = false
singularity.enabled = false
trace.enabled = false
report.enabled = false
timeline.enabled = false
dag.enabled = false
params.outdir = '/tmp/chimera-nxf-runtime/qc-results'
process {
    executor = 'local'
    withName: '.*' { cpus = 1; memory = '1GB'; scratch = false; errorStrategy = 'terminate' }
}
EOF
    if [[ "${1:-}" == final_hic_qc_wiring_check.nf ]]; then
        touch "$runtime/one_R1.fq" "$runtime/one_R2.fq" "$runtime/two_R1.fq" "$runtime/two_R2.fq"
        cat >> "$runtime/qc-check.config" <<EOF
params.qc_mode = '${2:-final_only}'
params.run_final_contact_maps = ${3:-true}
params.run_hic_balance = true
EOF
    fi
    extra=(-stub-run -c "$runtime/qc-check.config")
    mkdir -p "$runtime/qc-launch"
    cd "$runtime/qc-launch"
else
    script=reproduction.nf
fi
runtime_java=$(find "$runtime" -path '*/bin/java' -type f -print -quit)
# Match Nextflow 23.10.1's launcher permissions on Java 17. Missing these can
# silently prevent its asynchronous Kryo cache writer from persisting task entries.
java_opens=()
for package in java.lang java.io java.nio java.net java.util java.util.concurrent.locks \
    java.util.concurrent.atomic java.nio.file.spi sun.nio.ch sun.nio.fs \
    sun.net.www.protocol.http sun.net.www.protocol.https sun.net.www.protocol.ftp \
    sun.net.www.protocol.file jdk.internal.misc jdk.internal.vm java.util.regex; do
    java_opens+=("--add-opens=java.base/$package=ALL-UNNAMED")
done
run_nextflow() {
# The executable JAR can start a child JVM; pass options through to that JVM too.
JAVA_TOOL_OPTIONS="${java_opens[*]} --add-exports=java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED --add-opens=java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED" \
"$runtime_java" --add-exports java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    "${java_opens[@]}" -jar "$runtime/nextflow.jar" run "$script" "${extra[@]}" -cache true -ansi-log false
}
if ! run_nextflow; then
    find work -name .command.err -mmin -10 -exec sh -c 'echo "$1"; tail -20 "$1"' sh {} \;
    exit 1
fi
if [[ "${1:-}" == deletion_replay_checks.nf ]]; then
    deletion_out="$runtime/deletion-results/pangenome/1/deletion_filter"
    if [[ "${2:-}" == regional_patch_resume || "${2:-}" == regional_trace_resume ]]; then
        extra+=(-resume --pangenome_regional_patch_replay true)
        run_nextflow | tee regional-resume.log
        grep -q 'Cached process > PANGENOME_DELETION_DIAGNOSTIC:PANGENOME_DELETION_PREP' regional-resume.log
        grep -q 'Cached process > PANGENOME_DELETION_DIAGNOSTIC:PANGENOME_DELETION_FILTER' regional-resume.log
        test -s "$runtime/deletion-results/pangenome/1/deletion_filter/regional_assignment_report.md"
        test -s "$runtime/deletion-results/pangenome/1/deletion_filter/regional_coordinate_audit.tsv"
        test -s "$runtime/deletion-results/pangenome/tools/regional_patch/native_regressions.tsv"
        python3 - "$runtime/deletion-results/pangenome/1/deletion_filter" <<'PY'
import csv
from pathlib import Path
import sys
root=Path(sys.argv[1])
with (root/'regional_assignment_comparison.tsv').open() as f:
    rows={r['policy']:r for r in csv.DictReader(f,delimiter='\t')}
assert int(rows['native']['fixed_regional_bp'])==4000000,rows
for policy in ('chromosome_isolated','disabled'):
    assert int(rows[policy]['fixed_regional_bp'])==8000000,rows
    assert set(rows[policy]['fixed_chromosomes'].split(','))=={'chrA','chrB'},rows
with (root/'regional_coordinate_audit.tsv').open() as f:
    audits=list(csv.DictReader(f,delimiter='\t'))
assert len(audits)==3 and all(r['status']=='PASS' and int(r['duplicate_or_overlapping_bp'])==0 for r in audits)
print('PASS: both chromosome arms retained by the real corrected tool without duplicate source ownership')
PY
        mkdir -p "$repo/comparisons/core-graph-regional-validation"
        cp regional-resume.log "$repo/comparisons/core-graph-regional-validation/nextflow-resume.log"
        cp "$runtime/deletion-results/pangenome/1/deletion_filter/"regional_assignment* "$repo/comparisons/core-graph-regional-validation/"
        cp "$runtime/deletion-results/pangenome/tools/regional_patch/native_regressions.tsv" "$repo/comparisons/core-graph-regional-validation/"
        echo 'PASS: patched native regional assignment runs; preparation and deletion filters resume from cache'
        if [[ "${2:-}" == regional_trace_resume ]]; then
            extra+=(--pangenome_regional_gate_trace true)
            run_nextflow | tee regional-trace-resume.log
            grep -q 'Cached process > PANGENOME_DELETION_DIAGNOSTIC:PANGENOME_DELETION_FILTER' regional-trace-resume.log
            grep -q 'Cached process > .*PANGENOME_REGIONAL_BUILD' regional-trace-resume.log
            grep -q 'Cached process > .*PANGENOME_REGIONAL_ASSIGN_FIXED' regional-trace-resume.log
            test -s "$deletion_out/regional_gate_report.md"
            test -s "$deletion_out/trace_chromosome_isolated/fixed_events.tsv.gz"
            test -s "$runtime/deletion-results/pangenome/tools/regional_trace/native_trace.diff"
            test -s "$runtime/deletion-results/pangenome/tools/regional_trace/gate_regressions.tsv"
            python3 - "$runtime/deletion-results/pangenome/1/deletion_filter" <<'PY'
import csv
from pathlib import Path
import sys
root=Path(sys.argv[1])
with (root/'regional_observation_equivalence.tsv').open() as f:checks=list(csv.DictReader(f,delimiter='\t'))
assert len(checks)==6 and all(r['byte_identical']=='yes' for r in checks)
with (root/'regional_gate_summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert any(r['policy']=='chromosome_isolated' and r['variant']=='fixed' and r['role']=='secondary'
           and int(r['locally_accepted_runs'])==1 for r in rows),rows
assert {r['chromosome'] for r in rows if r['policy']=='chromosome_isolated' and r['variant']=='fixed'
        and int(r['final_assigned_span_bp'])==4000000}=={'chrA','chrB'}
print('PASS: six byte-equivalent observations, real regional gates, unchanged cached assignment/filter tasks')
PY
            cp regional-trace-resume.log "$repo/comparisons/core-graph-regional-validation/nextflow-trace-resume.log"
            cp "$deletion_out/"regional_gate* "$repo/comparisons/core-graph-regional-validation/"
        fi
    fi
    deletion_out="$runtime/deletion-results/pangenome/1/deletion_filter"
    test -s "$deletion_out/deletion_filter_report.md"
    test -s "$deletion_out/assignment_comparison.tsv"
    test -s "$deletion_out/review_chromosome_isolated/lanes.tsv"
    for policy in native chromosome_isolated disabled; do
        test -s "$deletion_out/review_$policy/filter_coverage.tsv"
        for assignment in baseline regional_1Mb; do
            test -s "$deletion_out/assignment_${policy}_${assignment}/surviving_intervals.tsv.gz"
        done
    done
    if find "$deletion_out" -type f \( -name '*.paf' -o -name '*.gfa' -o -name '*.vg' -o -name 'node_owners.tsv.gz' \) -print -quit | grep -q .; then
        echo 'Large deletion intermediates must not be published' >&2
        exit 1
    fi
    echo 'PASS: compact deletion evidence published; large intermediates stay in cache'
fi
if [[ "${1:-}" == assignment_replay_checks.nf ]]; then
    assignment_out="$runtime/assignment-results/pangenome/1/assignment"
    test -s "$assignment_out/assignment_report.md"
    test -s "$assignment_out/assignment_comparison.tsv"
    test -s "$assignment_out/coverage_trace.tsv"
    test -s "$assignment_out/coverage_bins.tsv"
    test -s "$assignment_out/coverage_changes.tsv.gz"
    test -s "$assignment_out/coverage_intervals.tsv.gz"
    test -s "$assignment_out/preparation/coverage_intervals.tsv.gz"
    test -s "$assignment_out/preparation/effective_config.xml"
    for assignment_case in baseline whole_relaxed regional_1Mb regional_3Mb; do
        test -s "$assignment_out/replay_$assignment_case/assignment_log.txt"
        test -s "$assignment_out/replay_$assignment_case/regions.tsv"
        test -s "$assignment_out/replay_$assignment_case/decisions.tsv"
        test -s "$assignment_out/replay_$assignment_case/surviving_intervals.tsv.gz"
    done
    if find "$assignment_out" -type f \( -name '*.gfa' -o -name '*.paf' \) -print -quit | grep -q .; then
        echo 'Large replay intermediates must not be published' >&2
        exit 1
    fi
    echo 'PASS: assignment logs/configuration published; large intermediates remain in cache'
fi
if [[ "${1:-}" == final_hic_qc_wiring_check.nf && "${2:-final_only}" != none ]]; then
    for checkpoint in final_raw_map final_filtered; do
        count=$(find "$runtime/qc-results/qc/hic/map/$checkpoint" -name '*.tsv' -type f | wc -l)
        [[ "$count" -eq 4 ]] || { echo "Expected four published $checkpoint tables; found $count" >&2; exit 1; }
    done
    echo 'PASS: final Hi-C standalone metric tables are published'
fi
if [[ "${3:-}" == cache ]]; then
    cp .nextflow.log "$repo/comparisons/misassembly-validation/first-nextflow.log"
    find work -name .exitcode -exec sh -c 'echo "$1"; cat "$1"' sh {} \; > "$repo/comparisons/misassembly-validation/cache-files.txt"
    extra+=(-resume)
    run_nextflow | tee cache-resume.log
    cp .nextflow.log "$repo/comparisons/misassembly-validation/resume-nextflow.log"
    grep -q 'Cached process > CHIMERA:MISASSEMBLY_ALIGN' cache-resume.log
    grep -q 'Cached process > CHIMERA:BREAK_CHIMERAS' cache-resume.log
    grep -q 'Cached process > CHIMERA:CHIMERA_REASSIGN_SPECIES' cache-resume.log
    grep -q 'Cached process > FINALIZE_ASSEMBLY' cache-resume.log
    grep -q 'Cached process > FINAL_READ_EVIDENCE:MAP_HIFI_FINAL' cache-resume.log
    echo 'PASS: discovery, evidence, cutting, reassignment, finishing and final read mapping resume from cache'
fi
if [[ "${1:-}" == core_contract_checks.nf && "${3:-}" == resume_build ]]; then
    extra+=(-resume)
    run_nextflow | tee core-resume.log
    grep -q 'Cached process > PANGENOME_CONSTRUCTION:CACTUS_PANGENOME' core-resume.log
    cp core-resume.log "$repo/comparisons/core-graph-resume-check.log"
    echo 'PASS: completed graph construction uses ordinary Nextflow resume'
fi
