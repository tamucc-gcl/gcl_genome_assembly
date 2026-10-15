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
Path(sys.argv[2]).write_text(source.replace('__REPO__',sys.argv[3]))
PY
cd "$runtime"
extra=()
if [[ "${1:-}" == chimera_cohort_wiring_check.nf ]]; then
    [[ "${3:-}" == resume ]] || python3 "$repo/tests/setup_misassembly_integration.py" "${2:-}"
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
elif [[ "${1:-}" == qc_post_cut_wiring_check.nf ]]; then
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
