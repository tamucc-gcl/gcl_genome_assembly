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
Path(sys.argv[2]).write_text(Path(sys.argv[1]).read_text().replace('__REPO__',sys.argv[3]))
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
    launch="$repo/comparisons/misassembly-validation/nxf-launch"
    mkdir -p "$launch"
    cd "$launch"
else
    script=reproduction.nf
fi
runtime_java=$(find "$runtime" -path '*/bin/java' -type f -print -quit)
run_nextflow() {
"$runtime_java" --add-exports java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED -jar "$runtime/nextflow.jar" run "$script" "${extra[@]}" -ansi-log false
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
    echo 'PASS: discovery, evidence, cutting and reassignment resume from cache'
fi
