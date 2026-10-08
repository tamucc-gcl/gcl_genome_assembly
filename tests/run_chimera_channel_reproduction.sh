#!/usr/bin/env bash
set -euo pipefail
runtime=/tmp/chimera-nxf-runtime
mkdir -p "$runtime"
downloads="$(dirname "$0")/../comparisons/nxf-runtime"
mkdir -p "$downloads"
[[ -s "$downloads/java.tar.gz" ]] || curl -fsSL https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jre/hotspot/normal/eclipse -o "$downloads/java.tar.gz"
cp "$downloads/java.tar.gz" "$runtime/java.tar.gz"
tar -xzf "$runtime/java.tar.gz" -C "$runtime"
[[ -s "$downloads/nextflow.jar" ]] || curl -fsSL https://github.com/nextflow-io/nextflow/releases/download/v23.10.1/nextflow-23.10.1-all -o "$downloads/nextflow.jar"
cp "$downloads/nextflow.jar" "$runtime/nextflow.jar"
printf 'id\trole\trole_reasons\nsample_hap1\tvoter\tqualified\n' > "$runtime/quality.tsv"
printf 'old_name\tnew_name\tchromosome_member\ns\tchr1_1\tyes\n' > "$runtime/names.tsv"
touch "$runtime/NO_PAIRS" "$runtime/a.fa" "$runtime/reads.fq" "$runtime/sample.hap1.p_ctg.gfa" "$runtime/lib.tsv" "$runtime/a.agp" "$runtime/a.pairs" "$runtime/calls.tsv"
cp "$(dirname "$0")/${1:-chimera_channel_reproduction.nf}" "$runtime/reproduction.nf"
cd "$runtime"
runtime_java=$(find "$runtime" -path '*/bin/java' -type f -print -quit)
"$runtime_java" --add-exports java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.base/java.lang=ALL-UNNAMED -jar "$runtime/nextflow.jar" run reproduction.nf
