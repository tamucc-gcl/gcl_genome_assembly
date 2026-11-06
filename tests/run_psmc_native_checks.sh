#!/usr/bin/env bash
# Native tools are supplied by an explicitly prepared test environment.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
runtime="${PSMC_TEST_ROOT:-/tmp/gcl-psmc-native}"
bin="${PSMC_TEST_BIN:?Set PSMC_TEST_BIN to the pinned PSMC environment bin directory}"
rbin="${PSMC_R_BIN:?Set PSMC_R_BIN to the tested R plotting environment bin directory}"
nxf="${PSMC_NEXTFLOW_JAR:-/tmp/chimera-nxf-runtime/nextflow.jar}"
java="${PSMC_JAVA:-$(find /tmp/chimera-nxf-runtime -path '*/bin/java' -type f -print -quit)}"
export PATH="$bin:$PATH"
mkdir -p "$runtime"
if [[ "${1:-}" != resume && "${1:-}" != scale ]]; then
    python3 "$repo/tests/setup_psmc_checks.py" "$runtime/fixtures"
fi
case="${PSMC_TEST_CASE:-enabled}"
if [[ "$case" == all_skipped ]]; then
    printf 'sample\tmode\nhaploid_ind\tauto\n' > "$runtime/fixtures/manifest.tsv"
fi
if [[ "$case" == paired_hifi ]]; then
    printf 'sample\tmode\nassembly_ind\thifi\n' > "$runtime/fixtures/manifest.tsv"
fi
python3 - "$repo" <<'PY'
from pathlib import Path
import sys
repo=Path(sys.argv[1])
(repo/'.psmc-native-check.nf').write_text((repo/'tests/psmc_checks.nf').read_text().replace('__REPO__',str(repo)))
PY
trap 'rm -f "$repo/.psmc-native-check.nf"' EXIT
cat > "$runtime/local.config" <<EOF
conda.enabled=false
singularity.enabled=false
docker.enabled=false
report.enabled=false
timeline.enabled=false
dag.enabled=false
trace.enabled=false
params.run_psmc=$([[ "$case" == disabled ]] && echo false || echo true)
params.psmc_test_case='$case'
params.psmc_fixture='$runtime/fixtures'
params.psmc_manifest='$runtime/fixtures/manifest.tsv'
params.psmc_min_contig_bp=100000
params.psmc_min_good_bases=10000
params.psmc_min_callable_bp=10000
params.psmc_min_hets=${PSMC_TEST_MIN_HETS:-20}
params.psmc_pattern='2+2+2'
params.psmc_iterations=3
params.psmc_plot_iteration=2
params.psmc_bootstraps=${PSMC_TEST_BOOTSTRAPS:-2}
params.psmc_bootstrap_block_bp=50000
params.psmc_bootstrap_max_forks=2
params.outdir='$runtime/results'
params.publish_dir_mode='copy'
process {
 executor='local'
 errorStrategy='terminate'
 scratch=false
 withName: '.*' { cpus=2; memory='2 GB' }
 withName: '.*PSMC_RENDER' { beforeScript='export PATH=$rbin:$bin:\$PATH' }
 withName: '.*PSMC_INDEX' { beforeScript='export PATH=$rbin:$bin:\$PATH' }
}
EOF
extra=()
if [[ "${1:-}" == resume || "${1:-}" == scale ]]; then extra+=(-resume); fi
if [[ "${1:-}" == scale ]]; then
    cat >> "$runtime/local.config" <<EOF
params.psmc_mutation_rate=1e-8
params.psmc_generation_time=0.5
EOF
fi
opens=()
for p in java.lang java.io java.nio java.net java.util java.util.concurrent.locks java.util.concurrent.atomic java.nio.file.spi sun.nio.ch sun.nio.fs sun.net.www.protocol.http sun.net.www.protocol.https sun.net.www.protocol.ftp sun.net.www.protocol.file jdk.internal.misc jdk.internal.vm java.util.regex; do
    opens+=("--add-opens=java.base/$p=ALL-UNNAMED")
done
cd "$runtime"
JAVA_TOOL_OPTIONS="${opens[*]} --add-exports=java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED --add-opens=java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED" \
    "$java" "${opens[@]}" --add-exports java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED \
    --add-opens java.management/com.sun.jmx.mbeanserver=ALL-UNNAMED -jar "$nxf" \
    -log "$runtime/nextflow-${1:-initial}.log" run "$repo/.psmc-native-check.nf" \
    -c "$runtime/local.config" -ansi-log false "${extra[@]}"
