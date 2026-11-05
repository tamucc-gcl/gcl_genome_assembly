nextflow.enable.dsl=2
include { RUN_PROVENANCE } from '__REPO__/workflows/run_provenance.nf'
process VERSION_RECORDS {
    output:
    path '*/versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    mkdir -p bin tidk_record cooler_record
    printf '#!/bin/bash\\nprintf "tidk\\\\nversion\\\\n0.2.65\\\\n"\\n' > bin/tidk
    printf '#!/bin/bash\\nprintf "cooler, version 0.10.4\\\\n"\\n' > bin/cooler
    chmod +x bin/*
    export PATH="\$PWD/bin:\$PATH"
    cd tidk_record
    __TIDK_VERSION_STATEMENT__
    cd ../cooler_record
    __COOLER_VERSION_STATEMENT__
    cd ..
    python3 - <<'PY'
import csv
from pathlib import Path
for tool,expected in [('tidk','0.2.65'),('cooler','0.10.4')]:
    rows=list(csv.reader((Path(tool+'_record')/'versions.tsv').open(),delimiter='\\t'))
    assert len(rows)==1 and len(rows[0])==2,rows
    assert rows[0][0]==tool and expected in rows[0][1],rows
PY
    """
}
workflow {
    VERSION_RECORDS()
    RUN_PROVENANCE(VERSION_RECORDS.out.versions.flatten())
    RUN_PROVENANCE.out.files.toList().view { files ->
        def all=files.flatten()
        def issues=all.find { it.name=='software_version_issues.tsv' }
        assert issues.readLines().size()==1
        'PASS: actual producer statements label single/multiline versions; real collector reports zero issues'
    }
}
