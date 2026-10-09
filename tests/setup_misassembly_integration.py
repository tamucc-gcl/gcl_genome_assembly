"""External-tool doubles for channel/CLI integration, not biological validation."""
from pathlib import Path
import sys,hashlib
root=Path('/tmp/chimera-nxf-runtime');root.mkdir(exist_ok=True)
evidence='evidence' in sys.argv or 'manual' in sys.argv
manual='manual' in sys.argv
(root/'a.fa').write_text('>s\n'+'A'*3000+('N'*100 if evidence else '')+'C'*3000+'\n')
(root/'b.fa').write_text('>x\n'+'A'*3000+'\n>y\n'+'C'*3000+'\n')
(root/'names.tsv').write_text('old_name\tnew_name\tchromosome_member\torient\tlength\tclass\ns\tchr1+chr2\tyes\t+\t'+str(6100 if evidence else 6000)+'\tchromosome\n')
(root/'b.names.tsv').write_text('old_name\tnew_name\tchromosome_member\torient\tlength\tclass\nx\tchr1\tyes\t+\t3000\tchromosome\ny\tchr2\tyes\t+\t3000\tchromosome\n')
(root/'quality.tsv').write_text('id\trole\trole_reasons\nsample_hap1\tvoter\tqualified\nreference_hap1\tvoter\tqualified\n')
bin=root/'bin';bin.mkdir(exist_ok=True)
shim=bin/'tools.py'
shim.write_text((Path(__file__).parent/'misassembly_tool_double.py').read_text())
shim.chmod(0o755)
for name in ('samtools','minimap2','tidk','date'):
    path=bin/name
    if path.is_symlink():path.unlink()
    path.symlink_to(shim)
(root/'check.config').write_text('''params.outdir='/tmp/chimera-nxf-runtime/new-published'
params.publish_dir_mode='copy'
params.chimera_detect=true
params.chimera_break=false
params.chimera_plots=false
params.chimera_hifi_context=false
params.chimera_min_arm_bp=1000
params.fixture_manual=false
params.fixture_evidence=false
params.harmonize_scaffold_names=true
conda.enabled=false
trace.enabled=true
trace.file='/tmp/chimera-nxf-runtime/integration-trace.tsv'
timeline.enabled=false
report.enabled=false
dag.enabled=false
process.executor='local'
process.cpus=1
process.memory='1 GB'
process.time='10m'
process.scratch=false
process.beforeScript=''
process {
  withName: '.*' {
    cpus=1
    memory='1 GB'
    errorStrategy='terminate'
}
}
env.PATH='/tmp/chimera-nxf-runtime/bin:/usr/local/bin:/usr/bin:/bin'
''')

(root/'reads.fq').write_text('@r\nA\n+\nI\n')
(root/'a.agp').write_text('s\t1\t3000\t1\tW\tx\t1\t3000\t+\ns\t3001\t3100\t2\tN\t100\tscaffold\tyes\tproximity_ligation\ns\t3101\t6100\t3\tW\ty\t1\t3000\t+\n')
(root/'lib.tsv').write_text('library_id\tqname_prefix\nA\tA_\nB\tB_\n')
(root/'a.pairs').write_text(''.join(f'{lib}_{chrom}_{i}\t{chrom}\t1\t{chrom}\t101\t+\t-\tUU\n' for lib in ('A','B') for chrom in ('x','y') for i in range(100)))
if evidence:
    with (root/'check.config').open('a') as handle:handle.write('\nparams.fixture_evidence=true\nparams.chimera_hifi_context=true\n')
if manual:
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
    from misassembly_report import FIELDS
    from misassembly_core import write_table
    write_table(root/'reviewed.tsv',[dict(selected='YES',review_disposition='CUT',id='manual01',assembly='sample_hap1',scaffold='s',cut_bp=3100,action='UNJOIN_UNSUPPORTED',gap_start=3000,gap_end=3100,localization_status='available_gap',reviewer='integration-test',reason='Software test of explicit manual override.',coordinate_stage='pre_finishing',assessment_sha256=hashlib.sha256((root/'a.fa').read_bytes()).hexdigest(),decision_source='review',evidence_packet_id='sample_hap1')],FIELDS)
    with (root/'check.config').open('a') as handle:handle.write("\nparams.chimera_break='/tmp/chimera-nxf-runtime/reviewed.tsv'\nparams.chimera_min_piece_bp=1000\nparams.fixture_manual=true\n")
