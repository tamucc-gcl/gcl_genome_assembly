process PAIRWISE_ALIGNMENT {
    tag "${id1}_vs_${id2}"
    label 'comparison_alignment'
    errorStrategy { task.exitStatus in [65,126,127] ? 'terminate' : 'retry' }
    maxRetries 1
    publishDir "${params.outdir}/pairwise_alignments", mode: params.publish_dir_mode
    input:
    tuple val(id1), path(fa1,stageAs:'ref/*'), path(fai1,stageAs:'ref/*'), val(id2), path(fa2,stageAs:'query/*'), path(fai2,stageAs:'query/*')
    path helper
    output:
    tuple val(id1),val(id2),path("${id1}_vs_${id2}.paf.gz"), emit: paf
    tuple val(id1),val(id2),path("${id1}_vs_${id2}.qc.tsv"), emit: qc
    tuple val(id1),val(id2),path("${id1}_vs_${id2}.log"), emit: log
    path "${id1}_vs_${id2}.versions.tsv", emit: versions
    script:
    def preset=params.pairwise_alignment_preset ?: 'asm5'
    def mq=params.pairwise_alignment_min_mapq != null ? params.pairwise_alignment_min_mapq : 5
    def minlen=params.pairwise_alignment_min_aln_bp != null ? params.pairwise_alignment_min_aln_bp : 150000
    if (!(preset in ['asm5','asm10','asm20'])) error "Unsupported assembly comparison preset '${preset}'"
    """
    set -euo pipefail
    # Multi-part indexing invalidates MAPQ; keep this reference in one index.
    REF_BASES=\$(awk '{n+=\$2} END {printf "%.0f",n}' '${fai1}')
    INDEX_BASES=\$((REF_BASES+1))
    minimap2 -t ${task.cpus} -x ${preset} -I \${INDEX_BASES} -c --eqx --secondary=no '${fa1}' '${fa2}' \\
        2> '${id1}_vs_${id2}.log' \\
      | python3 ${helper} filter --output '${id1}_vs_${id2}.paf.gz' \\
        --qc '${id1}_vs_${id2}.qc.tsv' --ref '${id1}' --query '${id2}' \\
        --preset '${preset}' --min-mapq ${mq} --min-length ${minlen}
    printf 'minimap2\\t%s\\n' "\$(minimap2 --version)" > '${id1}_vs_${id2}.versions.tsv'
    """
}

process PAIRWISE_PLOTS {
    tag "${id1}_vs_${id2}"
    label 'pairwise_plotting'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/pairwise_alignments", mode: params.publish_dir_mode
    input:
    tuple val(id1),val(id2),path(paf),path(fai1,stageAs:'ref/*'),path(fai2,stageAs:'query/*')
    path helpers
    output:
    tuple val(id1),val(id2),path("${id1}_vs_${id2}_dotplot.png"), emit: dotplot
    tuple val(id1),val(id2),path("${id1}_vs_${id2}_riparian.png"), emit: riparian
    tuple val(id1),val(id2),path("${id1}_vs_${id2}.plot_status.tsv"), emit: status
    path "${id1}_vs_${id2}.plot_versions.tsv", emit: versions
    script:
    def nz={ key,fallback -> params[key]!=null ? params[key] : fallback }
    """
    set -euo pipefail
    python3 plot_pairwise.py --paf '${paf}' --ref-fai '${fai1}' --query-fai '${fai2}' \\
      --ref '${id1}' --query '${id2}' --prefix '${id1}_vs_${id2}' \\
      --dot-min-seq ${nz('pairwise_dotplot_min_seq_bp',1000000)} \\
      --ribbon-min-seq ${nz('riparian_min_seq_bp',1000000)} --ribbon-min-aln ${nz('riparian_min_aln_bp',150000)} \\
      --ribbon-alpha ${nz('riparian_alpha',0.45)} \\
      --dot-width ${nz('pairwise_dotplot_width',10)} --dot-height ${nz('pairwise_dotplot_height',10)} \\
      --ribbon-width ${nz('riparian_width',14)} --ribbon-height ${nz('riparian_height',6)}
    python3 -c 'import matplotlib,sys; print("Python\\t"+sys.version.split()[0]); print("matplotlib\\t"+matplotlib.__version__)' > '${id1}_vs_${id2}.plot_versions.tsv'
    """
}

process COLLECT_PAIRWISE_RESULTS {
    tag 'collect_pairwise'
    cpus 1
    memory '1 GB'
    time '10m'
    conda 'conda-forge::python=3.11'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/pairwise_alignments", mode: params.publish_dir_mode
    input:
    path qc_files, stageAs: 'qc_files/*'
    path helper
    output:
    path 'pairwise_alignment_summary.tsv', emit: summary
    script:
    """
    python3 ${helper} collect --directory qc_files --output pairwise_alignment_summary.tsv
    """
}

process COMPARISON_PLAN {
    cpus 1
    memory '1 GB'
    time '10m'
    conda 'conda-forge::python=3.11'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/pairwise_alignments", mode: params.publish_dir_mode
    input:
    val rows
    path indexes, stageAs: 'indexes/*'
    val sources
    output:
    path 'comparison_plan.tsv', emit: plan
    path 'final_input_inventory.tsv', emit: inventory
    path 'final_indexes', emit: indexes
    script:
    def table=['taxid\tref_sample\tquery_sample\tref_id\tquery_id\tref_fasta\tquery_fasta\tcoordinate_stage\tmode']+rows.collect { it.join('\t') }
    def payload=groovy.json.JsonOutput.toJson(sources)
    """
    cat > comparison_plan.tsv <<'PLAN'
${table.join('\n')}
PLAN
    cat > sources.json <<'SOURCES'
${payload}
SOURCES
    python3 - <<'PY'
import csv,hashlib,json,shutil
from pathlib import Path
Path('final_indexes').mkdir()
fields=['id','sample','taxid','fasta','fai','coordinate_stage','fai_sha256','sequences','bases']
with open('final_input_inventory.tsv','w',newline='') as out:
    writer=csv.DictWriter(out,fields,delimiter='\\t');writer.writeheader()
    for source in json.loads(Path('sources.json').read_text()):
        index=Path('indexes')/source.pop('index_name')
        content=index.read_bytes();entries=[line.split('\\t') for line in content.decode().splitlines()]
        source.update(coordinate_stage='final',fai_sha256=hashlib.sha256(content).hexdigest(),
                      sequences=len(entries),bases=sum(int(row[1]) for row in entries))
        writer.writerow(source)
        shutil.copyfile(index,Path('final_indexes')/(source['id']+'.fasta.fai'))
PY
    """
}
