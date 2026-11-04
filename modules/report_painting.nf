process REPORT_PAINTING {
    tag "${meta.id}_vs_${reference_id}"
    label 'comparison_alignment'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/reports/chromosome_painting/${meta.taxid}", mode: params.publish_dir_mode
    input:
    tuple val(meta),val(reference_id),path(reference,stageAs:'reference/*'),path(ref_fai,stageAs:'reference/*'),path(query,stageAs:'query/*'),path(query_fai,stageAs:'query/*')
    path helper
    output:
    tuple val(meta),val(reference_id),path("${meta.id}.painting.tsv"),path("${meta.id}.painting.json"),emit: blocks
    tuple val(meta),val(reference_id),path("${meta.id}.painting.paf.gz"),path("${meta.id}.painting.log"),emit: evidence
    path "${meta.id}.painting.versions.tsv",emit: versions
    script:
    def self=meta.id==reference_id
    def preset=params.pairwise_alignment_preset ?: 'asm5'
    def mq=params.pairwise_alignment_min_mapq != null ? params.pairwise_alignment_min_mapq : 5
    def minlen=params.pairwise_alignment_min_aln_bp != null ? params.pairwise_alignment_min_aln_bp : 150000
    if (!(preset in ['asm5','asm10','asm20'])) error 'Unsupported painting alignment preset'
    """
    set -euo pipefail
    if ${self}; then
        printf '' | gzip > '${meta.id}.painting.paf.gz'
        printf 'Selected reference: identity labels, not an independent alignment.\\n' > '${meta.id}.painting.log'
    else
        REF_BASES=\$(awk '{n+=\$2} END {printf "%.0f",n}' '${ref_fai}')
        minimap2 -t ${task.cpus} -x ${preset} -I \$((REF_BASES+1)) --secondary=no '${reference}' '${query}' \\
          2> '${meta.id}.painting.log' | gzip > '${meta.id}.painting.paf.gz'
    fi
    python3 ${helper} --query-fai '${query_fai}' --ref-fai '${ref_fai}' \\
      --paf '${meta.id}.painting.paf.gz' --output '${meta.id}.painting.tsv' \\
      --metadata '${meta.id}.painting.json' --reference-id '${reference_id}' --query-id '${meta.id}' \\
      --min-mapq ${mq} --min-length ${minlen} ${self ? '--self-reference' : ''}
    printf 'minimap2\\t%s\\n' "\$(minimap2 --version)" > '${meta.id}.painting.versions.tsv'
    """
}
