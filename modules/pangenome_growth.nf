/*
========================================================================================
    PANGENOME GROWTH / OPENNESS MODULE  (workstream E)
========================================================================================
    Repo location: modules/pangenome_growth.nf

    Openness / scaling analysis with panacus, computed from the FINISHED clip GFA (no
    progressive rebuild). Grouping is by haplotype via PanSN, so the graph's many contig
    paths collapse to one unit per haploid genome. This process generates the DATA only;
    the styled curves, Heaps'-law fit (gamma), and confidence band are produced by the R
    plotting step (workstream D) from the histgrowth TSV. (panacus-visualize is deprecated
    and slated for removal, so it is intentionally not used.)

    Produces:
      - histgrowth TSV : pangenome (union), core, and soft-core growth curves. panacus
                         computes the EXACT expected growth (analytically averaged over all
                         orderings), i.e. the permutation-averaged curve by construction.
      - hist TSV       : coverage histogram (bp at each haplotype-coverage level) -> the
                         core / accessory / private partition.
      - core_accessory : parsed core/accessory/private bp for the report.

    Input : tuple(taxid, gfa)          (CACTUS_PANGENOME.out.gfa -- the clip GFA .gfa.gz)
    Output: histgrowth / hist / core_accessory / versions
========================================================================================
*/

process PANGENOME_GROWTH {
    tag "${taxid}"
    label 'panacus'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(gfa)

    output:
    tuple val(taxid), path("${taxid}.histgrowth.${cnt}.tsv"), emit: histgrowth
    tuple val(taxid), path("${taxid}.hist.${cnt}.tsv"),       emit: hist
    tuple val(taxid), path("${taxid}.core_accessory.tsv"),    emit: core_accessory
    path("versions.tsv"),                                     emit: versions

    script:
    cnt = params.pangenome_growth_count ?: 'bp'
    def cov = params.pangenome_growth_coverage ?: '1,1,1'
    def quo = params.pangenome_growth_quorum   ?: '0,1,0.9'
    """
    set -euo pipefail

    # panacus needs an uncompressed, blunt GFA with P/W lines (cactus clip GFA qualifies)
    zcat ${gfa} > graph.gfa

    # growth + core curves (exact expected; grouped by haplotype via PanSN) --------------
    panacus histgrowth -c ${cnt} -l ${cov} -q ${quo} --groupby-haplotype -t ${task.cpus} \\
        graph.gfa > ${taxid}.histgrowth.${cnt}.tsv

    # coverage histogram (bp per haplotype-coverage level) ------------------------------
    panacus hist -c ${cnt} --groupby-haplotype -t ${task.cpus} \\
        graph.gfa > ${taxid}.hist.${cnt}.tsv

    # core / accessory / private partition from the coverage histogram ------------------
    # hist TSV: skip comment/header lines; col1 = coverage level (1..N), col2 = count.
    awk 'BEGIN{FS=OFS="\\t"}
         /^#/ {next}
         \$1 ~ /^[0-9]+\$/ { cov[\$1]=\$2; if(\$1+0>N) N=\$1+0 }
         END{
           core=cov[N]+0; priv=cov[1]+0; acc=0
           for(k in cov){ if(k+0>1 && k+0<N) acc+=cov[k] }
           print "class","coverage_level","count"
           print "private",1,priv
           print "accessory","2.."(N-1),acc
           print "core",N,core
           print "n_haplotypes",N,N
         }' ${taxid}.hist.${cnt}.tsv > ${taxid}.core_accessory.tsv

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpanacus\\t%s\\n' "${task.process}" "\$(panacus --version 2>&1 | sed 's/panacus //')"
    } > versions.tsv

    rm -f graph.gfa
    """

    stub:
    cnt = params.pangenome_growth_count ?: 'bp'
    """
    printf '# panacus stub\\n' > ${taxid}.histgrowth.${cnt}.tsv
    printf '# panacus stub\\n' > ${taxid}.hist.${cnt}.tsv
    printf 'class\\tcoverage_level\\tcount\\n' > ${taxid}.core_accessory.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
