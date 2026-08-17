/*
========================================================================================
    PANGENOME MANIFEST MODULE
========================================================================================
    Repo location: modules/pangenome_manifest.nf

    Write pangenome_manifest.tsv (role, file, graph, label) enumerating the downstream-consumable
    products, published alongside them in the graph's output directory. This is the reuse
    seam: an external pipeline (short-read mapping / popgen) points at pangenome/<label>/,
    reads the manifest, and resolves files by role instead of hard-coding names. Filenames
    are basenames (relative to the manifest's own directory), so the bundle is portable.

    Read-mapping indexes (.dist/.min) are intentionally NOT here -- how the graph is indexed
    depends on the consuming project; it builds them from the GBZ + .hapl at map time.

    Input : tuple(taxid, gbz, gfa, og, snarls, hapl, vcf, vcf_tbi, ref_fa, ref_fai)
    Output: manifest
========================================================================================
*/

process PANGENOME_MANIFEST {
    tag "${taxid}"
    cpus 1
    memory '1 GB'
    time '10m'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    // full-graph products arrive as val NAMES (not paths): they are already declared and
    // published through CACTUS_PANGENOME's `all` catch-all, so picking them by name in the
    // workflow avoids adding output declarations to that (very expensive) process. An absent
    // product is passed as '' and its row is omitted.
    tuple val(taxid), path(gbz), path(gfa), path(og), path(snarls), path(hapl),
          path(vcf), path(vcf_tbi), path(ref_fa), path(ref_fai),
          val(full_gbz), val(full_gfa), val(full_og), val(full_snarls)

    output:
    tuple val(taxid), path("pangenome_manifest.tsv"), emit: manifest

    script:
    """
    {
      printf 'role\\tfile\\tgraph\\tlabel\\n'
      printf 'graph_gbz\\t%s\\tclip\\t%s\\n'        "${gbz}"     "${taxid}"
      printf 'graph_gfa\\t%s\\tclip\\t%s\\n'        "${gfa}"     "${taxid}"
      printf 'odgi\\t%s\\tclip\\t%s\\n'             "${og}"      "${taxid}"
      printf 'snarls\\t%s\\tclip\\t%s\\n'           "${snarls}"  "${taxid}"
      printf 'haplotype_index\\t%s\\tclip\\t%s\\n'  "${hapl}"    "${taxid}"
      printf 'variants_vcf\\t%s\\tclip\\t%s\\n'     "${vcf}"     "${taxid}"
      printf 'variants_vcf_tbi\\t%s\\tclip\\t%s\\n' "${vcf_tbi}" "${taxid}"
      printf 'reference_fasta\\t%s\\tclip\\t%s\\n'  "${ref_fa}"  "${taxid}"
      printf 'reference_fai\\t%s\\tclip\\t%s\\n'    "${ref_fai}" "${taxid}"
      if [ -n "${full_gbz}" ];    then printf 'graph_gbz\\t%s\\tfull\\t%s\\n' "${full_gbz}"    "${taxid}"; fi
      if [ -n "${full_gfa}" ];    then printf 'graph_gfa\\t%s\\tfull\\t%s\\n' "${full_gfa}"    "${taxid}"; fi
      if [ -n "${full_og}" ];     then printf 'odgi\\t%s\\tfull\\t%s\\n'      "${full_og}"     "${taxid}"; fi
      if [ -n "${full_snarls}" ]; then printf 'snarls\\t%s\\tfull\\t%s\\n'    "${full_snarls}" "${taxid}"; fi
    } > pangenome_manifest.tsv
    """

    stub:
    """
    printf 'role\\tfile\\tgraph\\tlabel\\n' > pangenome_manifest.tsv
    """
}
