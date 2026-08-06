/*
========================================================================================
    PANGENOME MANIFEST MODULE
========================================================================================
    Repo location: modules/pangenome_manifest.nf

    Write pangenome_manifest.tsv (role, file, label) enumerating the downstream-consumable
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
    tuple val(taxid), path(gbz), path(gfa), path(og), path(snarls), path(hapl),
          path(vcf), path(vcf_tbi), path(ref_fa), path(ref_fai)

    output:
    tuple val(taxid), path("pangenome_manifest.tsv"), emit: manifest

    script:
    """
    {
      printf 'role\\tfile\\tlabel\\n'
      printf 'graph_gbz\\t%s\\t%s\\n'        "${gbz}"     "${taxid}"
      printf 'graph_gfa\\t%s\\t%s\\n'        "${gfa}"     "${taxid}"
      printf 'odgi\\t%s\\t%s\\n'             "${og}"      "${taxid}"
      printf 'snarls\\t%s\\t%s\\n'           "${snarls}"  "${taxid}"
      printf 'haplotype_index\\t%s\\t%s\\n'  "${hapl}"    "${taxid}"
      printf 'variants_vcf\\t%s\\t%s\\n'     "${vcf}"     "${taxid}"
      printf 'variants_vcf_tbi\\t%s\\t%s\\n' "${vcf_tbi}" "${taxid}"
      printf 'reference_fasta\\t%s\\t%s\\n'  "${ref_fa}"  "${taxid}"
      printf 'reference_fai\\t%s\\t%s\\n'    "${ref_fai}" "${taxid}"
    } > pangenome_manifest.tsv
    """

    stub:
    """
    printf 'role\\tfile\\tlabel\\n' > pangenome_manifest.tsv
    """
}
