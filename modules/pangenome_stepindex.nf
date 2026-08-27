/*
========================================================================================
    PANGENOME STEP INDEX MODULE
========================================================================================
    Repo location: modules/pangenome_stepindex.nf

    Builds an odgi step index for ONE per-chromosome graph. Exists for two reasons.

    1. IT IS A WORKAROUND FOR A BUG. odgi v0.9.2 -- the build inside the cactus v3.1.4
       image -- crashes when `odgi untangle` builds its own step index:

         sdsl::int_vector::operator[]: Assertion `idx < this->size()' failed
         during "[odgi::algorithms::stepindex] Collecting Steps"

       It is not an ID-space problem: `odgi sort -O` compaction does not fix it and
       `odgi validate` passes on both graph flavours. But standalone `odgi stepindex`
       succeeds on the same graph, and `odgi untangle -a <stpidx>` then returns 0. So the
       index is built here and handed to untangle explicitly.

       If the cactus image is ever bumped and the internal builder is fixed, this module
       becomes a pure optimisation rather than a necessity -- harmless either way.

    2. IT IS INVARIANT TO EVERY UNTANGLE PARAMETER. The index depends only on the graph,
       while -e / -j / -n / -m will be tuned repeatedly. Measured on chr1_1 (94.7 Mb, the
       largest chromosome, 4 threads): stepindex 3m42s / 6.9 GB, untangle 10m26s / 7.7 GB.
       The index is 26% of the cost and ~900 MB on disk, so caching it separately is the
       same argument as splitting PANGENOME_VARIANTS from PANGENOME_CLASSIFY.

    Input : tuple(taxid, flavor, og)
    Output: tuple(taxid, flavor, og, stpidx) / versions
========================================================================================
*/

process PANGENOME_STEPINDEX {
    tag "${taxid}:${flavor}:${og.simpleName}"
    label 'pangenome_stepindex'

    input:
    tuple val(taxid), val(flavor), path(og)

    output:
    tuple val(taxid), val(flavor), path(og), path("${og.simpleName}.stpidx"), emit: stpidx
    path("versions.tsv"),                                                    emit: versions

    script:
    def base = og.simpleName
    """
    set -euo pipefail
    export HOME="\$PWD"

    odgi stepindex -i ${og} -o ${base}.stpidx -t ${task.cpus} -P

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    : > ${og.simpleName}.stpidx
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
