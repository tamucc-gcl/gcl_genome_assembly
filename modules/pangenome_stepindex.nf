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

    3. IT COMPACTS THE NODE ID SPACE, which is a second and unrelated odgi bug. Cactus
       compacts IDs for the FULL per-chromosome graphs but not the CLIP ones -- chr8_1.og
       holds 9,395,486 nodes with IDs to 142,488,970 (15.2x) while chr8_1.full.og is 1.00x.
       odgi's atomic bitvector is sized by node COUNT and indexed by node ID, so untangle
       aborts in "set target nodes" with `atomic_bitvector.hpp:159 Assertion idx < _size`.
       All 15 clip tasks failed this way; all 15 full tasks passed. `odgi sort -O` fixes it.

    Input : tuple(taxid, flavor, og)
    Output: tuple(taxid, flavor, COMPACTED og, stpidx) / versions
========================================================================================
*/

process PANGENOME_STEPINDEX {
    tag "${taxid}:${flavor}:${og.simpleName}"
    label 'pangenome_stepindex'

    input:
    tuple val(taxid), val(flavor), path(og)

    output:
    // the COMPACTED graph is emitted, not the input: the step index is built against
    // renumbered node IDs and is meaningless applied to the original.
    tuple val(taxid), val(flavor),
          path("${og.simpleName}.opt.og"),
          path("${og.simpleName}.opt.stpidx"), emit: stpidx
    path("versions.tsv"),                      emit: versions

    script:
    def base = og.simpleName
    """
    set -euo pipefail
    export HOME="\$PWD"

    # ---- compact the node ID space ---------------------------------------------------
    # MANDATORY for clip-flavour graphs. Cactus compacts IDs for the full per-chromosome
    # graphs but not the clip ones, so a clip graph can hold 9.4M nodes with IDs running to
    # 142M (measured on chr8_1: 15.2x). odgi's atomic bitvector is sized by node COUNT and
    # indexed by node ID, so untangle aborts in "set target nodes" with
    # `atomic_bitvector.hpp:159 Assertion idx < _size failed`. All 15 clip tasks failed this
    # way while all 15 full tasks succeeded.
    #
    # Applied unconditionally: detecting whether it is needed costs a full `odgi view -g`
    # scan for the max node id, which is more expensive than just doing it. Cheap and
    # idempotent on an already-dense graph.
    #
    # Renumbering is invisible downstream -- untangle reports in PATH coordinates, so no
    # node ID escapes this process.
    odgi sort -i ${og} -o ${base}.opt.og -O -t ${task.cpus} -P

    # index MUST be built on the compacted graph
    odgi stepindex -i ${base}.opt.og -o ${base}.opt.stpidx -t ${task.cpus} -P

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    : > ${og.simpleName}.opt.og
    : > ${og.simpleName}.opt.stpidx
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
