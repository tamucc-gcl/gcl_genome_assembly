/*
========================================================================================
    PANGENOME VARIANT CATALOG MODULE
========================================================================================
    Repo location: modules/pangenome_variants.nf

    Turn the raw decomposed VCF (vg deconstruct output, with LV/PS/AT tags) into a report
    variant catalog: keep top-level, position-anchored sites (vcfbub --max-level 0 + a
    reference-allele-length cap), REALIGN AND DECOMPOSE the surviving block alleles with
    vcfwave, split multiallelics, then classify every allele relative to the reference path
    by REF/ALT length, with an SV size spectrum and a type breakdown.

    vcfwave is not optional bookkeeping. `--max-level 0` keeps only top-level bubbles and
    discards nested ones, so without decomposition a bubble carrying internal variation is
    reported as one whole-block REF/ALT pair and the SNPs and indels inside it are lost. On
    Spratelloides that was 401,642 of 1,567,141 "SV" records (53% of the SV bp) with both
    sides >= 50 bp; sampled and self-aligned, they were >= 95% identical near-unique
    sequence differing by an embedded indel. vcfwave splits those into minimal variants.

    The pre-vcfwave block-level VCF is published as ${taxid}.blocks.vcf.gz so the effect of
    the decomposition stays auditable.

    vcfwave realigns with BiWFA and DOES sketch for inverted alignments (-I 64 / -K 17 by
    default), so unlike the vcfbub-only version some inversions are now resolved rather than
    left as substitutions. Decomposed records carry an ORIGIN INFO tag pointing back at the
    source record position, so any variant here can be traced to the block it came from.

    Two consequences worth knowing downstream: decomposed deletions produce haploid/missing
    genotypes at overlapping sites, and records above params.pangenome_vcfwave_args -L pass
    through undecomposed. What vcfwave cannot decompose lands in SV_COMPLEX (both sides
    >= minsv, similar length) or SV_BLOCKSUB (both sides >= minsv, lengths differ by
    >= minsv); a non-trivial SV_BLOCKSUB count is the signal that alleles are genuinely
    unrelated at that locus rather than an artifact of representation.

    Input : tuple(taxid, raw_vcf)          (CACTUS_PANGENOME.out.raw_vcf)
    Output: parents vcf (+tbi) / filtered vcf (+tbi) / blocks vcf (+tbi) / tier_audit /
            bcftools_stats / versions

    NOTE: variant_summary and sv_sizes have MOVED to PANGENOME_CLASSIFY, which classifies
    from graph traversals rather than allele string lengths and emits richer versions
    (topology labels, exclusive primary_class, AC/AN, and three distinct bp measures -- per
    allele, merged reference footprint, and pangenome node footprint).
========================================================================================
*/

process PANGENOME_VARIANTS {
    tag "${taxid}"
    label 'pangenome_variants'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(raw_vcf)

    output:
    tuple val(taxid), path("${taxid}.variants.vcf.gz"),     emit: vcf
    tuple val(taxid), path("${taxid}.variants.vcf.gz.tbi"), emit: tbi
    tuple val(taxid), path("${taxid}.blocks.vcf.gz"),       emit: blocks_vcf
    tuple val(taxid), path("${taxid}.blocks.vcf.gz.tbi"),   emit: blocks_tbi
    // parent tier: LV == 0, no vcfbub popping. The ONLY view where AT is interpretable,
    // because vcfwave + `bcftools norm -m` rewrite REF/ALT while AT is inherited from the
    // parent record, so allele i stops corresponding to traversal i+1.
    tuple val(taxid), path("${taxid}.parents.vcf.gz"),       emit: parents_vcf
    tuple val(taxid), path("${taxid}.parents.vcf.gz.tbi"),   emit: parents_tbi
    tuple val(taxid), path("${taxid}.tier_audit.tsv"),       emit: tier_audit
    tuple val(taxid), path("${taxid}.variants.bcftools_stats.txt"), emit: bcftools_stats
    path("versions.tsv"),                                   emit: versions

    script:
    def minsv = params.pangenome_sv_min_bp ?: 50
    def maxref = params.pangenome_vcfbub_max_ref ?: 100000
    def usewave = params.pangenome_vcfwave != false
    def wavemax = params.pangenome_vcfwave_max_len ?: 100000
    def nchunk = params.pangenome_vcfwave_chunks ?: (task.cpus * 4)
    def waveargs = params.pangenome_vcfwave_args ?: ''
    def wavecmd = usewave ? "vcfwave -L ${wavemax} ${waveargs}" : 'cat'
    """
    set -euo pipefail

    # top-level, position-anchored sites (drop nested + oversized bubbles). Published as-is
    # so the block-level view survives and the decomposition below can be diffed against it.
    vcfbub --input ${raw_vcf} --max-level 0 --max-ref-length ${maxref} \\
        | bgzip > ${taxid}.blocks.vcf.gz
    tabix -p vcf ${taxid}.blocks.vcf.gz

    WAVE_CMD='${wavecmd}'
    export WAVE_CMD

    # WAVE_CMD is built in the Groovy prelude, where a backslash-dollar escapes Groovy
    # interpolation rather than protecting a shell dollar. Getting that wrong leaves an
    # un-interpolated placeholder that only fails later inside eval, once per chunk.
    case "\$WAVE_CMD" in
      *'\${'*) echo "ERROR: WAVE_CMD holds an unexpanded Groovy placeholder: \$WAVE_CMD" >&2
                echo '       the script: prelude must interpolate, not backslash-escape.' >&2
                exit 1 ;;
    esac

    # vcflib < ~1.0.9 has no vcfwave, and the env would otherwise build fine and then break
    # mid-pipe hours in. Check up front, with a message that says what is actually wrong.
    if [ "${usewave}" = "true" ]; then
      command -v vcfwave >/dev/null 2>&1 || {
        echo "ERROR: vcfwave not found. vcflib is too old (>=1.0.9 required; 1.0.1 ships" >&2
        echo "       only vcfallelicprimitives). Pin bioconda::vcflib in the" >&2
        echo "       pangenome_variants label, or set params.pangenome_vcfwave = false." >&2
        exit 1; }
    fi

    # vcfwave is single-threaded by design (its own help says -t usually makes it slower),
    # so parallelism has to come from chunking. Decomposition is record-LOCAL: splitting a
    # 23,215-record slice in half and concatenating gave byte-identical sorted output, so
    # chunk boundaries cannot change the result and chunks need not respect contigs.
    # Equal-RECORD chunks with more chunks than cores lets xargs schedule dynamically, which
    # matters because BiWFA cost tracks allele LENGTH, not record count -- a few long alleles
    # in one chunk would otherwise stall the whole task.
    mkdir -p wchunks
    bcftools view -h ${taxid}.blocks.vcf.gz > wchunks/hdr.vcf
    bcftools view -H ${taxid}.blocks.vcf.gz > wchunks/body.vcf
    NREC=\$(wc -l < wchunks/body.vcf)
    NCHUNK=${nchunk}
    [ "\$NREC" -lt "\$NCHUNK" ] && NCHUNK=1
    echo "decomposing \$NREC block records in \$NCHUNK chunks across ${task.cpus} cores" >&2

    if [ "\$NREC" -eq 0 ]; then
      cp wchunks/hdr.vcf wchunks/all.vcf
    else
      split -d -a 5 -n l/\$NCHUNK wchunks/body.vcf wchunks/part_
      for f in wchunks/part_*; do cat wchunks/hdr.vcf "\$f" > "\$f.vcf"; rm -f "\$f"; done

      # eval so WAVE_CMD's flags word-split correctly; one file per invocation
      printf '%s\n' '#!/bin/sh' 'set -eu' 'eval "\$WAVE_CMD \"\$1\"" > "\$1.wave"' \
          > wchunks/run_one.sh
      chmod +x wchunks/run_one.sh
      ls wchunks/part_*.vcf | xargs -P ${task.cpus} -n1 wchunks/run_one.sh

      # awk not grep: grep -v exits 1 on a header-only chunk and pipefail would kill the task.
      # Header comes from a vcfwave output so its added INFO tags (ORIGIN) are declared.
      HDR=\$(ls wchunks/part_*.vcf.wave | head -1)
      { awk '/^#/' "\$HDR"; awk '!/^#/' wchunks/part_*.vcf.wave; } > wchunks/all.vcf
    fi

    # split multiallelics, THEN sort: vcfwave emits records out of order after splitting a
    # block, and tabix needs them sorted. bcftools sort takes a file, not stdin.
    bcftools norm -m -any -Oz -o wave.norm.vcf.gz wchunks/all.vcf
    bcftools sort -T ./bcfsort -Oz -o ${taxid}.variants.vcf.gz wave.norm.vcf.gz
    rm -rf wchunks wave.norm.vcf.gz
    tabix -p vcf ${taxid}.variants.vcf.gz

    # ---- parent tier -----------------------------------------------------------------
    # LV == 0 is exactly the top-level bubble set: on Spratelloides this yields 23,877,797
    # records, identical to `vcfbub --max-level 0`. Derived by filter rather than a second
    # vcfbub run, and NOT OR-ed with a size threshold -- vcfbub's hierarchy already handles
    # nesting, so OR-ing a span floor would count a large nested bubble twice, once inside
    # its parent and once as itself.
    bcftools view -i 'INFO/LV=0' -Oz -o ${taxid}.parents.vcf.gz ${raw_vcf}
    tabix -p vcf ${taxid}.parents.vcf.gz

    # ---- tier audit ------------------------------------------------------------------
    # Emitted every run. The vcfbub cap is resolution recovery, not a ceiling: it ADDS
    # records by popping oversized top-level bubbles into their children (31,365,160 with
    # the cap vs 23,877,797 without on this graph, and 0 records at --max-ref-length 0,
    # which is a literal zero limit). Reporting the three counts side by side keeps that
    # from being rediscovered, and makes any change in the representation visible.
    {
      printf 'tier\\trecords\\tnote\\n'
      printf 'raw\\t%s\\tall nesting levels, pre-vcfbub\\n' \\
        "\$(bcftools index -n ${raw_vcf} 2>/dev/null || bcftools view -H ${raw_vcf} | wc -l)"
      printf 'parent_LV0\\t%s\\ttop-level only; AT interpretable here\\n' \\
        "\$(bcftools index -n ${taxid}.parents.vcf.gz)"
      printf 'blocks_vcfbub\\t%s\\tmax-level 0, max-ref-length ${maxref}\\n' \\
        "\$(bcftools index -n ${taxid}.blocks.vcf.gz)"
      printf 'fine_decomposed\\t%s\\tafter vcfwave + norm -m; AT NO LONGER VALID\\n' \\
        "\$(bcftools index -n ${taxid}.variants.vcf.gz)"
    } > ${taxid}.tier_audit.tsv

    # Classification is NOT done here. PANGENOME_CLASSIFY derives classes from AT traversals
    # rather than REF/ALT string lengths; the old awk discarded AT, LV/PS, AC/AN and every
    # sample column before deciding, and its SV_COMPLEX / SV_BLOCKSUB were not variant
    # classes but the buckets it used when it could not tell.

    # bcftools stats on the filtered catalog -> MultiQC bcftools module
    bcftools stats ${taxid}.variants.vcf.gz > ${taxid}.variants.bcftools_stats.txt

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tbcftools\\t%s\\n' "${task.process}" "\$(bcftools --version 2>&1 | awk 'NR==1{print \$2}')"
      printf '%s\\tvcfbub\\t%s\\n'   "${task.process}" "\$(vcfbub --help 2>&1 | awk 'NR==1{print \$NF}')"
      printf '%s\\tvcfwave\\t%s\\n'  "${task.process}" "\$(vcfwave -h 2>&1 | awk '/version/{print \$NF; exit}' || echo NA)"
    } > versions.tsv
    """

    stub:
    """
    : > ${taxid}.variants.vcf.gz
    : > ${taxid}.variants.vcf.gz.tbi
    : > ${taxid}.parents.vcf.gz
    : > ${taxid}.parents.vcf.gz.tbi
    : > ${taxid}.blocks.vcf.gz
    : > ${taxid}.blocks.vcf.gz.tbi
    printf 'tier\\trecords\\tnote\\n' > ${taxid}.tier_audit.tsv
    : > ${taxid}.variants.bcftools_stats.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
