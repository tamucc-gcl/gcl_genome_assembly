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
    Output: filtered vcf (+tbi) / blocks vcf (+tbi) / variant_summary / sv_sizes / versions
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
    tuple val(taxid), path("${taxid}.variant_summary.tsv"), emit: summary
    tuple val(taxid), path("${taxid}.sv_sizes.tsv"),        emit: sv_sizes
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

    bcftools query -f '%REF\\t%ALT\\n' ${taxid}.variants.vcf.gz > refalt.tsv

    # SNP + INDEL + SV is the exhaustive top-level partition (the report sums those three);
    # SV_INS + SV_DEL + SV_COMPLEX + SV_BLOCKSUB sums to SV. m = the SHORTER allele: a clean
    # presence/absence event has m ~ 1, whereas m >= minsv means both haplotypes carry
    # substantial sequence at the locus and it was never a simple insertion or deletion.
    awk -v minsv=${minsv} 'BEGIN{FS=OFS="\\t"}
        { rl=length(\$1); al=length(\$2); d=(rl>al?rl-al:al-rl); m=(rl<al?rl:al)
          if(rl==1 && al==1)              snp++
          else if(d<minsv && m<minsv)     indel++
          else if(d<minsv)              { sv++; cpx++ }
          else if(m>=minsv)             { sv++; blk++ }
          else { sv++; if(al>rl) ins++; else del++ } }
        END{ print "class","count"
             print "SNP",snp+0; print "INDEL",indel+0; print "SV",sv+0
             print "SV_INS",ins+0; print "SV_DEL",del+0
             print "SV_COMPLEX",cpx+0; print "SV_BLOCKSUB",blk+0 }' \\
        refalt.tsv > ${taxid}.variant_summary.tsv

    # SV size spectrum (bp) + type, for the histogram. ref_len/alt_len are carried so the
    # shorter allele is recoverable: sv_size_bp is only |REF|-|ALT|, which for a block
    # substitution silently discards the kilobases sitting on the other side.
    awk -v minsv=${minsv} 'BEGIN{FS=OFS="\\t"; print "sv_size_bp","sv_type","ref_len","alt_len"}
        { rl=length(\$1); al=length(\$2); d=(rl>al?rl-al:al-rl); m=(rl<al?rl:al)
          if(rl==1 && al==1) next
          if(d<minsv && m>=minsv)      { print d,"COMPLEX",rl,al }
          else if(d>=minsv && m>=minsv){ print d,"BLOCKSUB",rl,al }
          else if(d>=minsv)            { print d,(al>rl?"INS":"DEL"),rl,al } }' \\
        refalt.tsv > ${taxid}.sv_sizes.tsv

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
    printf 'class\\tcount\\n' > ${taxid}.variant_summary.tsv
    printf 'sv_size_bp\\tsv_type\\n' > ${taxid}.sv_sizes.tsv
    : > ${taxid}.variants.bcftools_stats.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
