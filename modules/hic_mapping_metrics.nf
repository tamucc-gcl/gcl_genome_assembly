/*
========================================================================================
    Hi-C mapping metrics (NO remapping)
========================================================================================
    - BAM metrics: mapping %, primary mapped, etc.
    - pairs.gz metrics: cis/trans, trans:cis
    - optional: scaffold-space cis/trans via AGP (contig->scaffold relabel)
========================================================================================
*/

process HIC_BAM_METRICS {
    tag "${haplotype_id}_${checkpoint}"
    label 'hic_qc'

    //publishDir "${params.outdir}/qc/hic/map/${checkpoint}", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), val(checkpoint), path(bam), path(bai)

    output:
    tuple val(haplotype_id), val(checkpoint), path("${haplotype_id}_${checkpoint}.bam_metrics.tsv"), emit: metrics

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    # Pull a few stable numbers from flagstat output (avoid parsing huge text later)
    # Notes:
    # - "mapped" includes secondary/supp; we also compute primary-only.
    # - For primary-only we exclude 0x100 (secondary) + 0x800 (supp) = 0x900 (2304)
    total=\$(samtools view -c ${bam})
    mapped=\$(samtools view -c -F 4 ${bam})
    primary_total=\$(samtools view -c -F 2304 ${bam})
    primary_mapped=\$(samtools view -c -F 2308 ${bam})  # 2304 + 4 = exclude secondary/supp + unmapped

    mapped_pct=\$(awk -v m="\$mapped" -v t="\$total" 'BEGIN{ if(t>0) printf("%.4f", 100*m/t); else print "NA"}')
    primary_mapped_pct=\$(awk -v m="\$primary_mapped" -v t="\$primary_total" 'BEGIN{ if(t>0) printf("%.4f", 100*m/t); else print "NA"}')

    {
        echo -e "haplotype_id\\tcheckpoint\\tbam_total_align\\tbam_mapped_align\\tbam_mapped_pct\\tbam_primary_align\\tbam_primary_mapped\\tbam_primary_mapped_pct"
        echo -e "${haplotype_id}\\t${checkpoint}\\t\$total\\t\$mapped\\t\$mapped_pct\\t\$primary_total\\t\$primary_mapped\\t\$primary_mapped_pct"
    } > ${haplotype_id}_${checkpoint}.bam_metrics.tsv
    """
}

process HIC_PAIRS_METRICS {
    tag "${haplotype_id}_${checkpoint}"
    label 'hic_qc'

    //publishDir "${params.outdir}/qc/hic/map/${checkpoint}", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id),
        val(checkpoint),
        path(pairs_gz),
        path(agp),
        path(parse_stats),
        path(dedup_stats)

    output:
    tuple val(haplotype_id), val(checkpoint), path("${haplotype_id}_${checkpoint}.pairs_metrics.tsv"), emit: metrics

    script:
    // Check if AGP is provided (not empty)
    def have_agp = agp.size() > 0
    """
    set -euo pipefail
    export LC_ALL=C

    # Build contig->scaffold map if AGP provided
    # AGP columns: object(1) ... component_id(6)
    HAVE_AGP=0
    if [[ "${have_agp}" == "true" && -s "${agp}" ]]; then
        awk 'BEGIN{FS="\\t"} !/^#/ && \$5 ~ /[NW]/ {print \$6"\\t"\$1}' "${agp}" > contig_to_scaffold.tsv
        HAVE_AGP=1
    else
        touch contig_to_scaffold.tsv
    fi

    # Count pairs, cis/trans in contig-space and (optionally) scaffold-space.
    # pairs format (pairtools): readID chr1 pos1 chr2 pos2 ...
    # chr1 = col2, chr2 = col4
    zcat ${pairs_gz} \\
        | awk -v FS="\\t" -v OFS="\\t" -v mapfile="contig_to_scaffold.tsv" -v have_agp="\$HAVE_AGP" '
        BEGIN {
            # load contig->scaffold mapping if present
            if (have_agp==1) {
                while ((getline < mapfile) > 0) { m[\$1]=\$2 }
                close(mapfile)
            }
        }
        /^#/ { next }
        {
            n++
            c1=\$2; c2=\$4

            # contig-space cis/trans
            if (c1==c2) cis_contig++; else trans_contig++

            # scaffold-space cis/trans (relabel only)
            if (have_agp==1) {
                s1 = (c1 in m ? m[c1] : c1)
                s2 = (c2 in m ? m[c2] : c2)
                if (s1==s2) cis_scaf++; else trans_scaf++
            }
        }
        END{
            # Avoid div-by-zero
            tc = (cis_contig>0 ? trans_contig/cis_contig : "NA")
            if (have_agp==1) ts = (cis_scaf>0 ? trans_scaf/cis_scaf : "NA"); else ts="NA"

            print "pairs_total", n
            print "cis_pairs_contig", cis_contig+0
            print "trans_pairs_contig", trans_contig+0
            print "trans_to_cis_contig", tc
            print "cis_pairs_scaffold", (have_agp==1 ? cis_scaf+0 : "NA")
            print "trans_pairs_scaffold", (have_agp==1 ? trans_scaf+0 : "NA")
            print "trans_to_cis_scaffold", ts
        }' > pairs_metrics.kv

    # Optional: compute "retention" using pairtools parse stats total pairs if available.
    # This is "UU kept" / "total input pairs seen by parse" (pair-level retention).
    total_in="NA"
    if [[ -s "${parse_stats}" ]]; then
        # pairtools --output-stats is typically a 2-column TSV: key <tab> value (or key value)
        total_in=\$(awk '
            BEGIN{FS="[\\t ]+"}
            tolower(\$1) ~ /(total.*pairs|total|total_pairs|total_read_pairs)/ {print \$2; exit}
        ' ${parse_stats} || true)
        [[ -z "\$total_in" ]] && total_in="NA"
    fi

    pairs_kept=\$(awk 'BEGIN{FS="\\t"} \$1=="pairs_total"{print \$2}' pairs_metrics.kv)
    retention="NA"
    if [[ "\$total_in" != "NA" && "\$total_in" -gt 0 ]]; then
        retention=\$(awk -v k="\$pairs_kept" -v t="\$total_in" 'BEGIN{printf("%.4f", 100*k/t)}')
    fi

    # Write tidy TSV (one row)
    cis_contig=\$(awk 'BEGIN{FS="\\t"} \$1=="cis_pairs_contig"{print \$2}' pairs_metrics.kv)
    trans_contig=\$(awk 'BEGIN{FS="\\t"} \$1=="trans_pairs_contig"{print \$2}' pairs_metrics.kv)
    tc_contig=\$(awk 'BEGIN{FS="\\t"} \$1=="trans_to_cis_contig"{print \$2}' pairs_metrics.kv)

    cis_scaf=\$(awk 'BEGIN{FS="\\t"} \$1=="cis_pairs_scaffold"{print \$2}' pairs_metrics.kv)
    trans_scaf=\$(awk 'BEGIN{FS="\\t"} \$1=="trans_pairs_scaffold"{print \$2}' pairs_metrics.kv)
    tc_scaf=\$(awk 'BEGIN{FS="\\t"} \$1=="trans_to_cis_scaffold"{print \$2}' pairs_metrics.kv)

    {
        echo -e "haplotype_id\\tcheckpoint\\tpairs_total\\tcis_pairs_contig\\ttrans_pairs_contig\\ttrans_to_cis_contig\\tcis_pairs_scaffold\\ttrans_pairs_scaffold\\ttrans_to_cis_scaffold\\tparse_total_pairs\\tretention_pct"
        echo -e "${haplotype_id}\\t${checkpoint}\\t\${pairs_kept}\\t\${cis_contig}\\t\${trans_contig}\\t\${tc_contig}\\t\${cis_scaf}\\t\${trans_scaf}\\t\${tc_scaf}\\t\${total_in}\\t\${retention}"
    } > ${haplotype_id}_${checkpoint}.pairs_metrics.tsv

    rm -f pairs_metrics.kv contig_to_scaffold.tsv
    """
}